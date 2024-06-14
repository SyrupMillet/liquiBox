!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,nx,Lx
   use fft3d_class,       only: fft3d
   use ddadi_class,       only: ddadi
   use multivdscalar_class, only: multivdscalar
   use lpt_class,         only: lpt
   use lowmach_class,     only: lowmach
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single-phase incompressible flow solver, pressure and implicit solvers, and a time tracker
   type(fft3d),       public :: ps
   type(lowmach),      public :: fs
   type(timetracker), public :: time
   type(lpt),         public :: lp
   type(multivdscalar), public :: msc  !< scalar (Composition) solver
   type(ddadi),         public :: mss

   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,lptfile,hitfile,cvgfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi, resRHO
   real(WP), dimension(:,:,:), allocatable :: srcUlp,srcVlp,srcWlp
   real(WP), dimension(:,:,:,:), allocatable :: SR
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   real(WP), dimension(:,:,:), allocatable :: rho0          !> density field

   !> Multiscalar solver work arrays
   real(WP), dimension(:,:,:,:), allocatable :: resMSC,srcMSC,mscTmp    !< Residuals, source, temp solution for multiscalar solver
   logical, dimension(:,:,:,:), allocatable :: bqflag                   !< Flag for bquick scheme

   !> gas density array for computing mean density
   real(WP), dimension(:), allocatable :: speciesFluidDensities
   integer:: speciesNum

   !> Fluid, forcing, and particle parameters
   real(WP) :: visc,meanU,meanV,meanW
   real(WP) :: Urms0,TKE0,EPS0,Re_max
   real(WP) :: TKE,URMS,EPSp
   real(WP) :: tauinf,G,Gdtau,Gdtaui,dx
   
   !> For monitoring
   real(WP) :: EPS
   real(WP) :: Re_L,Re_lambda
   real(WP) :: eta,ell
   real(WP) :: dx_eta,ell_Lx,Re_ratio,eps_ratio,tke_ratio,nondtime
   

contains
   
   
   !> Compute turbulence stats
   subroutine compute_stats()
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      real(WP) :: myTKE,myEPS,rhoMax
      integer :: i,j,k,ierr
      
      ! Compute mean velocities
      call fs%cfg%integrate(A=Ui,integral=meanU); meanU=meanU/fs%cfg%vol_total
      call fs%cfg%integrate(A=Vi,integral=meanV); meanV=meanV/fs%cfg%vol_total
      call fs%cfg%integrate(A=Wi,integral=meanW); meanW=meanW/fs%cfg%vol_total
      
      ! Compute strainrate and grad(u)
      call fs%get_strainrate(SR=SR)
      call fs%get_gradu(gradu=gradU)
      
      ! Compute current TKE and dissipation rate
      myTKE=0.0_WP
      myEPS=0.0_WP
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               myTKE=myTKE+0.5_WP*((Ui(i,j,k)-meanU)**2+(Vi(i,j,k)-meanV)**2+(Wi(i,j,k)-meanW)**2)*fs%cfg%vol(i,j,k)
               myEPS=myEPS+2.0_WP*fs%visc(i,j,k)*fs%cfg%vol(i,j,k)*(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))/fs%rho(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(myTKE,TKE,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); TKE=TKE/fs%cfg%vol_total
      call MPI_ALLREDUCE(myEPS,EPS,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); EPS=EPS/fs%cfg%vol_total
      
      call msc%get_max()
      rhoMax = msc%rhomax

      ! Compute standard parameters for HIT
      Urms=sqrt(2.0_WP/3.0_WP*TKE)
      Re_L=TKE**2.0_WP/EPS/(visc/rhoMax)
      Re_lambda=sqrt(20.0_WP*Re_L/3.0_WP)
      eta=((visc/rhoMax)**3.0_WP/EPS)**0.25_WP
      ell=(0.6667_WP*TKE)**1.5_WP/EPS
      
      ! Some more useful info
      nondtime =time%t/tauinf
      dx_eta   =dx/eta
      eps_ratio=EPS/EPS0
      tke_ratio=TKE/TKE0
      ell_Lx   =ell/Lx
      Re_ratio =Re_lambda/Re_max
      
   end subroutine compute_stats


   !> Compute Density and Turbulent Diffusivity for MSC solver
   subroutine compute_msc_density()
      integer :: i,j,k
      real(WP) :: turDiff
      do k=msc%cfg%kmino_,msc%cfg%kmaxo_
         do j=msc%cfg%jmino_,msc%cfg%jmaxo_
            do i=msc%cfg%imino_,msc%cfg%imaxo_
               rho0(i,j,k)=dot_product(msc%SC(i,j,k,:),speciesFluidDensities)
               msc%rho(i,j,k)=dot_product(msc%SC(i,j,k,:),speciesFluidDensities)*(1-lp%VF(i,j,k))
            end do
         end do
      end do

   end subroutine compute_msc_density

   subroutine compute_msc_diffusivity()
      real(WP) :: turDiff

      ! Compute turbulent diffusivity
      turDiff=(0.09_WP/0.7_WP)*(TKE**2.0_WP)/EPS
      msc%diff=turDiff

   end subroutine compute_msc_diffusivity
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! read species properties
      call param_read('Species Numbers',speciesNum)
      allocate(speciesFluidDensities(speciesNum))
      call param_read('Species gas phase densities',speciesFluidDensities)
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resRHO(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcUlp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcVlp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcWlp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         allocate(SR  (1:6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      
         allocate(rho0(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      
         !< Scalar Solver
         allocate(resMSC(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,speciesNum))
         allocate(srcMSC(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,speciesNum)) ; srcMSC=0.0_WP
         allocate(mscTmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,speciesNum))
         allocate(bqflag(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,speciesNum))
      end block allocate_work_arrays
       
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         ! Create flow solver
         fs=lowmach(cfg=cfg,name='NS solver')
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Prepare and configure pressure solver
         ps=fft3d(cfg=cfg,name='Pressure',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)
      end block create_and_initialize_flow_solver

      ! Create a multiscalar solver
      create_multiscalar: block
         use multivdscalar_class, only: dirichlet,neumann,quick,bquick,upwind
         real(WP) :: diffusivity
         ! Create scalar solver
         msc=multivdscalar(cfg=cfg,scheme=bquick,nscalar=speciesNum,name='Variable density multi scalar')
         call param_read('Species names',msc%SCname)

         ! Configure implicit scalar solver
         mss=ddadi(cfg=cfg,name='multiScalar',nst=13)
         ! Setup the solver
         call msc%setup(implicit_solver=mss)

      end block create_multiscalar
      
      ! Initialize LPT solver
      initialize_lpt: block
         use random, only: random_uniform
         real(WP) :: dp, px, py, pz, pp
         integer :: i,np, np_axis
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get drag model from the inpit
         call param_read('Drag model',lp%drag_model,default='Schiller-Naumann')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Get particle diameter from the input
         call param_read('Particle diameter',dp)
         ! Get number of particles
         call param_read('Number of particles per axis',np)
         np_axis = np
         pp = Lx/real(np,WP)
         px = pp/2.0_WP ; py = pp/2.0_WP ; pz = pp/2.0_WP
         np = np**3
         ! Root process initializes np particles randomly
         if (lp%cfg%amRoot) then
            call lp%resize(np)
            do i=1,np
               ! Give id
               lp%p(i)%id=int(i,8)
               ! Set the diameter
               lp%p(i)%d=dp
               ! assign position in the domain
               lp%p(i)%pos= [px+mod(i-1,np_axis)*pp,py+mod((i-1)/np_axis,np_axis)*pp,pz+mod((i-1)/(np_axis*np_axis),np_axis)*pp]
               ! Give zero velocity
               lp%p(i)%vel=0.0_WP
               lp%p(i)%angVel=0.0_WP
               ! Zero out collision forces
               lp%p(i)%Acol=0.0_WP
               lp%p(i)%Tcol=0.0_WP
               ! Give zero dt
               lp%p(i)%dt=0.0_WP
               ! Locate the particle on the mesh
               lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
               ! Activate the particle
               lp%p(i)%flag=0
            end do
         end if
         ! Distribute particles
         call lp%sync()
         ! Get initial particle volume fraction
         call lp%update_VF()
         ! Set collision timescale
         call param_read('Collision timescale',lp%tau_col,default=15.0_WP*time%dt)
         ! Set coefficient of restitution
         call param_read('Coefficient of restitution',lp%e_n)
         call param_read('Wall restitution',lp%e_w,default=lp%e_n)
         call param_read('Friction coefficient',lp%mu_f,default=0.0_WP)
         ! Set gravity
         call param_read('Gravity',lp%gravity)
         if (lp%cfg%amRoot) then
            print*,"===== Particle Setup Description ====="
            print*,'Number of particles', np
         end if

      end block initialize_lpt

      ! Initialize our multi mixture fraction field
      initialize_multiscalar: block
         use multivdscalar_class, only: bcond
         integer :: n,i,j,k,isc
         type(bcond), pointer :: mybc
         real(WP), dimension(speciesNum) :: initComp
         call param_read('Gas phase inital composition',initComp)
         ! set initial field
         do k=msc%cfg%kmino_, msc%cfg%kmaxo_
            do j=msc%cfg%jmino_, msc%cfg%jmaxo_
               do i=msc%cfg%imino_, msc%cfg%imaxo_
                  msc%SC(i,j,k,:)=initComp
               end do
            end do
         end do

         ! Compute density
         call compute_msc_density()
         call msc%get_max()
         call msc%get_int()

         call msc%cfg%sync(msc%rho)
         call msc%cfg%sync(fs%visc)
         do isc=1,msc%nscalar
            call msc%cfg%sync(msc%SC  (:,:,:,isc))
            call msc%cfg%sync(msc%diff(:,:,:,isc))
         end do

      end block initialize_multiscalar
      
      ! Prepare initial velocity field
      initialize_velocity: block
         use random,    only: random_normal
         use param,     only: param_exists
         use mathtools, only: Pi
         use messager,  only: die
         use, intrinsic :: iso_fortran_env, only: output_unit
         integer :: i,j,k
         real(WP) :: myKE,taueta, rhoMax
         ! Read in forcing, grid, and initial velocity field parameters
         call param_read('Forcing constant',G)
         dx=Lx/real(nx,WP)

         fs%rho = msc%rho
         rhoMax = msc%rhomax

         if (param_exists('Steady-state TKE')) then
            ! A target TKE was provided
            call param_read('Steady-state TKE',TKE0)
            EPS0=5.0_WP*(0.6667_WP*TKE0)**1.5_WP/Lx
         else
            ! Force to maximum Re_lambda
            EPS0=(visc/rhoMax)**3*(Pi*cfg%nx/(1.5_WP*Lx))**4
            TKE0=1.5_WP*(0.2_WP*Lx*EPS0)**(0.6667_WP)
         end if
         Re_max=sqrt(15.0_WP*sqrt(0.6667_WP*TKE0)*0.2_WP*Lx/(visc/rhoMax))
         tauinf=2.0_WP*TKE0/(3.0_WP*EPS0)
         taueta=sqrt((visc/rhoMax)/EPS0)
         Gdtau =G/tauinf
         Gdtaui=1.0_WP/Gdtau
         if (Gdtaui.lt.time%dt) call die('[linear_forcing] Controller time constant less than timestep')
         Urms0=sqrt(0.6667_WP*TKE0)
         ! Print out some expected turbulence properties
         if (fs%cfg%amRoot) then
            write(output_unit,'("Expected turbulence properties:")')
            write(output_unit,'("Re_lambda = ",es12.5)') Re_max
            write(output_unit,'("tau_eddy  = ",es12.5)') tauinf
            write(output_unit,'("Urms      = ",es12.5)') Urms0
            write(output_unit,'("tau_eta   = ",es12.5)') taueta
         end if

         ! Gaussian initial field
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  fs%U(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
                  fs%V(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
                  fs%W(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
               end do
            end do
         end do
         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         call fs%cfg%sync(fs%W)
         ! Compute mean and remove it from the velocity field to obtain <U>=0
         call fs%cfg%integrate(A=fs%U,integral=meanU); meanU=meanU/fs%cfg%vol_total
         call fs%cfg%integrate(A=fs%V,integral=meanV); meanV=meanV/fs%cfg%vol_total
         call fs%cfg%integrate(A=fs%W,integral=meanW); meanW=meanW/fs%cfg%vol_total
         fs%U=fs%U-meanU
         fs%V=fs%V-meanV
         fs%W=fs%W-meanW

         ! Initialize density field residual to zero
         resRHO=0.0_WP
         ! Project to ensure divergence-free
         call fs%get_div(drhodt=resRHO)
         fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
         fs%psolv%sol=0.0_WP
         call fs%psolv%solve()
         call fs%shift_p(fs%psolv%sol)
         call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
         fs%P=fs%P+fs%psolv%sol
         fs%U=fs%U-time%dt*resU/fs%rho
         fs%V=fs%V-time%dt*resV/fs%rho
         fs%W=fs%W-time%dt*resW/fs%rho

         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(drhodt=resRHO)

         ! Compute turbulence stats
         call compute_stats()
      end block initialize_velocity

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         pmesh=partmesh(nvar=0,nvec=0,name='lpt')
         call lp%update_partmesh(pmesh)
      end block create_pmesh

      
      ! Add Ensight output
      create_ensight: block
         integer :: i
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='HIT')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_particle('particles',pmesh)
         call ens_out%add_scalar('particleVF',lp%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('density',rho0)
         do i=1,speciesNum
            call ens_out%add_scalar("frac_"//trim(msc%SCname(i)),msc%SC(:,:,:,i))
         end do
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create LPT monitor
         call lp%get_max()
         lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
         call lptfile%add_column(time%n,'Timestep number')
         call lptfile%add_column(time%t,'Time')
         call lptfile%add_column(lp%np,'Particle number')
         call lptfile%add_column(lp%Umin,'Particle Umin')
         call lptfile%add_column(lp%Umax,'Particle Umax')
         call lptfile%add_column(lp%Vmin,'Particle Vmin')
         call lptfile%add_column(lp%Vmax,'Particle Vmax')
         call lptfile%add_column(lp%Wmin,'Particle Wmin')
         call lptfile%add_column(lp%Wmax,'Particle Wmax')
         call lptfile%add_column(lp%dmin,'Particle dmin')
         call lptfile%add_column(lp%dmax,'Particle dmax')
         call lptfile%write()
         ! Create hit monitor
         hitfile=monitor(fs%cfg%amRoot,'hit')
         call hitfile%add_column(time%n,'Timestep number')
         call hitfile%add_column(time%t,'Time')
         call hitfile%add_column(Re_L,'Re_L')
         call hitfile%add_column(Re_lambda,'Re_lambda')
         call hitfile%add_column(eta,'eta')
         call hitfile%add_column(TKE,'TKE')
         call hitfile%add_column(URMS,'Urms')
         call hitfile%add_column(EPS,'EPS')
         call hitfile%add_column(ell,'L')
         call hitfile%write()
         ! Create hit convergence monitor
         cvgfile=monitor(fs%cfg%amRoot,'convergence')
         call cvgfile%add_column(time%n,'Timestep number')
         call cvgfile%add_column(time%t,'Time')
         call cvgfile%add_column(nondtime,'Time/t_int')
         call cvgfile%add_column(Re_ratio,'Re_ratio')
         call cvgfile%add_column(eps_ratio,'EPS_ratio')
         call cvgfile%add_column(tke_ratio,'TKE_ratio')
         call cvgfile%add_column(dx_eta,'dx/eta')
         call cvgfile%add_column(ell_Lx,'ell/Lx')
         call cvgfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Time integrate our problem
   subroutine simulation_run
      implicit none
      real(WP) :: cfl
      integer :: isc,i,j,k
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
         call fs%get_cfl(time%dt,cfl)
         time%cfl=max(time%cfl,cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old velocity
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW

         ! Remember old density
         msc%rhoold=msc%rho
         msc%SCold =msc%SC

         call fs%get_div_stress(resU,resV,resW)

         ! =================== Particle Solver ===================
         call lp%collide(dt=time%dtmid)
         call lp%advance(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W,rho=fs%rho,visc=fs%visc,stress_x=resU,stress_y=resV,stress_z=resW,&
             srcU=srcUlp,srcV=srcVlp,srcW=srcWlp)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! ============= MULTI SCALAR SOLVER =======================
            call compute_msc_diffusivity()
            
            !< Reset metrics for bquick
            call msc%metric_reset()

            ! Build mid-time scalar
            msc%SC=0.5_WP*(msc%SC+msc%SCold)

            ! Get rhoSC time derivative
            call msc%get_drhoSCdt(resMSC,fs%rhoU,fs%rhoV,fs%rhoW)

            ! Assemble explicit residual
            do isc=1,msc%nscalar
               resMSC(:,:,:,isc)=time%dt*resMSC(:,:,:,isc)&
               &   -2.0_WP*msc%rho*msc%SC(:,:,:,isc)+(msc%rho+msc%rhoold)*msc%SCold(:,:,:,isc)&
               &   +time%dt*srcMSC(:,:,:,isc)  !< Source term
               !< Get temperary solution for bquick
               mscTmp(:,:,:,isc)=2.0_WP*msc%SC(:,:,:,isc)-msc%SCold(:,:,:,isc)+resMSC(:,:,:,isc)/msc%rho
            end do

            !< For each species, if <0 or >1, set flag
            bqflag = .false.
            do isc=1,msc%nscalar
               do k=msc%cfg%kmino_,msc%cfg%kmaxo_
                  do j=msc%cfg%jmino_,msc%cfg%jmaxo_
                     do i=msc%cfg%imino_,msc%cfg%imaxo_
                        if ((mscTmp(i,j,k,isc).le.0.0_WP).or.(mscTmp(i,j,k,isc).ge.1.0_WP)) then
                           bqflag(i,j,k,isc)=.true.
                        end if
                     end do
                  end do
               end do
            end do

            ! Adjust metrics
            call msc%metric_adjust(mscTmp,bqflag)

            ! re-Assemble explicit residual
            call msc%get_drhoSCdt(resMSC,fs%rhoU,fs%rhoV,fs%rhoW)
            do isc=1,msc%nscalar
               resMSC(:,:,:,isc) = time%dt*resMSC(:,:,:,isc)&
               &   -2.0_WP*msc%rho*msc%SC(:,:,:,isc)+(msc%rho+msc%rhoold)*msc%SCold(:,:,:,isc)&
               &   +time%dt*srcMSC(:,:,:,isc)  !< Source term
            end do

            call msc%solve_implicit(time%dt,resMSC,fs%rhoU,fs%rhoV,fs%rhoW)

            ! Apply these residuals, reconstruct the n+1 field
            msc%SC=2.0_WP*msc%SC-msc%SCold+resMSC

            call compute_msc_density()
            call msc%get_max()
            call msc%get_int()

            ! =================== FLOW SOLVER =======================

            fs%rho=0.5_WP*(msc%rho+msc%rhoold)
            
            ! Build mid-time velocity and momentum
            fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW
            
            ! Add linear forcing term based on Bassenne et al. (2016)
            linear_forcing: block
               use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
               use parallel, only: MPI_REAL_WP
               real(WP) :: myTKE,A,myEPSp,rhoMax
               integer :: i,j,k,ierr

               call msc%get_max()
               rhoMax = msc%rhomax
               ! Calculate mean velocity
               call fs%cfg%integrate(A=fs%U,integral=meanU); meanU=meanU/fs%cfg%vol_total
               call fs%cfg%integrate(A=fs%V,integral=meanV); meanV=meanV/fs%cfg%vol_total
               call fs%cfg%integrate(A=fs%W,integral=meanW); meanW=meanW/fs%cfg%vol_total
               ! Calculate TKE and pseudo-EPS
               call fs%interp_vel(Ui,Vi,Wi)
               call fs%get_gradu(gradu=gradU)
               myTKE=0.0_WP; myEPSp=0.0_WP
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        myTKE =myTKE + 0.5_WP*((Ui(i,j,k)-meanU)**2+(Vi(i,j,k)-meanV)**2+(Wi(i,j,k)-meanW)**2)*fs%cfg%vol(i,j,k)
                        myEPSp=myEPSp+fs%cfg%vol(i,j,k)*fs%visc(i,j,k)*(gradU(1,1,i,j,k)**2+gradU(1,2,i,j,k)**2+gradU(1,3,i,j,k)**2+&
                        &                                               gradU(2,1,i,j,k)**2+gradU(2,2,i,j,k)**2+gradU(2,3,i,j,k)**2+&
                        &                                               gradU(3,1,i,j,k)**2+gradU(3,2,i,j,k)**2+gradU(3,3,i,j,k)**2)&
                        &                                               /fs%rho(i,j,k)
                     end do
                  end do
               end do
               call MPI_ALLREDUCE(myTKE ,TKE ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); TKE =TKE /fs%cfg%vol_total
               call MPI_ALLREDUCE(myEPSp,EPSp,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); EPSp=EPSp/fs%cfg%vol_total
               A=(EPSp-Gdtau*(TKE-TKE0))/(2.0_WP*TKE)*rhoMax
               resU=resU+time%dt*(fs%U-meanU)*A
               resV=resV+time%dt*(fs%V-meanV)*A
               resW=resW+time%dt*(fs%W-meanW)*A
            end block linear_forcing

            ! Add momentum source term from lpt
            add_lpt_src: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        resU(i,j,k)=resU(i,j,k)+sum(fs%itpr_x(:,i,j,k)*srcUlp(i-1:i,j,k))
                        resV(i,j,k)=resV(i,j,k)+sum(fs%itpr_y(:,i,j,k)*srcVlp(i,j-1:j,k))
                        resW(i,j,k)=resW(i,j,k)+sum(fs%itpr_z(:,i,j,k)*srcWlp(i,j,k-1:k))
                     end do
                  end do
               end do
            end block add_lpt_src

            ! TODO
            ! Form implicit residuals ?
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU/fs%rho
            fs%V=2.0_WP*fs%V-fs%Vold+resV/fs%rho
            fs%W=2.0_WP*fs%W-fs%Wold+resW/fs%rho
            
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            call fs%rho_multiply()
            call fs%apply_bcond(time%tmid,time%dtmid)
            
            ! Solve Poisson equation
            call msc%get_drhodt(dt=time%dt,drhodt=resRHO)
            call fs%correct_mfr(drhodt=resRHO)
            call fs%get_div(drhodt=resRHO)

            fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct momentum and rebuild velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dtmid*resU
            fs%rhoV=fs%rhoV-time%dtmid*resV
            fs%rhoW=fs%rhoW-time%dtmid*resW
            
            call fs%rho_divide()

            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call msc%get_drhodt(dt=time%dt,drhodt=resRHO)
         call fs%get_div(drhodt=resRHO)
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call lp%update_partmesh(pmesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call compute_stats()
         call fs%get_max()
         call mfile%write()
         call cflfile%write()
         call lp%get_max()
         call lptfile%write()
         call hitfile%write()
         call cvgfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,SR,gradU)
      
   end subroutine simulation_final
   
   
   
   
   
end module simulation
