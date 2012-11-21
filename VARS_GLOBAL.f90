!#####################################################################
!     Program  : VARS_GLOBAL
!     TYPE     : Module
!     PURPOSE  : Defines the global variables used thru all the code
!     AUTHORS  : Adriano Amaricci
!#####################################################################
MODULE VARS_GLOBAL
  !Local:
  USE CONTOUR_GF
  !SciFor library
  USE PARSE_CMD
  USE COMMON_VARS
  USE GREENFUNX
  USE TIMER
  USE VECTORS
  USE SQUARE_LATTICE
  USE IOTOOLS
  USE FFTGF
  USE SPLINE
  USE TOOLS
  USE MPI
  implicit none

  !Version revision
  include "revision.inc"

  !Gloabl  variables
  !=========================================================
  integer           :: L             !a big number
  integer           :: Ltau          !Imaginary time slices
  integer           :: Nstep         !Number of Time steps
  integer           :: Lk            !total lattice  dimension
  integer           :: Lkreduced     !reduced lattice dimension
  integer           :: Nx,Ny         !lattice grid dimensions
  integer           :: eqnloop       !Number of dmft loop for the equilibrium solution
  real(8)           :: beta0,xmu0,U0 !quench variables
  logical           :: iquench       !initial quench
  logical           :: update_wfftw  !update_wfftw=TT update WF using FFTw (iff Efield=0)
  logical           :: solve_wfftw   !solve_wfftw=TT solve Kadanof-Baym equations using FFTw (iff Efield=0)
  character(len=6)  :: method        !choose the perturbation theory method: IPT,SPT
  character(len=16) :: bath_type     !choose the shape of the BATH
  character(len=16) :: field_profile !choose the profile of the electric field
  real(8)           :: Vbath
  real(8)           :: Wbath         !Width of the BATH DOS
  real(8)           :: eps_error     !convergence error threshold
  integer           :: Nsuccess      !number of convergence success
  real(8)           :: weight        !mixing weight parameter
  real(8)           :: wmin,wmax     !min/max frequency
  real(8)           :: tmin,tmax     !min/max time
  logical           :: plot3D,fchi
  logical           :: solveEQ
  integer           :: size_cutoff
  logical           :: g0loc_guess    !use non-interacting local GF as guess.
  integer           :: iloop,nloop    !dmft loop variables
  real(8)           :: ts             !n.n./n.n.n. hopping amplitude
  real(8)           :: u              !local,non-local interaction 
  real(8)           :: dt,dtau        !time step
  real(8)           :: fmesh          !freq. step
  real(8)           :: beta           !inverse temperature
  real(8)           :: eps            !broadening
  integer           :: P,Q            !Uniform Power Mesh parameters
  !

  !FILES and VARS TO RESTART
  !=========================================================
  character(len=32)              :: irdFILE


  !FREQS & TIME ARRAYS:
  !=========================================================  
  real(8),dimension(:),allocatable :: wr,t,wm
  real(8),dimension(:),allocatable :: tau
  real(8),dimension(:),allocatable :: ftau


  !LATTICE (weight & dispersion) ARRAYS:
  !=========================================================  
  real(8),dimension(:),allocatable :: wt,epsik,sorted_epsik
  integer,allocatable,dimension(:) :: sorted_ik


  !ELECTRIC FIELD VARIABLES (& NML):
  !=========================================================  
  type(vect2D)    :: Ak,Ek         !Electric field vector potential and vector
  real(8)         :: Efield        !Electric field strength
  real(8)         :: Ex,Ey         !Electric field vectors as input
  real(8)         :: t0,t1         !turn on/off time, t0 also center of the pulse
  real(8)         :: w0,tau0       !parameters for pulsed light
  real(8)         :: omega0        !parameter for the Oscilatting field
  real(8)         :: E1            !Electric field strenght for the AC+DC case (tune to resonate)

  !EQUILIUBRIUM (and Wigner transformed) GREEN'S FUNCTION 
  !=========================================================
  !Frequency domain:
  type(keldysh_equilibrium_gf)        :: gf0
  type(keldysh_equilibrium_gf)        :: gf
  type(keldysh_equilibrium_gf)        :: sf
  real(8),dimension(:),allocatable    :: exa


  !MATSUBARA GREEN'S FUNCTION and n(k)
  !=========================================================
  complex(8),allocatable,dimension(:)  :: eq_G0iw
  real(8),allocatable,dimension(:)     :: eq_G0tau
  complex(8),allocatable,dimension(:)  :: eq_Siw
  real(8),allocatable,dimension(:)     :: eq_Stau
  complex(8),allocatable,dimension(:)  :: eq_Giw
  real(8),allocatable,dimension(:)     :: eq_Gtau
  real(8),allocatable,dimension(:)     :: eq_nk
  !


  !NON-EQUILIBRIUM FUNCTIONS:
  !=========================================================  
  !WEISS-FIELDS
  type(kbm_contour_gf) :: G0
  !SELF-ENERGY
  type(kbm_contour_gf) :: Sigma
  !LOCAL GF
  type(kbm_contour_gf) :: locG
  !Bath SELF-ENERGY
  type(kbm_contour_gf) :: S0
  !MOMENTUM-DISTRIBUTION
  real(8),allocatable,dimension(:,:)    :: nk



  !SUSCEPTIBILITY ARRAYS (in KADANOFF-BAYM)
  !=========================================================  
  real(8),allocatable,dimension(:,:,:,:) :: chi
  real(8),allocatable,dimension(:,:,:,:) :: chi_pm
  real(8),allocatable,dimension(:,:,:)   :: chi_dia


  !DATA DIRECTORY:
  !=========================================================
  character(len=32) :: data_dir,plot_dir



  !NAMELISTS:
  !=========================================================
  namelist/variables/dt,beta,U,Efield,Vbath,ts,nstep,nloop,eqnloop,& !global parameters
       eps_error,nsuccess,&                                          !convergence
       weight,&                                                      !mix
       Ex,Ey,t0,t1,tau0,w0,omega0,E1,field_profile,&                 !field
       Nx,Ny,&                                                       !k-points 
       method,update_wfftw,solve_wfftw,plot3D,g0loc_guess,solveEQ,&  !flags
       L,p,q,Lkreduced,Wbath,bath_type,eps,&                        !other parameters
       data_dir,plot_dir,fchi,&       
       irdFILE,&
       iquench,beta0,xmu0,U0


contains

  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine read_input_init(inputFILE)
    character(len=*)               :: inputFILE
    character(len=256),allocatable :: help_buffer(:)
    integer                        :: i
    logical                        :: control

    call version(revision)

    allocate(help_buffer(60))
    help_buffer=([&
         'NAME',&
         '  neqDMFT',&
         '',&
         'DESCRIPTION',&
         '  Run the non-equilibrium DMFT in presence of an external electric field E. ',&
         '  The field is controlled by few flags in the nml/cmd options. It can be ',&
         '  constant, pulse or switched off smoothly. Many more fields can be added by  ',&
         '  simply coding them in ELECTRIC_FIELD.f90. The executable read the file',&
         '  *inputFILE.ipt, if not found dump default values to a defualt file.',&
         '  ',&
         '  The output consist of very few data files that contain all the information,',&
         '  these are eventually read by a second program *get_data_neqDMFT to extract ',&
         '  all the relevant information.',&
         ' ',&
         '  In this version the impurity solver is: IPT',&
         ' ',&
         'OPTIONS (important)',&
         ' dt=[0.157080]            -- Time step for solution of KB equations',&
         ' beta=[100.0]             -- Inverse temperature ',&
         ' U=[6]                    -- Hubbard local interaction value',&
         ' Efield=[0]               -- Strenght of the electric field',&
         ' Vbath=[0]                -- Strenght of the coupling to bath (Lambda=Vbath^2/Wbath)',&
         ' Wbath=[10.0]             -- Bandwidth of the fermionic thermostat',&
         ' ts=[1]                   -- Hopping parameter',&
         ' nstep=[50]               -- Number of time steps: T_max = dt*nstep',&
         ' nloop=[30]               -- Maximum number of DMFT loops allowed (then exit)',&
         ' eps_error=[1.d-4]        -- Tolerance on convergence',&
         ' weight=[0.9]             -- Mixing parameter',&
         ' Nsuccess =[2]            -- Number of consecutive success for convergence to be true',&
         ' Ex=[1]                   -- X-component of the Electric field vector',&
         ' Ey=[1]                   -- Y-component of the Electric field vector',&
         ' t0=[0]                   -- Switching on time parameter for the Electric field',&
         ' t1=[10^6]                -- Switching off time parameter for the Electric field',&
         ' tau0=[1]                 -- Width of gaussian packect envelope for the impulsive Electric field',&
         ' w0=[20]                  -- Frequency of the of the impulsive Electric field',&
         ' field_profile=[dc]       -- Type of electric field profile (dc,ac,ac+dc,etc..)',&
         ' bath_type=[constant]     -- Fermionic thermostat type (constant,gaussian,bethe,etc..)',&         
         ' data_dir=[DATAneq]       -- Name of the directory containing data files',&
         ' plot_dir=[PLOT]          -- Name of the directory containing plot files',&
         ' fchi=[F]                 -- Flag for the calculation of the optical response',&
         ' L=[1024]                 -- A large number for whatever reason',&
         ' P=[5]                    -- Uniform Power mesh power-mesh parameter',&
         ' Q=[5]                    -- Uniform Power mesh uniform-mesh parameter',&
         ' eps=[0.05d0]             -- Broadening on the real-axis',&
         ' Nx=[50]                  -- Number of k-points along x-axis ',&
         ' Ny=[50]                  -- Number of k-points along y-axis ',&    
         '  '])
    call parse_cmd_help(help_buffer)


    dt=0.1d0
    beta= 100.d0
    U=4.d0
    Efield=0.d0
    Vbath=0.d0
    ts=1.d0
    nstep=100
    nloop=30
    eqnloop=30
    p=5
    q=5
    !CONVERGENCE:
    eps_error=1.d-4
    Nsuccess=2
    !MIX:
    weight=1.d0
    !FIELD:
    Ex=1.d0 
    Ey=0.d0
    E1=0.d0
    t0=0.d0 
    t1=1.d9
    tau0=0.d0
    w0=0.d0 
    field_profile='dc'
    !GRID k-POINTS:
    Nx=25 
    Ny=25
    !FLAGS:
    method='ipt'
    solveEQ=.false. 
    update_wfftw= .false.
    solve_wfftw= .false.
    plot3D= .false.
    fchi= .false.
    g0loc_guess= .false.
    !PARAMETERS:
    L=1024
    Lkreduced=200
    wbath=20.d0 
    bath_type='constant'
    eps=0.01d0
    !FILES&DIR:
    irdFILE='restartSigma'
    data_dir='DATAneq'
    plot_dir='PLOT'
    !QUENCH:
    iquench=.false.
    beta0=50.d0
    U0=6.d0
    xmu0=0.d0 

    inquire(file=adjustl(trim(inputFILE)),exist=control)
    if(control)then
       open(10,file=adjustl(trim(inputFILE)))
       read(10,nml=variables)
       close(10)
    else
       print*,"Can not find INPUT file"
       print*,"Dumping a default version in default."//trim(inputFILE)
       call dump_input_file("default.")
       call error("Can not find INPUT file, dumping a default version in default."//trim(inputFILE))
    endif

    ! !GLOBAL
    call parse_cmd_variable(dt ,"DT")
    call parse_cmd_variable(beta ,"BETA")
    call parse_cmd_variable(U ,"U")
    call parse_cmd_variable(Efield ,"EFIELD")
    call parse_cmd_variable(Vbath ,"VBATH")
    call parse_cmd_variable(ts ,"TS")
    call parse_cmd_variable(p ,"P")
    call parse_cmd_variable(q ,"Q")
    call parse_cmd_variable(nstep ,"NSTEP")
    call parse_cmd_variable(nloop ,"NLOOP")
    call parse_cmd_variable(eqnloop ,"EQNLOOP")
    !CONVERGENCE:
    call parse_cmd_variable(eps_error ,"EPS_ERROR")
    call parse_cmd_variable(Nsuccess ,"NSUCCESS")
    !MIX:
    call parse_cmd_variable(weight ,"WEIGHT")
    !FIELD:
    call parse_cmd_variable(Ex ,"EX")
    call parse_cmd_variable(Ey ,"EY")
    call parse_cmd_variable(t0 ,"T0")
    call parse_cmd_variable(t1 ,"T1")
    call parse_cmd_variable(tau0 ,"TAU0")
    call parse_cmd_variable(w0 ,"W0")
    call parse_cmd_variable(field_profile ,"FIELD_PROFILE")
    !GRID k-POINTS:
    call parse_cmd_variable(Nx ,"NX")
    call parse_cmd_variable(Ny ,"NY")
    !FLAGS:
    call parse_cmd_variable(method ,"METHOD")
    call parse_cmd_variable(solveEQ ,"SOLVEEQ")
    call parse_cmd_variable(update_wfftw ,"UPDATE_WFFTW")
    call parse_cmd_variable(solve_wfftw ,"SOLVE_WFFTW")
    call parse_cmd_variable(plot3D ,"PLOT3D")
    call parse_cmd_variable(fchi ,"FCHI")
    call parse_cmd_variable(g0loc_guess ,"G0LOC_GUESS")
    !PARAMETERS:
    call parse_cmd_variable(L ,"L")
    call parse_cmd_variable(Lkreduced ,"LKREDUCED")
    call parse_cmd_variable(wbath ,"WBATH")
    call parse_cmd_variable(bath_type ,"BATH_TYPE")
    call parse_cmd_variable(eps ,"EPS")
    !FILES&DIR:
    call parse_cmd_variable(data_dir,"DATA_DIR")
    call parse_cmd_variable(plot_dir,"PLOT_DIR")
    !QUENCH:
    call parse_cmd_variable(iquench ,"IQUENCH")
    call parse_cmd_variable(beta0 ,"BETA0")
    call parse_cmd_variable(U0 ,"U0")
    call parse_cmd_variable(xmu0 ,"XMU0") 

    Ltau=2*P*Q


    if(mpiID==0)then
       write(*,*)"CONTROL PARAMETERS"
       write(*,nml=variables)
       write(*,*)"--------------------------------------------"
       write(*,*)""
       call dump_input_file("used.")
    endif

    call create_data_dir(trim(data_dir))
    if(plot3D)call create_data_dir(trim(plot_dir))

  contains

    subroutine dump_input_file(prefix)
      character(len=*) :: prefix
      if(mpiID==0)then
         open(10,file=trim(adjustl(trim(prefix)))//adjustl(trim(inputFILE)))
         write(10,nml=variables)
         close(10)
      endif
    end subroutine dump_input_file
  end subroutine read_input_init
  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+----------------------------------------------------------------+
  !PURPOSE  : massive allocation of work array
  !+----------------------------------------------------------------+
  subroutine global_memory_allocation()
    integer          :: i
    real(8)          :: ex
    call msg("Allocating the memory")
    !Weiss-fields:
    call allocate_kbm_contour_gf(G0,Nstep,Ltau)
    !Interaction self-energies:
    call allocate_kbm_contour_gf(Sigma,Nstep,Ltau)
    !Local Green's functions:
    call allocate_kbm_contour_gf(locG,Nstep,Ltau)
    !Bath self-energies:
    call allocate_kbm_contour_gf(S0,Nstep,Ltau)
    ! allocate(S0gtr(-nstep:nstep),S0less(-nstep:nstep))
    ! allocate(S0gmix(0:Ltau,0:nstep),S0lmix(0:nstep,0:Ltau))
    !Momentum-distribution:
    allocate(nk(0:nstep,Lk),eq_nk(Lk))
    !Equilibrium/Wigner rotated Green's function
    call allocate_gf(gf0,nstep)
    call allocate_gf(gf,nstep)
    call allocate_gf(sf,nstep)
    !Matsubara Green's functions:
    allocate(eq_Stau(-Ltau:Ltau),eq_G0tau(-Ltau:Ltau),eq_Gtau(-Ltau:Ltau))
    allocate(eq_Siw(L),eq_G0iw(L),eq_Giw(L))
    !Susceptibility/Optical response
    if(fchi)allocate(chi(2,2,0:nstep,0:nstep))
    !Other:
    allocate(exa(-nstep:nstep))
    ex=-1.d0       
    do i=-nstep,nstep
       ex=-ex
       exa(i)=ex
    enddo
  end subroutine global_memory_allocation

  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine fftgf_iw2tau_upm(wm,gw,tm,gt,beta)
    integer                             :: i,j,N,L
    real(8),dimension(:)                :: wm
    complex(8),dimension(size(wm))      :: gw
    real(8),dimension(:)                :: tm
    real(8),dimension(size(tm))         :: gt
    real(8),dimension(2*size(wm))       :: tmpGw
    complex(8)                          :: tail,fg
    real(8)                             :: tau,beta,mues,At,foo
    !
    L=size(wm)
    N=size(tm)
    !
    mues =-real(gw(L),8)*wm(L)**2
    !
    tmpGw=(0.d0,0.d0)
    do i=1,L
       tail=-cmplx(mues,wm(i),8)/(mues**2+wm(i)**2)
       fg  = (0.d0,0.d0)
       if(i<=N)fg  = gw(i)-tail
       tmpGw(2*i)  = dimag(fg)
       tmpGw(2*i-1)= dreal(fg)
    enddo
    do i=1,N-1
       tau=tm(i)
       if(mues > 0.d0)then
          if((mues*beta) > 30.d0)then
             At = -exp(-mues*tau)
          else
             At = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
          endif
       else
          if((mues*beta) < -30.d0)then
             At = -exp(mues*(beta-tau))
          else
             At = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
          endif
       endif
       foo=0.d0
       do j=1,L
          foo=foo + sin(wm(j)*tau)*tmpGw(2*j) + cos(wm(j)*tau)*tmpGw(2*j-1)
       enddo
       gt(i) = foo*2.d0/beta + At
    enddo
    gt(N)=-(gt(1)+1.d0)
  end subroutine fftgf_iw2tau_upm

end module VARS_GLOBAL

