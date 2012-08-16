!#####################################################################
!     Program  : VARS_GLOBAL
!     TYPE     : Module
!     PURPOSE  : Defines the global variables used thru all the code
!     AUTHORS  : Adriano Amaricci
!#####################################################################
MODULE VARS_GLOBAL
  !SciFor library
  USE COMMON_VARS
  USE GREENFUNX
  USE CHRONOBAR
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
  integer,protected :: Lmu           !# of bath energies
  integer,protected :: Ltau          !Imaginary time slices
  integer           :: L             !a big number
  integer           :: Lk            !total lattice  dimension
  integer           :: Lkreduced     !reduced lattice dimension
  integer           :: Nx,Ny         !lattice grid dimensions
  integer           :: nstep         !Number of Time steps
  real(8)           :: beta0,xmu0,U0 !quench variables
  logical           :: iquench       !quench flag
  logical           :: Equench       !initial condition with (T) or without (F) electric field
  logical           :: irdeq         !irdeq=TT read inputs from equilbrium solution
  logical           :: update_wfftw  !update_wfftw=TT update WF using FFTw (iff Efield=0)
  logical           :: solve_wfftw   !solve_wfftw=TT solve Kadanof-Baym equations using FFTw (iff Efield=0)
  character(len=6)  :: method        !choose the perturbation theory method: IPT,SPT
  character(len=16) :: bath_type     !choose the shape of the BATH
  character(len=16) :: field_profile !choose the profile of the electric field
  real(8)           :: Wbath         !Width of the BATH DOS
  real(8)           :: eps_error
  integer           :: Nsuccess
  real(8)           :: weight        !mixing weight parameter
  real(8)           :: wmin,wmax     !
  character(len=32) :: irdG0wfile,irdnkfile !the bath GF file
  logical           :: plotVF,plot3D,fchi
  integer           :: size_cutoff
  logical           :: solve_eq
  logical           :: g0loc_guess





  !FREQS & TIME ARRAYS:
  !=========================================================  
  real(8),dimension(:),allocatable :: wr,t,wm,tau,taureal
  real(8)                          :: dtaureal

  !LATTICE (weight & dispersion) ARRAYS:
  !=========================================================  
  real(8),dimension(:),allocatable :: wt,epsik


  !ELECTRIC FIELD VARIABLES (& NML):
  !=========================================================  
  type(vect2D)    :: Ek            !Electric field vector
  real(8)         :: Efield        !Electric field strength
  real(8)         :: Ex,Ey         !Electric field vectors as input
  real(8)         :: t0,t1         !turn on/off time, t0 also center of the pulse
  real(8)         :: w0,tau0       !parameters for pulsed light


  !EQUILIUBRIUM (and Wigner transformed) GREEN'S FUNCTION 
  !=========================================================
  type(keldysh_equilibrium_gf)        :: gf0,gf,sf
  real(8),dimension(:),allocatable    :: exa


  !INITIAL CONDITIONS: BATH DOS, N(\e(k)), Matsubara Self-energy
  !=========================================================
  real(8),allocatable,dimension(:)     :: eq_nk
  real(8),allocatable,dimension(:)     :: eq_Stau
  complex(8),allocatable,dimension(:)  :: eq_Siw
  complex(8),allocatable,dimension(:)  :: eq_G0w


  !NON-EQUILIBRIUM FUNCTIONS:
  !=========================================================  
  !WEISS-FIELDS
  complex(8),allocatable,dimension(:,:) :: G0gtr,G0gmix
  complex(8),allocatable,dimension(:,:) :: G0less,G0lmix
  real(8),allocatable,dimension(:,:)    :: G0mat

  !MOMENTUM-DISTRIBUTION
  real(8),allocatable,dimension(:,:)    :: nk

  !SELF-ENERGIES
  !Bath:
  !Sigma_bath^X(t,t`) X=>,<,\lmix,\gmix,Matsubara
  complex(8),allocatable,dimension(:)   :: S0gtr,S0less
  complex(8),allocatable,dimension(:,:) :: S0gmix,S0lmix
  !Interaction:
  !Sigma_U^X(t,t`) X=>,<,\lmix,\gmix
  complex(8),allocatable,dimension(:,:) :: Sgtr,Sgmix
  complex(8),allocatable,dimension(:,:) :: Sless,Slmix
  real(8),allocatable,dimension(:,:)    :: Smat

  !LOCAL GF
  complex(8),allocatable,dimension(:,:) :: locGgtr,locGgmix
  complex(8),allocatable,dimension(:,:) :: locGless,locGlmix
  real(8),allocatable,dimension(:,:)    :: locGmat

  ! !IMPURITY
  ! complex(8),allocatable,dimension(:,:) :: impGgtr,impGless
  ! complex(8),allocatable,dimension(:,:) :: impGgmix,impGlmix



  !SUSCEPTIBILITY ARRAYS (in KADANOFF-BAYM)
  !=========================================================  
  real(8),allocatable,dimension(:,:,:,:) :: chi
  real(8),allocatable,dimension(:,:,:,:) :: chi_pm
  real(8),allocatable,dimension(:,:,:)   :: chi_dia


  !FLAGS
  !=========================================================  



  !NAMELISTS:
  !=========================================================
  namelist/variables/dt,beta,U,Efield,Vpd,ts,nstep,nloop,eps_error,nsuccess,weight,&
       Ex,Ey,t0,t1,tau0,w0,field_profile,Nx,Ny,&
       L,Ltau,Lmu,Lkreduced,Wbath,bath_type,eps,&
       method,irdeq,update_wfftw,solve_wfftw,plotVF,plot3D,fchi,equench,&
       solve_eq,g0loc_guess,&
       irdNkfile,irdG0wfile,&
       iquench,beta0,xmu0,U0


contains

  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine read_input_init(inputFILE,printf)
    character(len=*) :: inputFILE
    integer          :: i
    logical,optional :: printf
    logical          :: lprint
    logical          :: control

    lprint=.false. ; if(present(printf))lprint=printf
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
         'OPTIONS',&
         ' dt=[0.157080]            -- Time step for solution of KB equations',&
         ' beta=[100.0]             -- Inverse temperature ',&
         ' U=[6]                    -- Hubbard local interaction value',&
         ' Efield=[0]               -- Strenght of the electric field',&
         ' Vpd=[0]                  -- Strenght of the coupling to bath (Lambda=Vpd^2/Wbath)',&
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
         ' field_profile=[constant] -- Type of electric field profile (constant,gaussian,ramp)',&
         ' irdeq=[F]        -- ',&
         ' method=[ipt]     -- ',&
         ' update_wfftw=[F] -- ',&
         ' solve_wfftw =[F] -- ',&
         ' plotVF=[F]       -- ',&
         ' plot3D=[F]       -- ',&
         ' fchi=[F]         -- ',&
         ' equench=[F]      -- ',&
         ' L=[1024]         -- ',&
         ' Ltau=[32]        -- ',&
         ' Lmu=[2048]       -- ',&
         ' Lkreduced=[200]  -- ',&
         ' wbath=[10.0]     -- ',&
         ' bath_type=[constant] -- ',&
         ' eps=[0.05d0]         -- ',&
         ' irdnkfile =[eqnk.restart]-- ',&
         ' irdG0wfile=[eqG0w.restart]-- ',&
         ' Nx=[50]      -- ',&
         ' Ny=[50]      -- ',&    
         ' iquench=[F]  -- ',&
         ' beta0=[100]  -- ',&
         ' U0=[6]       -- ',&
         ' xmu0=[0]     -- ',& 
         '  '])
    call parse_cmd_help(help_buffer)

    !DEFAULT NML VARIABLES VALUES:
    dt      = 0.157080
    beta    = 100.0
    U       = 6.0
    Efield  = 0.0
    Vpd     = 0.0
    ts      = 1.0
    nstep   = 50
    nloop   = 30
    eps_error= 1.d-4
    Nsuccess = 2
    weight  = 0.9d0

    Ex  = 1.d0
    Ey  = 1.d0
    t0  = 0.d0
    t1  = 1000000.d0              !infinite time SHIT!!
    tau0= 1.d0
    w0  = 20.d0
    field_profile='constant'

    method        = 'ipt'
    irdeq         = .false.
    update_wfftw  = .false.
    solve_wfftw   = .false.
    plotVF        = .false.
    plot3D        = .false.
    fchi          = .false.
    equench       = .false.
    solve_eq      = .true.
    g0loc_guess   = .false.
    iquench       = .false.

    L          = 1024
    Ltau       = 32
    Lmu        = 2048
    Lkreduced  = 200
    wbath      = 20.0
    bath_type  = "constant"
    eps        = 0.04d0
    irdnkfile  = "eqnk.restart"
    irdG0wfile = "eqG0w.restart"

    Nx = 25
    Ny = 25    

    beta0   = 100.0
    U0      = 6.0
    xmu0    = 0.0    

    omp_num_threads =1

    inquire(file=adjustl(trim(inputFILE)),exist=control)
    if(control)then
       open(10,file=adjustl(trim(inputFILE)))
       read(10,nml=variables)
       close(10)
    else
       print*,"Can not find INPUT file"
       print*,"Dumping a default version in default."//trim(inputFILE)
       call dump_input_file("default.")
       call abort("Can not find INPUT file, dumping a default version in default."//trim(inputFILE))
    endif
    include "nml_read_cml.f90"

    write(*,*)"CONTROL PARAMETERS"
    write(*,nml=variables)
    write(*,*)"--------------------------------------------"
    write(*,*)""
    if(lprint)call dump_input_file("used.")

    return
  contains
    subroutine dump_input_file(prefix)
      character(len=*) :: prefix
      open(10,file=trim(adjustl(trim(prefix)))//adjustl(trim(inputFILE)))
      write(10,nml=variables)
      close(10)
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
    allocate(G0gtr(0:nstep,0:nstep),G0less(0:nstep,0:nstep))
    allocate(G0gmix(0:Ltau,0:nstep),G0lmix(0:nstep,0:Ltau))
    allocate(G0mat(0:Ltau,0:Ltau))

    !Bath self-energies:
    allocate(S0gtr(-nstep:nstep),S0less(-nstep:nstep))
    allocate(S0gmix(0:Ltau,0:nstep),S0lmix(0:nstep,0:Ltau))

    !Interaction self-energies:
    allocate(Sgtr(0:nstep,0:nstep),Sless(0:nstep,0:nstep))
    allocate(Sgmix(0:Ltau,0:nstep),Slmix(0:nstep,0:Ltau))
    allocate(Smat(0:Ltau,0:Ltau))

    !Local Green's functions:
    allocate(locGgtr(0:nstep,0:nstep),locGless(0:nstep,0:nstep))
    allocate(locGgmix(0:Ltau,0:nstep),locGlmix(0:nstep,0:Ltau))
    allocate(locGmat(0:Ltau,0:Ltau))

    !Momentum-distribution:
    allocate(nk(0:nstep,Lk))

    !Equilibrium/Wigner rotated Green's function
    call allocate_gf(gf0,nstep)
    call allocate_gf(gf,nstep)
    call allocate_gf(sf,nstep)

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



  pure function fermi0(x,beta)
    real(8),intent(in) :: x,beta
    real(8) :: fermi0
    fermi0=fermi(x,beta)
  end function fermi0


end module VARS_GLOBAL

