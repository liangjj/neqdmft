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
  implicit none

  !Version revision
  include "revision.inc"

  !Gloabl  variables
  !=========================================================
  integer,protected :: Lmu           !# of bath energies
  integer,protected :: Ltau          !Imaginary time slices
  integer,protected :: L             !a big number
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
  real(8)           :: wmin,wmax     !
  character(len=32) :: irdG0file,irdNkfile


  !FREQS & TIME ARRAYS:
  !=========================================================  
  real(8),dimension(:),allocatable :: wr,t,wm,tau


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


  !EQUILIUBRIUM/WIGNER TRANSFORMED GREEN'S FUNCTION 
  !=========================================================  
  !Equilibrium initial conditions: Bath DOS, n(\e(k))
  real(8),allocatable,dimension(:)     :: irdNk
  complex(8),allocatable,dimension(:)  :: irdG0w

  !Frequency domain:
  type(keldysh_equilibrium_gf)        :: gf0,gf,sf

  real(8),dimension(:),allocatable    :: exa


  !NON-EQUILIBRIUM GREEN'S FUNCTION: 4 = G^<,G^>
  !=========================================================  
  !NON-INTERACTING
  complex(8),allocatable,dimension(:,:) :: G0gtr
  complex(8),allocatable,dimension(:,:) :: G0less
  real(8),allocatable,dimension(:,:)    :: nk

  !SELF-ENERGIES
  !=========================================================  
  !Sigma^V_k,mu(t,t`)
  complex(8),allocatable,dimension(:)   :: S0gtr
  complex(8),allocatable,dimension(:)   :: S0less
  !Sigma^U(t,t`)
  complex(8),allocatable,dimension(:,:) :: Sgtr
  complex(8),allocatable,dimension(:,:) :: Sless

  !INTERACTING
  !=========================================================  
  complex(8),allocatable,dimension(:,:) :: locGgtr
  complex(8),allocatable,dimension(:,:) :: locGless

  !IMPURITY
  !=========================================================  
  complex(8),allocatable,dimension(:,:) :: impGgtr
  complex(8),allocatable,dimension(:,:) :: impGless


  !SUSCEPTIBILITY ARRAYS (in KADANOFF-BAYM)
  !=========================================================  
  real(8),allocatable,dimension(:,:,:,:) :: chi
  real(8),allocatable,dimension(:,:,:,:) :: chi_pm
  real(8),allocatable,dimension(:,:,:)   :: chi_dia


  !FLAGS
  !=========================================================  
  logical              :: plotVF,plot3D,fchi
  integer              :: size_cutoff


  !NAMELISTS:
  !=========================================================
  namelist/variables/dt,beta,U,Efield,Vpd,ts,nstep,nloop,eps_error,nsuccess,&
       Ex,Ey,t0,t1,tau0,w0,field_profile,Nx,Ny,&
       L,Ltau,Lmu,Lkreduced,Wbath,bath_type,eps,irdG0file,irdnkfile,omp_num_threads,&
       method,irdeq,update_wfftw,solve_wfftw,plotVF,plot3D,fchi,equench,&
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
    logical          :: control

    call version(revision)

    allocate(help_buffer(59))
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
         ' dt=[0.157080]            -- ',&
         ' beta=[100.0]             -- ',&
         ' U=[6]                    -- ',&
         ' Efield=[0]               -- ',&
         ' Vpd=[0]                  -- ',&
         ' ts=[1]                   -- ',&
         ' nstep=[50]               -- ',&
         ' nloop=[30]               -- ',&
         ' eps_error=[1.d-4]        -- ',&
         ' Nsuccess =[2]            -- ',&
         ' Ex=[1]                   -- ',&
         ' Ey=[1]                   -- ',&
         ' t0=[0]                   -- ',&
         ' t1=[dt*(nstep+10)]       -- ',&
         ' tau0=[1]                 -- ',&
         ' w0=[20]                  -- ',&
         ' field_profile=[constant] -- ',&
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
         ' irdG0file=[eqG0w.restart]-- ',&
         ' irdnkfile =[eqnk.restart]-- ',&
         ' omp_num_threads=[1]      -- ',&
         ' Nx=[50]      -- ',&
         ' Ny=[50]      -- ',&    
         ' iquench=[F]  -- ',&
         ' beta0=[100]  -- ',&
         ' U0=[6]       -- ',&
         ' xmu0=[0]     -- ',& 
         '  '])
    call parse_cmd_help(help_buffer)

    include "nml_default_values.f90"
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
    if(present(printf).AND.printf.eq..true.)call dump_input_file("used.")

    !SET OMP THREADS NUMBER
    call omp_set_num_threads(omp_num_threads)

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
    allocate(G0gtr(0:nstep,0:nstep),G0less(0:nstep,0:nstep))
    allocate(S0gtr(-nstep:nstep),S0less(-nstep:nstep))
    allocate(Sgtr(0:nstep,0:nstep),Sless(0:nstep,0:nstep))
    allocate(locGless(0:nstep,0:nstep),locGgtr(0:nstep,0:nstep))
    allocate(impGless(0:nstep,0:nstep),impGgtr(0:nstep,0:nstep))
    allocate(nk(0:nstep,Lk))

    call allocate_gf(gf0,nstep)
    call allocate_gf(gf,nstep)
    call allocate_gf(sf,nstep)

    if(fchi)then
       allocate(chi(2,2,0:nstep,0:nstep))
    endif

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




  function fermi0(x,beta)
    real(8),intent(in) :: x,beta
    real(8) :: fermi0
    fermi0=fermi(x,beta)
  end function fermi0


end module VARS_GLOBAL

