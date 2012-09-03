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
  real(8)           :: eps_error     !convergence error threshold
  integer           :: Nsuccess      !number of convergence success
  real(8)           :: weight        !mixing weight parameter
  real(8)           :: wmin,wmax     !min/max frequency
  logical           :: solve_eq      !Solve equilibrium Flag:
  logical           :: g0loc_guess   !use non-interacting local GF as guess.
  logical           :: plotVF,plot3D,fchi
  integer           :: size_cutoff


  !FILES TO RESTART
  !=========================================================
  character(len=32) :: irdG0file,irdNkfile,irdSlfile,irdSgfile


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
  real(8)         :: omega0        !parameter for the Oscilatting field

  !EQUILIUBRIUM/WIGNER TRANSFORMED GREEN'S FUNCTION 
  !=========================================================
  !Equilibrium initial conditions: Bath DOS, n(\e(k))
  real(8),allocatable,dimension(:)     :: irdNk,irdG0tau
  complex(8),allocatable,dimension(:)  :: irdG0w,irdG0iw

  !Frequency domain:
  type(keldysh_equilibrium_gf)        :: gf0
  type(keldysh_equilibrium_gf)        :: gf
  type(keldysh_equilibrium_gf)        :: sf


  !NON-EQUILIBRIUM GREEN'S FUNCTION: 4 = G^<,G^>
  !=========================================================  
  type keldysh_contour_gf
     complex(8),dimension(:,:),pointer  :: less,gtr
  end type keldysh_contour_gf

  interface assignment(=)
     module procedure keldysh_contour_gf_equality,keldysh_contour_gf_equality_
  end interface assignment(=)

  !MOMENTUM-DISTRIBUTION
  !=========================================================  
  real(8),allocatable,dimension(:,:)    :: nk

  !NON-INTERACTING
  !=========================================================  
  type(keldysh_contour_gf) :: G0

  !SELF-ENERGIES
  !=========================================================  
  !Sigma^V_k,mu(t,t`)
  type(keldysh_contour_gf) :: S0
  !Sigma^U(t,t`)
  type(keldysh_contour_gf) :: Sig

  !INTERACTING
  !=========================================================  
  type(keldysh_contour_gf) :: locG

  !IMPURITY
  !=========================================================  
  type(keldysh_contour_gf) :: impG


  !SUSCEPTIBILITY ARRAYS (in KADANOFF-BAYM)
  !=========================================================  
  real(8),allocatable,dimension(:,:,:,:) :: chi
  real(8),allocatable,dimension(:,:,:,:) :: chi_pm
  real(8),allocatable,dimension(:,:,:)   :: chi_dia



  !Other:
  real(8),dimension(:),allocatable    :: exa


  !NAMELISTS:
  !=========================================================
  namelist/variables/dt,beta,U,Efield,Vpd,ts,nstep,nloop,eps_error,nsuccess,weight,&
       Ex,Ey,t0,t1,tau0,w0,omega0,field_profile,Nx,Ny,&
       L,Ltau,Lmu,Lkreduced,Wbath,bath_type,eps,&
       method,irdeq,update_wfftw,solve_wfftw,plotVF,plot3D,fchi,equench,&
       iquench,beta0,xmu0,U0,irdG0file,irdnkfile,irdSlfile,irdSgfile,&
       solve_eq,g0loc_guess

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
         ' omega0=[pi]              -- Frequency of the of the Oscillating Electric field',&        
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
         ' irdG0file=[eqG0w.restart]-- ',&
         ' irdnkfile =[eqnk.restart]-- ',&
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
    call allocate_keldysh_contour_gf(G0,nstep)
    call allocate_keldysh_contour_gf(S0,nstep)
    call allocate_keldysh_contour_gf(Sig,nstep)
    call allocate_keldysh_contour_gf(locG,nstep)
    call allocate_keldysh_contour_gf(impG,nstep)
    allocate(nk(0:nstep,Lk),irdnk(Lk))

    call allocate_gf(gf0,nstep)
    call allocate_gf(gf,nstep)
    call allocate_gf(sf,nstep)

    if(fchi)allocate(chi(2,2,0:nstep,0:nstep))

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


  subroutine keldysh_contour_gf_equality(G1,G2)
    type(keldysh_contour_gf),intent(inout) :: G1
    type(keldysh_contour_gf),intent(in)    :: G2
    G1%less = G2%less
    G1%gtr = G2%gtr
  end subroutine keldysh_contour_gf_equality

  subroutine keldysh_contour_gf_equality_(G1,C)
    type(keldysh_contour_gf),intent(inout) :: G1
    complex(8),intent(in) :: C
    G1%less = C
    G1%gtr = C
  end subroutine keldysh_contour_gf_equality_



  !******************************************************************
  !******************************************************************
  !******************************************************************



  subroutine allocate_keldysh_contour_gf(G,N)
    type(keldysh_contour_gf) :: G
    integer                  :: i,j,N
    nullify(G%less,G%gtr)
    allocate(G%less(0:N,0:N),G%gtr(0:N,0:N))
    G%less=zero
    G%gtr =zero
  end subroutine allocate_keldysh_contour_gf



  !******************************************************************
  !******************************************************************
  !******************************************************************


  function build_keldysh_matrix_gf(G,N) result(matG)
    type(keldysh_contour_gf)              :: G
    complex(8),dimension(0:2*N+1,0:2*N+1) :: matG
    integer                               :: i,j,N
    forall(i=0:N,j=0:N)
       matG(i,j)         = step(t(i)-t(j))*G%gtr(i,j) + step(t(j)-t(i))*G%less(i,j)
       matG(i,N+1+j)     =-G%less(i,j)
       matG(N+1+i,j)     = G%gtr(i,j)
       matG(N+1+i,N+1+j) =-(step(t(i)-t(j))*G%less(i,j)+ step(t(j)-t(i))*G%gtr(i,j))
    end forall
  end function build_keldysh_matrix_gf


  !******************************************************************
  !******************************************************************
  !******************************************************************



  pure function fermi0(x,beta)
    real(8),intent(in) :: x,beta
    real(8) :: fermi0
    fermi0=fermi(x,beta)
  end function fermi0


end module VARS_GLOBAL

