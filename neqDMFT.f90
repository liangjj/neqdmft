!###################################################################
!PURPOSE  : Solve conduction band electrons driven 
! by electric field interacting with a resevoir of free 
! electrons at temperature T
!AUTHORS  : Adriano Amaricci && Cedric Weber
!###################################################################
program neqDMFT
  USE VARS_GLOBAL                 !global variables, calls to 3rd library 
  USE ELECTRIC_FIELD              !contains electric field init && routines
  USE BATH                        !contains bath inizialization
  USE EQUILIBRIUM                 !solves the equilibrium problem w/ IPT
  USE IPT_NEQ                     !performs the non-eq. IPT. Write Sigma
  USE UPDATE_WF                   !contains routines for WF update and printing.
  USE KADANOFBAYM                 !solves KB equations numerically to get k-sum
  implicit none

  logical :: converged

  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)

  !READ THE INPUT FILE (in vars_global):
  call read_input_init("inputFILE.in")

  !BUILD THE TIME,FREQUENCY GRIDS:
  include "grid_setup.f90"

  !BUILD THE 2D-SQUARE LATTICE STRUCTURE (in lib/square_lattice):
  Lk   = square_lattice_dimension(Nx,Ny)
  allocate(epsik(Lk),wt(Lk))
  wt   = square_lattice_structure(Lk,Nx,Ny)
  epsik= square_lattice_dispersion_array(Lk,ts)
  allocate(sorted_epsik(Lk),sorted_ik(Lk))
  sorted_epsik=epsik ; call sort_array(sorted_epsik,sorted_ik)
  if(mpiID==0)call get_free_dos(epsik,wt)

  !SET THE ELECTRIC FIELD (in electric_field):
  call set_efield_vector()

  !ALLOCATE FUNCTIONS IN THE MEMORY (in vars_global):
  call global_memory_allocation

  !BUILD THE  DISSIPATIVE BATH FUNCTIONS (in bath):
  call get_thermostat_bath()

  !SOLVE THE EQUILIBRIUM PROBLEM WITH IPT (in equilibrium):
  if(solveEQ)call solve_equilibrium_ipt()



  !START DMFT LOOP SEQUENCE:
  !==============================================================
  !initialize the run by guessing/reading the self-energy functions (in IPT_NEQ.f90):
  call neq_init_run

  iloop=0;converged=.false.
  do while(.not.converged);iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop",unit=6)
     !
     call neq_get_localgf        !-|(in kadanoff-baym)
     if(iloop==2)stop
     call neq_update_weiss_field !-|SELF-CONSISTENCY (in funx_neq)
     !
     call neq_solve_ipt          !-|IMPURITY SOLVER (in ipt_neq)
     !
     ! call print_observables      !(in funx_neq)
     converged = convergence_check()
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     !
     call end_loop()
  enddo
  call msg("BRAVO")
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)

end PROGRAM neqDMFT
