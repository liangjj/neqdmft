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
  USE FUNX_NEQ                    !contains routines for WF update and printing.
  USE KADANOFBAYM                 !solves KB equations numerically to get k-sum
  implicit none

  logical :: converged


  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)

  call read_input_init("inputFILE.in",printf=.true.)
  include "grid_setup.f90"
  include "build_square_lattice.f90"

  !SET THE ELECTRIC FIELD:use constant field by default
  Ek = set_efield_vector(Ex,Ey)
  if(mpiID==0)call print_Afield_form(t(0:nstep))

  !STARTS THE REAL WORK:
  call global_memory_allocation !allocate functions in the memory

  if(solve_eq)call solve_equilibrium_ipt()

  call get_thermostat_bath()    !get the dissipative bath functions
  !call neq_guess_weiss_field   !guess/read the first Weiss field:

  call neq_init_sigma           !initialize (read) the first Sigma. Guess.

  iloop=0;converged=.false.
  do while (.not.converged);iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop",unit=6)
     !
     call neq_get_localgf !-|
     call neq_update_weiss_field  !-|SELF-CONSISTENCY
     !
     call neq_solve_ipt           !-|IMPURITY SOLVER
     !
     call print_observables
     converged = convergence_check()
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     !
     call end_loop()
  enddo
  call msg("BRAVO")
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)

end PROGRAM neqDMFT
