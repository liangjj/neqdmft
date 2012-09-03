MODULE CONTOUR_GF
  USE COMMON_VARS
  USE TOOLS
  USE IOTOOLS
  USE MPI
  implicit none
  private

  type :: keldysh_contour_gf
     complex(8),dimension(:,:),pointer  :: less,gtr
     logical                            :: status
     integer                            :: N
  end type keldysh_contour_gf

  interface keldysh_contour_gf_sum
     module procedure keldysh_contour_gf_sum_d,keldysh_contour_gf_sum_z
  end interface keldysh_contour_gf_sum

  public :: keldysh_contour_gf
  public :: allocate_keldysh_contour_gf
  public :: deallocate_keldysh_contour_gf
  public :: write_keldysh_contour_gf,read_keldysh_contour_gf,plot_keldysh_contour_gf
  public :: mpi_reduce_keldysh_contour_gf
  public :: mpi_bcast_keldysh_contour_gf
  public :: keldysh_contour_gf_sum



  type :: kbm_contour_gf
     complex(8),dimension(:,:),pointer  :: less,gtr
     complex(8),dimension(:,:),pointer  :: lmix,gmix
     real(8),dimension(:,:),pointer     :: mats
     logical                            :: status=.false.
     integer                            :: N=0,L=0
  end type kbm_contour_gf

  interface operator(*)
     module procedure &
          keldysh_contour_gf_scalarL_d,keldysh_contour_gf_scalarL_c,&
          keldysh_contour_gf_scalarR_d,keldysh_contour_gf_scalarR_c,&
          kbm_contour_gf_scalarL_d,kbm_contour_gf_scalarL_c,&
          kbm_contour_gf_scalarR_d,kbm_contour_gf_scalarR_c
  end interface operator(*)

  interface assignment(=)
     module procedure &
          keldysh_contour_gf_equality,&
          keldysh_contour_gf_equality_,&
          kbm_contour_gf_equality_
  end interface assignment(=)

  interface kbm_contour_gf_sum
     module procedure kbm_contour_gf_sum_d,kbm_contour_gf_sum_z
  end interface kbm_contour_gf_sum

  public :: kbm_contour_gf
  public :: allocate_kbm_contour_gf
  public :: deallocate_kbm_contour_gf
  public :: write_kbm_contour_gf,read_kbm_contour_gf,plot_kbm_contour_gf
  public :: mpi_reduce_kbm_contour_gf
  public :: mpi_bcast_kbm_contour_gf
  public :: assignment(=)
  public :: operator(*)
  public :: kbm_contour_gf_sum

contains


  !################################################################################
  !########### KELDYSH CONTOUR GREEN'S FUNCTION (REAL-TIME ONLY) ##################
  !################################################################################

  subroutine keldysh_contour_gf_equality(G1,G2)
    type(keldysh_contour_gf),intent(inout) :: G1
    type(keldysh_contour_gf),intent(in)    :: G2
    G1%less = G2%less
    G1%gtr = G2%gtr
  end subroutine keldysh_contour_gf_equality
  !--------------------------------------------

  !--------------------------------------------
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
    G%N=N
    allocate(G%less(0:N,0:N),G%gtr(0:N,0:N))
    G%less=zero
    G%gtr =zero
    G%status=.true.
  end subroutine allocate_keldysh_contour_gf


  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine deallocate_keldysh_contour_gf(G)
    type(keldysh_contour_gf) :: G
    deallocate(G%less,G%gtr)
    G%N=0
    G%status=.false.
  end subroutine deallocate_keldysh_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine write_keldysh_contour_gf(G,file)
    type(keldysh_contour_gf) :: G
    character(len=*)         :: file
    integer                  :: N
    N=G%N+1
    if( (size(G%less)/=N**2) .OR. (size(G%gtr)/=N**2) )&
         call error("ERROR contour_gf/write_keldysh_contour_gf: wrong dimensions")
    call splot(trim(file)//"_less.data",G%less(0:,0:))
    call splot(trim(file)//"_gtr.data",G%gtr(0:,0:))
  end subroutine write_keldysh_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine read_keldysh_contour_gf(G,file)
    type(keldysh_contour_gf) :: G
    character(len=*)     :: file
    integer              :: N,L
    N=G%N+1
    if( (size(G%less)/=N**2) .OR. (size(G%gtr)/=N**2) )&
         call error("ERROR contour_gf/write_keldysh_contour_gf: wrong dimensions")
    call sread(trim(file)//"_less.data",G%less(0:,0:))
    call sread(trim(file)//"_gtr.data",G%gtr(0:,0:))
  end subroutine read_keldysh_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine plot_keldysh_contour_gf(G,t,file)
    type(keldysh_contour_gf)  :: G
    character(len=*)      :: file
    real(8),dimension(0:) :: t
    integer               :: N
    N=G%N+1
    if( (size(G%less)/=N**2) .OR. (size(G%gtr)/=N**2) )&
         call error("ERROR contour_gf/plot_keldysh_contour_gf: wrong dimensions")
    call splot(reg_filename(file)//"_less_t_t",t(0:),t(0:),G%less(0:,0:))
    call splot(reg_filename(file)//"_gtr_t_t",t(0:),t(0:),G%gtr(0:,0:))
  end subroutine plot_keldysh_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine keldysh_contour_gf_sum_d(Ga,Gb,C)
    type(keldysh_contour_gf)  :: Ga
    type(keldysh_contour_gf)  :: Gb
    real(8)               :: C
    Ga%less(0:,0:) = Ga%less(0:,0:) + Gb%less(0:,0:)*C
    Ga%gtr(0:,0:)  = Ga%gtr(0:,0:)  + Gb%gtr(0:,0:)*C
  end subroutine keldysh_contour_gf_sum_d

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine keldysh_contour_gf_sum_z(Ga,Gb,C)
    type(keldysh_contour_gf)  :: Ga
    type(keldysh_contour_gf)  :: Gb
    complex(8)            ::C
    Ga%less(0:,0:) = Ga%less(0:,0:) + Gb%less(0:,0:)*C
    Ga%gtr(0:,0:)  = Ga%gtr(0:,0:)  + Gb%gtr(0:,0:)*C
  end subroutine keldysh_contour_gf_sum_z

  !******************************************************************
  !******************************************************************
  !******************************************************************

  function keldysh_contour_gf_scalarL_d(C,G) result(F)
    real(8),intent(in) :: C
    type(keldysh_contour_gf),intent(in) :: G
    type(keldysh_contour_gf) :: F
    F%less(0:,0:)= C*G%less(0:,0:)
    F%gtr(0:,0:) = C*G%gtr(0:,0:)
  end function keldysh_contour_gf_scalarL_d

  !******************************************************************
  !******************************************************************
  !******************************************************************

  function keldysh_contour_gf_scalarL_c(C,G) result(F)
    complex(8),intent(in) :: C
    type(keldysh_contour_gf),intent(in) :: G
    type(keldysh_contour_gf) :: F
    F%less(0:,0:)=C*G%less(0:,0:)
    F%gtr(0:,0:)=C*G%gtr(0:,0:)
  end function keldysh_contour_gf_scalarL_c

  !******************************************************************
  !******************************************************************
  !******************************************************************

  function keldysh_contour_gf_scalarR_d(G,C) result(F)
    real(8),intent(in) :: C
    type(keldysh_contour_gf),intent(in) :: G
    type(keldysh_contour_gf) :: F
    F%less(0:,0:)=G%less(0:,0:)*C
    F%gtr(0:,0:)=G%gtr(0:,0:)*C
  end function keldysh_contour_gf_scalarR_d

  !******************************************************************
  !******************************************************************
  !******************************************************************

  function keldysh_contour_gf_scalarR_c(G,C) result(F)
    complex(8),intent(in) :: C
    type(keldysh_contour_gf),intent(in) :: G
    type(keldysh_contour_gf) :: F
    F%less(0:,0:)=G%less(0:,0:)*C
    F%gtr(0:,0:)=G%gtr(0:,0:)*C
  end function keldysh_contour_gf_scalarR_c

  !******************************************************************
  !******************************************************************
  !******************************************************************

  !This a rough implementation: just a shortcut:
  subroutine mpi_reduce_keldysh_contour_gf(tmpG,G)
    type(keldysh_contour_gf) :: tmpG
    type(keldysh_contour_gf) :: G
    if( .not.g%status )call error("ERROR contour_gf/mpi_reduce_keldysh_contour_gf: object function not allocated.")
    if( (G%N /= tmpG%N) )call error("ERROR contour_gf/mpi_reduce_keldysh_contour_gf: wrong dimensions.")
    call MPI_REDUCE(tmpG%less(0:,0:),G%less(0:,0:),size(tmpG%less(0:,0:)),&
         MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_REDUCE(tmpG%gtr(0:,0:),G%gtr(0:,0:),size(tmpG%gtr(0:,0:)),&
         MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
  end subroutine mpi_reduce_keldysh_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine mpi_bcast_keldysh_contour_gf(G)
    type(keldysh_contour_gf),intent(inout) :: G
    if( .not.g%status )call error("ERROR contour_gf/mpi_bcast_keldysh_contour_gf: object function not allocated.")    
    call MPI_BCAST(G%less(0:,0:),size(G%less(0:,0:)),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(G%gtr(0:,0:),size(G%gtr(0:,0:)),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
  end subroutine mpi_bcast_keldysh_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************





  !################################################################################
  !############ KADANOFF-BAYM-MATSUBARA CONTOUR GREEN'S FUNCTION ##################
  !################################################################################
  subroutine allocate_kbm_contour_gf(G,N,L)
    type(kbm_contour_gf)  :: G
    integer                 :: i,j,N,L
    nullify(G%less,G%gtr,G%lmix,G%gmix,G%mats)
    G%N=N
    G%L=L
    allocate(G%less(0:N,0:N)) ; G%less=zero
    allocate(G%gtr(0:N,0:N))  ; G%gtr=zero
    allocate(G%lmix(0:N,0:L)) ; G%lmix=zero
    allocate(G%gmix(0:L,0:N)) ; G%gmix=zero
    allocate(G%mats(0:L,0:L)) ; G%mats=zero
    G%status=.true.
  end subroutine allocate_kbm_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine deallocate_kbm_contour_gf(G)
    type(kbm_contour_gf) :: G
    deallocate(G%less,G%gtr,G%lmix,G%gmix,G%mats)
    G%N=0
    G%L=0
    G%status=.false.
  end subroutine deallocate_kbm_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine write_kbm_contour_gf(G,file)
    type(kbm_contour_gf) :: G
    character(len=*)     :: file
    integer              :: N,L
    N=G%N+1 ; L=G%L+1
    if( (size(G%less)/=N**2) .OR. (size(G%gtr)/=N**2) )&
         call error("ERROR contour_gf/write_kbm_contour_gf: wrong dimensions 1")
    if( (size(G%lmix)/=N*L)  .OR. (size(G%gmix)/=N*L) )&
         call error("ERROR contour_gf/write_kbm_contour_gf: wrong dimensions 2")
    if( size(G%mats)/=L*L)&
         call error("ERROR contour_gf/write_kbm_contour_gf: wrong dimensions 3")
    call splot(trim(file)//"_less.data",G%less(0:,0:))
    call splot(trim(file)//"_gtr.data",G%gtr(0:,0:))
    call splot(trim(file)//"_lmix.data",G%lmix(0:,0:))
    call splot(trim(file)//"_gmix.data",G%gmix(0:,0:))
    call splot(trim(file)//"_mats.data",G%mats(0:,0:))
  end subroutine write_kbm_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine read_kbm_contour_gf(G,file)
    type(kbm_contour_gf) :: G
    character(len=*)     :: file
    integer              :: N,L
    N=G%N+1 ; L=G%L+1
    if( (size(G%less)/=N**2) .OR. (size(G%gtr)/=N**2) )&
         call error("ERROR contour_gf/write_kbm_contour_gf: wrong dimensions 1")
    if( (size(G%lmix)/=N*L)  .OR. (size(G%gmix)/=N*L) )&
         call error("ERROR contour_gf/write_kbm_contour_gf: wrong dimensions 2")
    if( size(G%mats)/=L*L)&
         call error("ERROR contour_gf/write_kbm_contour_gf: wrong dimensions 3")
    call sread(trim(file)//"_less.data",G%less(0:,0:))
    call sread(trim(file)//"_gtr.data",G%gtr(0:,0:))
    call sread(trim(file)//"_lmix.data",G%lmix(0:,0:))
    call sread(trim(file)//"_gmix.data",G%gmix(0:,0:))
    call sread(trim(file)//"_mats.data",G%mats(0:,0:))
  end subroutine read_kbm_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine plot_kbm_contour_gf(G,t,tau,file)
    type(kbm_contour_gf)  :: G
    character(len=*)      :: file
    real(8),dimension(0:) :: t,tau
    integer               :: N,L
    N=G%N+1 ; L=G%L+1
    if( (size(G%less)/=N**2) .OR. (size(G%gtr)/=N**2) )&
         call error("ERROR contour_gf/plot_kbm_contour_gf: wrong dimensions 1")
    if( (size(G%lmix)/=N*L)  .OR. (size(G%gmix)/=N*L) )&
         call error("ERROR contour_gf/plot_kbm_contour_gf: wrong dimensions 2")
    if( size(G%mats)/=L*L)&
         call error("ERROR contour_gf/plot_kbm_contour_gf: wrong dimensions 3")
    call splot(reg_filename(file)//"_less_t_t",t(0:),t(0:),G%less(0:,0:))
    call splot(reg_filename(file)//"_gtr_t_t",t(0:),t(0:),G%gtr(0:,0:))
    call splot(reg_filename(file)//"_lmix_t_tau",t(0:),tau(0:),G%lmix(0:,0:))
    call splot(reg_filename(file)//"_gmix_tau_t",tau(0:),t(0:),G%gmix(0:,0:))
    call splot(reg_filename(file)//"_mats_tau_tau",tau(0:),tau(0:),G%mats(0:,0:))
  end subroutine plot_kbm_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine kbm_contour_gf_equality_(G1,C)
    type(kbm_contour_gf),intent(inout) :: G1
    complex(8),intent(in) :: C
    G1%less(0:,0:) = C
    G1%gtr(0:,0:) = C
    G1%lmix(0:,0:) = C
    G1%gmix(0:,0:) = C
    G1%mats(0:,0:) = C
  end subroutine kbm_contour_gf_equality_

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine kbm_contour_gf_sum_d(Ga,Gb,C)
    type(kbm_contour_gf)  :: Ga
    type(kbm_contour_gf)  :: Gb
    real(8)               :: C
    Ga%less(0:,0:) = Ga%less(0:,0:) + Gb%less(0:,0:)*C
    Ga%gtr(0:,0:)  = Ga%gtr(0:,0:)  + Gb%gtr(0:,0:)*C
    Ga%lmix(0:,0:) = Ga%lmix(0:,0:) + Gb%lmix(0:,0:)*C
    Ga%gmix(0:,0:) = Ga%gmix(0:,0:) + Gb%gmix(0:,0:)*C
    Ga%mats(0:,0:) = Ga%mats(0:,0:) + Gb%mats(0:,0:)*C
  end subroutine kbm_contour_gf_sum_d

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine kbm_contour_gf_sum_z(Ga,Gb,C)
    type(kbm_contour_gf)  :: Ga
    type(kbm_contour_gf)  :: Gb
    complex(8)            ::C
    Ga%less(0:,0:) = Ga%less(0:,0:) + Gb%less(0:,0:)*C
    Ga%gtr(0:,0:)  = Ga%gtr(0:,0:)  + Gb%gtr(0:,0:)*C
    Ga%lmix(0:,0:) = Ga%lmix(0:,0:) + Gb%lmix(0:,0:)*C
    Ga%gmix(0:,0:) = Ga%gmix(0:,0:) + Gb%gmix(0:,0:)*C
    Ga%mats(0:,0:) = Ga%mats(0:,0:) + Gb%mats(0:,0:)*C
  end subroutine kbm_contour_gf_sum_z

  !******************************************************************
  !******************************************************************
  !******************************************************************

  function kbm_contour_gf_scalarL_d(C,G) result(F)
    real(8),intent(in) :: C
    type(kbm_contour_gf),intent(in) :: G
    type(kbm_contour_gf) :: F
    F%less(0:,0:)= C*G%less(0:,0:)
    F%gtr(0:,0:) = C*G%gtr(0:,0:)
    F%lmix(0:,0:)= C*G%lmix(0:,0:)
    F%gmix(0:,0:)= C*G%gmix(0:,0:)
    F%mats(0:,0:)= C*G%mats(0:,0:)
  end function kbm_contour_gf_scalarL_d

  !******************************************************************
  !******************************************************************
  !******************************************************************

  function kbm_contour_gf_scalarL_c(C,G) result(F)
    complex(8),intent(in) :: C
    type(kbm_contour_gf),intent(in) :: G
    type(kbm_contour_gf) :: F
    F%less(0:,0:)=C*G%less(0:,0:)
    F%gtr(0:,0:)=C*G%gtr(0:,0:)
    F%lmix(0:,0:)=C*G%lmix(0:,0:)
    F%gmix(0:,0:)=C*G%gmix(0:,0:)
    F%mats(0:,0:)=C*G%mats(0:,0:)
  end function kbm_contour_gf_scalarL_c

  !******************************************************************
  !******************************************************************
  !******************************************************************

  function kbm_contour_gf_scalarR_d(G,C) result(F)
    real(8),intent(in) :: C
    type(kbm_contour_gf),intent(in) :: G
    type(kbm_contour_gf) :: F
    F%less(0:,0:)=G%less(0:,0:)*C
    F%gtr(0:,0:)=G%gtr(0:,0:)*C
    F%lmix(0:,0:)=G%lmix(0:,0:)*C
    F%gmix(0:,0:)=G%gmix(0:,0:)*C
    F%mats(0:,0:)=G%mats(0:,0:)*C
  end function kbm_contour_gf_scalarR_d

  !******************************************************************
  !******************************************************************
  !******************************************************************

  function kbm_contour_gf_scalarR_c(G,C) result(F)
    complex(8),intent(in) :: C
    type(kbm_contour_gf),intent(in) :: G
    type(kbm_contour_gf) :: F
    F%less(0:,0:)=G%less(0:,0:)*C
    F%gtr(0:,0:)=G%gtr(0:,0:)*C
    F%lmix(0:,0:)=G%lmix(0:,0:)*C
    F%gmix(0:,0:)=G%gmix(0:,0:)*C
    F%mats(0:,0:)=G%mats(0:,0:)*C  
  end function kbm_contour_gf_scalarR_c

  !******************************************************************
  !******************************************************************
  !******************************************************************

  !This a rough implementation: just a shortcut:
  subroutine mpi_reduce_kbm_contour_gf(tmpG,G)
    type(kbm_contour_gf) :: tmpG
    type(kbm_contour_gf) :: G
    if( .not.g%status )call error("ERROR contour_gf/mpi_reduce_kbm_contour_gf: object function not allocated.")
    if( (G%N /= tmpG%N) .OR. (G%L /=tmpG%L) )call error("ERROR contour_gf/mpi_reduce_kbm_contour_gf: wrong dimensions.")
    call MPI_REDUCE(tmpG%less(0:,0:),G%less(0:,0:),size(tmpG%less(0:,0:)),&
         MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_REDUCE(tmpG%gtr(0:,0:),G%gtr(0:,0:),size(tmpG%gtr(0:,0:)),&
         MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_REDUCE(tmpG%lmix(0:,0:),G%lmix(0:,0:),size(tmpG%lmix(0:,0:)),&
         MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_REDUCE(tmpG%gmix(0:,0:),G%gmix(0:,0:),size(tmpG%gmix(0:,0:)),&
         MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_REDUCE(tmpG%mats(0:,0:),G%mats(0:,0:),size(tmpG%mats(0:,0:)),&
         MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
  end subroutine mpi_reduce_kbm_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine mpi_bcast_kbm_contour_gf(G)
    type(kbm_contour_gf),intent(inout) :: G
    if( .not.g%status )call error("ERROR contour_gf/mpi_bcast_kbm_contour_gf: object function not allocated.")    
    call MPI_BCAST(G%less(0:,0:),size(G%less(0:,0:)),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(G%gtr(0:,0:),size(G%gtr(0:,0:)),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(G%lmix(0:,0:),size(G%lmix(0:,0:)),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(G%gmix(0:,0:),size(G%gmix(0:,0:)),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(G%mats(0:,0:),size(G%mats(0:,0:)),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
  end subroutine mpi_bcast_kbm_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************


END MODULE CONTOUR_GF
