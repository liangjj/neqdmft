MODULE CONTOUR_GF
  USE COMMON_VARS
  USE TOOLS
  USE IOTOOLS
  USE SPLINE
  USE MPI
  implicit none
  private


  !##################################################################
  ! KELDYSH CONTOUR GREEN'S FUNCTIONS
  !##################################################################
  type :: keldysh_contour_gf
     complex(8),dimension(:,:),pointer  :: less,gtr
     logical                            :: status
     integer                            :: N
  end type keldysh_contour_gf

  public :: keldysh_contour_gf
  public :: allocate_keldysh_contour_gf
  public :: deallocate_keldysh_contour_gf
  public :: write_keldysh_contour_gf
  public :: read_keldysh_contour_gf
  public :: plot_keldysh_contour_gf
  public :: inquire_keldysh_contour_gf

  !##################################################################
  ! KELDYSH-BAYM-MATSUBARS CONTOUR GREEN'S FUNCTIONS:
  !##################################################################
  type :: kbm_contour_gf
     complex(8),dimension(:,:),pointer  :: less,gtr
     complex(8),dimension(:,:),pointer  :: lmix,gmix
     real(8),dimension(:,:),pointer     :: mats
     logical                            :: status=.false.
     integer                            :: N=0,L=0
  end type kbm_contour_gf


  public :: kbm_contour_gf
  public :: allocate_kbm_contour_gf
  public :: deallocate_kbm_contour_gf
  public :: write_kbm_contour_gf
  public :: inquire_kbm_contour_gf
  public :: read_kbm_contour_gf
  public :: plot_kbm_contour_gf


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
          kbm_contour_gf_equality,&
          kbm_contour_gf_equality_
  end interface assignment(=)
  public :: assignment(=)
  public :: operator(*)


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

  function inquire_keldysh_contour_gf(file) result(check)
    logical          :: check,bool1,bool2
    character(len=*) :: file
    inquire(file=trim(file)//"_less.data",exist=bool1)
    if(.not.bool1)inquire(file=trim(file)//"_less.data.gz",exist=bool1)
    if(.not.bool1)call warning("Can not read "//trim(file)//"_less.data")
    inquire(file=trim(file)//"_gtr.data",exist=bool2)
    if(.not.bool2)inquire(file=trim(file)//"_gtr.data.gz",exist=bool2)
    if(.not.bool2)call warning("Can not read "//trim(file)//"_gtr.data")
    check=bool1.AND.bool2
  end function inquire_keldysh_contour_gf

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
    call splot(trim(file)//"_less_t_t",t(0:),t(0:),G%less(0:,0:))
    call splot(trim(file)//"_gtr_t_t",t(0:),t(0:),G%gtr(0:,0:))
  end subroutine plot_keldysh_contour_gf

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

  function inquire_kbm_contour_gf(file) result(check)
    integer          :: i
    logical          :: check,bool(5)
    character(len=*) :: file
    character(len=16),dimension(5)  :: ctype=(['less','gtr','lmix','gmix','mats'])
    check=.true.
    do i=1,5
       inquire(file=trim(file)//"_"//trim(ctype(i))//".data",exist=bool(i))
       if(.not.bool(i))inquire(file=trim(file)//"_"//trim(ctype(i))//".data.gz",exist=bool(i))
       !if(.not.bool(i))call warning("Can not read "//trim(file)//"_"//trim(ctype(i))//".data")
       check=check.AND.bool(i)
    enddo
  end function inquire_kbm_contour_gf


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
    integer               :: i,j,N,L,M
    !real(8),allocatable   :: ftau(:),upmGtau(:),Gtau(:),Gmats(:,:)
    N=G%N+1 ; L=G%L+1
    if( (size(G%less)/=N**2) .OR. (size(G%gtr)/=N**2) )&
         call error("ERROR contour_gf/plot_kbm_contour_gf: wrong dimensions 1")
    if( (size(G%lmix)/=N*L)  .OR. (size(G%gmix)/=N*L) )&
         call error("ERROR contour_gf/plot_kbm_contour_gf: wrong dimensions 2")
    if( size(G%mats)/=L*L)&
         call error("ERROR contour_gf/plot_kbm_contour_gf: wrong dimensions 3")
    call splot(trim(file)//"_less_t_t",t(0:),t(0:),G%less(0:,0:))
    call splot(trim(file)//"_gtr_t_t",t(0:),t(0:),G%gtr(0:,0:))
    call splot(trim(file)//"_lmix_t_tau",t(0:),tau(0:),G%lmix(0:,0:))
    call splot(trim(file)//"_gmix_tau_t",tau(0:),t(0:),G%gmix(0:,0:))
    call splot(trim(file)//"_mats_tau_tau",tau(0:),tau(0:),G%mats(0:,0:))
    ! L=G%L
    ! M=min(2*L,200)
    ! allocate(ftau(0:M),upmGtau(0:L),Gtau(-M:M),Gmats(0:M,0:M))
    ! ftau(0:) = linspace(minval(tau(0:)),maxval(tau(0:)),M+1)
    ! forall(i=0:L)upmGtau(i) = G%mats(i,0)
    ! call cubic_spline(upmGtau(0:L),tau(0:L),Gtau(0:M),ftau(0:M))
    ! forall(i=1:M)Gtau(-i)=-Gtau(M-i)
    ! forall(i=0:M,j=0:M)Gmats(i,j)=Gtau(i-j)
    ! call splot(trim(file)//"_mats_tau_tau",ftau(0:),ftau(0:),Gmats(0:,0:))
    ! deallocate(ftau,upmGtau,Gtau,Gmats)
  end subroutine plot_kbm_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine kbm_contour_gf_equality(G1,C)
    type(kbm_contour_gf),intent(inout) :: G1
    complex(8),intent(in) :: C
    G1%less(0:,0:) = C
    G1%gtr(0:,0:) = C
    G1%lmix(0:,0:) = C
    G1%gmix(0:,0:) = C
    G1%mats(0:,0:) = C
  end subroutine kbm_contour_gf_equality

  subroutine kbm_contour_gf_equality_(G1,G2)
    type(kbm_contour_gf),intent(inout) :: G1
    type(kbm_contour_gf),intent(in)    :: G2
    if(G1%N/=G2%N .OR. G1%L/=G2%L)&
         call error("ERROR contour_gf/kbm_contour_equality_: wrong dimensions")
    G1%less(0:,0:)=G2%less(0:,0:)
    G1%gtr(0:,0:)=G2%gtr(0:,0:)
    G1%lmix(0:,0:)=G2%lmix(0:,0:)
    G1%gmix(0:,0:)=G2%gmix(0:,0:)
    G1%mats(0:,0:)=G2%mats(0:,0:)
  end subroutine kbm_contour_gf_equality_

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



END MODULE CONTOUR_GF
