MODULE CONTOUR_GF
  USE COMMON_VARS
  USE IOTOOLS
  implicit none
  private

  type :: kbm_contour_gf
     complex(8),dimension(:,:),pointer  :: less,gtr
     complex(8),dimension(:,:),pointer  :: lmix,gmix
     real(8),dimension(:,:),pointer     :: mats
     logical                            :: status=.false.
     integer                            :: N=0,L=0
  end type kbm_contour_gf

  interface assignment(=)
     module procedure kbm_contour_gf_equality,kbm_contour_gf_equality_
  end interface assignment(=)

  interface operator(*)
     module procedure kbm_contour_gf_scalarL_d,kbm_contour_gf_scalarL_c,&
          kbm_contour_gf_scalarR_d,kbm_contour_gf_scalarR_c
  end interface operator(*)

  public :: kbm_contour_gf
  public :: allocate_kbm_contour_gf
  public :: deallocate_kbm_contour_gf
  public :: write_kbm_contour_gf,read_kbm_contour_gf
  public :: assignment(=)
  public :: operator(*)

contains

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

  subroutine deallocate_kbm_contour_gf(G)
    type(kbm_contour_gf) :: G
    deallocate(G%less,G%gtr,G%lmix,G%gmix,G%mats)
    G%N=0
    G%L=0
    G%status=.false.
  end subroutine deallocate_kbm_contour_gf


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
    call splot(trim(file)//"_less.ipt",G%less(0:,0:))
    call splot(trim(file)//"_gtr.ipt",G%gtr(0:,0:))
    call splot(trim(file)//"_lmix.ipt",G%lmix(0:,0:))
    call splot(trim(file)//"_gmix.ipt",G%gmix(0:,0:))
    call splot(trim(file)//"_mats.ipt",G%mats(0:,0:))
  end subroutine write_kbm_contour_gf

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
    call sread(trim(file)//"_less.ipt",G%less(0:,0:))
    call sread(trim(file)//"_gtr.ipt",G%gtr(0:,0:))
    call sread(trim(file)//"_lmix.ipt",G%lmix(0:,0:))
    call sread(trim(file)//"_gmix.ipt",G%gmix(0:,0:))
    call sread(trim(file)//"_mats.ipt",G%mats(0:,0:))
  end subroutine read_kbm_contour_gf

  subroutine kbm_contour_gf_equality(G1,G2)
    type(kbm_contour_gf),intent(inout) :: G1
    type(kbm_contour_gf),intent(in)    :: G2
    G1%less = G2%less
    G1%gtr = G2%gtr
    G1%lmix = G2%lmix
    G1%gmix = G2%gmix
    G1%mats = G2%mats
  end subroutine kbm_contour_gf_equality

  subroutine kbm_contour_gf_equality_(G1,C)
    type(kbm_contour_gf),intent(inout) :: G1
    complex(8),intent(in) :: C
    G1%less = C
    G1%gtr = C
    G1%lmix = C
    G1%gmix = C
    G1%mats = C
  end subroutine kbm_contour_gf_equality_

  function kbm_contour_gf_scalarL_d(C,G)
    real(8),intent(in) :: C
    type(kbm_contour_gf),intent(in) :: G
    type(kbm_contour_gf) :: kbm_contour_gf_scalarL_d
    kbm_contour_gf_scalarL_d%less=C*G%less
    kbm_contour_gf_scalarL_d%gtr=C*G%gtr
    kbm_contour_gf_scalarL_d%lmix=C*G%lmix
    kbm_contour_gf_scalarL_d%gmix=C*G%gmix
    kbm_contour_gf_scalarL_d%mats=C*G%mats    
  end function kbm_contour_gf_scalarL_d
  function kbm_contour_gf_scalarL_c(C,G)
    complex(8),intent(in) :: C
    type(kbm_contour_gf),intent(in) :: G
    type(kbm_contour_gf) :: kbm_contour_gf_scalarL_c
    kbm_contour_gf_scalarL_c%less=C*G%less
    kbm_contour_gf_scalarL_c%gtr=C*G%gtr
    kbm_contour_gf_scalarL_c%lmix=C*G%lmix
    kbm_contour_gf_scalarL_c%gmix=C*G%gmix
    kbm_contour_gf_scalarL_c%mats=C*G%mats    
  end function kbm_contour_gf_scalarL_c

  function kbm_contour_gf_scalarR_d(G,C)
    real(8),intent(in) :: C
    type(kbm_contour_gf),intent(in) :: G
    type(kbm_contour_gf) :: kbm_contour_gf_scalarR_d
    kbm_contour_gf_scalarR_d%less=G%less*C
    kbm_contour_gf_scalarR_d%gtr=G%gtr*C
    kbm_contour_gf_scalarR_d%lmix=G%lmix*C
    kbm_contour_gf_scalarR_d%gmix=G%gmix*C
    kbm_contour_gf_scalarR_d%mats=G%mats*C
  end function kbm_contour_gf_scalarR_d
  function kbm_contour_gf_scalarR_c(G,C)
    complex(8),intent(in) :: C
    type(kbm_contour_gf),intent(in) :: G
    type(kbm_contour_gf) :: kbm_contour_gf_scalarR_c
    kbm_contour_gf_scalarR_c%less=G%less*C
    kbm_contour_gf_scalarR_c%gtr=G%gtr*C
    kbm_contour_gf_scalarR_c%lmix=G%lmix*C
    kbm_contour_gf_scalarR_c%gmix=G%gmix*C
    kbm_contour_gf_scalarR_c%mats=G%mats*C  
  end function kbm_contour_gf_scalarR_c

END MODULE CONTOUR_GF
