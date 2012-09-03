!###############################################################
!     PURPOSE  : A non-equilibrium IPT solver module. 
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_NEQ
  USE VARS_GLOBAL
  USE EQUILIBRIUM
  USE ELECTRIC_FIELD
  USE MATRIX
  implicit none
  private

  public  :: neq_init_run
  public  :: neq_solve_ipt


contains


  subroutine neq_init_run()
    integer                          :: i,j,ik
    integer                          :: iselect,irdL
    real(8)                          :: en,intE,A,irdfmesh
    real(8)                          :: nless,ngtr,xmu_,beta_
    complex(8)                       :: peso
    logical                          :: checkS,checkS1,checkS2
    logical                          :: checkGN,checkG0,checkNk
    real(8),allocatable,dimension(:) :: irdwr
    type(matsubara_gf)               :: fg0m,sm

    !Initial selection: no file exist
    iselect=0

    !Check if Sigma^<,> files exist:
    inquire(file=trim(irdSlfile),exist=checkS1)
    if(.not.checkS1)inquire(file=trim(irdSlfile)//".gz",exist=checkS1)
    inquire(file=trim(irdSgfile),exist=checkS2)
    if(.not.checkS2)inquire(file=trim(irdSgfile)//".gz",exist=checkS2)

    !Check if G0(w) file exists:
    inquire(file=trim(irdG0file),exist=checkG0)
    if(.not.checkG0)inquire(file=trim(irdG0file)//".gz",exist=checkG0)

    !Check if n(k) file exists:
    inquire(file=trim(irdnkfile),exist=checkNk)
    if(.not.checkNk)inquire(file=trim(irdNkfile)//".gz",exist=checkNk)

    checkS=checkS1.AND.checkS2.AND.checkNk

    checkGN=checkG0.AND.checkNk

    if(checkS)then              !S^<,>(t,t') AND n(k) files exist
       iselect=1
    elseif(checkGN)then         !G0(w) AND n(k) files exist
       iselect=2
    endif



    select case(iselect)

    case default                !No files are given:

       if(mpiID==0)then
          !Get non-interacting n(k):
          xmu_=xmu   ; if(iquench)xmu_=xmu0
          beta_=beta ; if(iquench)beta_=beta0
          do ik=1,Lk
             irdNk(ik)=fermi0((epsik(ik)-xmu_),beta_)
          enddo
          call splot("guessnkVSepsk.ipt",epsik,irdnk)

          !Guess G0-->Sigma^(0)
          if(g0loc_guess)then
             if(equench)then
                call msg("Using G0_loc + electric field for self-energy guess",lines=1)
                do ik=1,Lk
                   en   = epsik(ik)
                   nless= fermi0(en,beta)
                   ngtr = fermi0(en,beta)-1.d0
                   do j=0,nstep
                      do i=0,nstep
                         intE=int_Ht(ik,i,j)
                         peso=exp(-xi*intE)
                         G0%less(i,j)= G0%less(i,j) + xi*nless*peso*wt(ik)
                         G0%gtr(i,j) = G0%gtr(i,j)  + xi*ngtr*peso*wt(ik)
                      enddo
                   enddo
                enddo
             else
                call msg("Using G0_loc for self-energy guess",lines=1)
                do ik=1,Lk
                   en   = epsik(ik)
                   nless= fermi0(en,beta)
                   ngtr = fermi0(en,beta)-1.d0
                   A    = wt(ik)
                   do i=-nstep,nstep
                      peso=exp(-xi*en*t(i))
                      gf0%less%t(i)=gf0%less%t(i) + xi*nless*A*peso
                      gf0%gtr%t(i) =gf0%gtr%t(i)  + xi*ngtr*A*peso
                   enddo
                enddo
                forall(i=0:nstep,j=0:nstep)
                   G0%less(i,j)=gf0%less%t(i-j)
                   G0%gtr(i,j) =gf0%gtr%t(i-j)
                end forall
             endif

          else

             call msg("Using Hartree-Fock for self-energy guess")
             call msg("G0less=G0gtr=zero",lines=1)
             G0=zero

          endif

          call splot("guessG0less.data",G0%less(0:nstep,0:nstep))
          call splot("guessG0gtr.data",G0%gtr(0:nstep,0:nstep))
       endif


       call MPI_BCAST(G0%less,(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
       call MPI_BCAST(G0%gtr,(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
       call MPI_BCAST(irdNk,Lk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
       call neq_solve_ipt()


    case(1)
       call msg("Reading self-energy guess and n(k) from input files",lines=1)
       call sread(trim(irdSlfile),Sig%less(0:nstep,0:nstep))
       call sread(trim(irdSgfile),Sig%gtr(0:nstep,0:nstep))
       call read_nkfile(trim(irdnkfile))


    case(2)
       call msg("Reading G0(w) and n(k) from input files")
       call msg("Using G0(w) for self-energy guess",lines=1)
       if(mpiID==0)then
          call read_nkfile(trim(irdnkfile))
          irdL=file_length(trim(irdG0file))
          allocate(irdG0w(irdL),irdwr(irdL))
          call sread(trim(irdG0file),irdwr,irdG0w)
          !1)
          ! call linear_spline(irdG0w,irdwr,gf0%ret%w,wr)
          ! gf0%less%w = less_component_w(gf0%ret%w,wr,beta)
          ! gf0%gtr%w  = gtr_component_w(gf0%ret%w,wr,beta)
          ! call fftgf_rw2rt(gf0%less%w,gf0%less%t,nstep) ; gf0%less%t=fmesh/pi2*gf0%less%t
          ! call fftgf_rw2rt(gf0%gtr%w, gf0%gtr%t,nstep)  ; gf0%gtr%t =fmesh/pi2*gf0%gtr%t
          !2)
          irdfmesh=abs(irdwr(2)-irdwr(1)) !Get G0 mesh:
          do ik=1,irdL
             en   = irdwr(ik)
             nless= fermi0(en,beta)
             ngtr = fermi0(en,beta)-1.d0
             A    = -aimag(irdG0w(ik))/pi*irdfmesh
             do i=-nstep,nstep
                peso=exp(-xi*en*t(i))
                gf0%less%t(i)=gf0%less%t(i) + xi*nless*A*peso
                gf0%gtr%t(i) =gf0%gtr%t(i)  + xi*ngtr*A*peso
             enddo
          enddo
          forall(i=0:nstep,j=0:nstep)
             G0%less(i,j)=gf0%less%t(i-j)
             G0%gtr(i,j) =gf0%gtr%t(i-j)
          end forall
          call splot("guessG0less.data",G0%less(0:nstep,0:nstep))
          call splot("guessG0gtr.data",G0%gtr(0:nstep,0:nstep))
       endif
       call MPI_BCAST(G0%less,(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
       call MPI_BCAST(G0%gtr,(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
       call MPI_BCAST(irdNk,Lk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)

       call neq_solve_ipt()


       ! case(3)
       !    call msg("Reading G0(w) from input file")
       !    call msg("Using G0(w) to build n(k) and for self-energy guess",lines=1)
       !    if(mpiID==0)then
       !       irdL=file_length(trim(irdG0file))
       !       allocate(irdG0w(irdL),irdwr(irdL))
       !       call sread(trim(irdG0file),irdwr,irdG0w)
       !       !1)
       !       ! call linear_spline(irdG0w,irdwr,gf0%ret%w,wr)
       !       ! gf0%less%w = less_component_w(gf0%ret%w,wr,beta)
       !       ! gf0%gtr%w  = gtr_component_w(gf0%ret%w,wr,beta)
       !       ! call fftgf_rw2rt(gf0%less%w,gf0%less%t,nstep) ; gf0%less%t=fmesh/pi2*gf0%less%t
       !       ! call fftgf_rw2rt(gf0%gtr%w, gf0%gtr%t,nstep)  ; gf0%gtr%t =fmesh/pi2*gf0%gtr%t
       !       !2)
       !       irdfmesh=abs(irdwr(2)-irdwr(1)) !Get G0 mesh:
       !       do ik=1,irdL
       !          en   = irdwr(ik)
       !          nless= fermi0(en,beta)
       !          ngtr = fermi0(en,beta)-1.d0
       !          A    = -aimag(irdG0w(ik))/pi*irdfmesh
       !          do i=-nstep,nstep
       !             peso=exp(-xi*en*t(i))
       !             gf0%less%t(i)=gf0%less%t(i) + xi*nless*A*peso
       !             gf0%gtr%t(i) =gf0%gtr%t(i)  + xi*ngtr*A*peso
       !          enddo
       !       enddo
       !       forall(i=0:nstep,j=0:nstep)
       !          G0%less(i,j)=gf0%less%t(i-j)
       !          G0%gtr(i,j) =gf0%gtr%t(i-j)
       !       end forall
       !       call splot("guessG0less.data",G0%less(0:nstep,0:nstep))
       !       call splot("guessG0gtr.data",G0%gtr(0:nstep,0:nstep))


       !       !Get n(k) within IPT approximation!!
       !       call msg("Getting n(k) within using IPT approximation!!")
       !       call allocate_gf(fg0m,L)
       !       call allocate_gf(sm,L)
       !       call get_matsubara_gf_from_DOS(irdwr,irdG0w,fg0m%iw,beta)
       !       call fftgf_iw2tau(fg0m%iw,fg0m%tau,beta)
       !       forall(i=0:L)sm%tau(i)=U**2*(fg0m%tau(i))**2*fg0m%tau(L-i)
       !       call fftgf_tau2iw(sm%tau,sm%iw,beta)
       !       do ik=1,Lk
       !          fg0m%iw=one/(xi*wm - epsik(ik) - sm%iw)
       !          call fftgf_iw2tau(fg0m%iw,fg0m%tau,beta)
       !          irdnk(ik)=-fg0m%tau(L)
       !       enddo
       !       call splot("guessnkVSepsk.ipt",epsik,irdnk)

       !    endif
       !    call MPI_BCAST(G0%less,(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
       !    call MPI_BCAST(G0%gtr,(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
       !    call MPI_BCAST(irdNk,Lk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
       !    call neq_solve_ipt()

    end select


  contains

    function int_Ht(ik,it,jt)
      real(8)      :: int_Ht
      integer      :: i,j,ii,ik,it,jt,sgn
      type(vect2D) :: kt,Ak
      int_Ht=0.d0 ; if(it==jt)return
      sgn=1 ; if(jt > it)sgn=-1
      i=ik2ix(ik); j=ik2iy(ik)
      do ii=jt,it,sgn
         Ak=Afield(t(ii),Ek)
         kt=kgrid(i,j) - Ak
         int_Ht=int_Ht + sgn*square_lattice_dispersion(kt)*dt
      enddo
    end function int_Ht

    subroutine read_nkfile(file)
      character(len=*)    :: file
      integer             :: redLk
      real(8),allocatable :: rednk(:),redek(:)
      integer,allocatable :: orderk(:)
      real(8),allocatable :: uniq_rednk(:),uniq_redek(:)
      logical,allocatable :: maskk(:)
      !n(k): A lot of work here to reshape the array
      redLk=file_length(file)
      allocate(rednk(redLk),redek(redLk),orderk(redLk))
      call sread(file,redek,rednk)
      !
      !work on the read arrays:
      !1 - sorting: sort the energies (X-axis), mirror on occupation (Y-axis) 
      !2 - delete duplicates energies (X-axis), mirror on occupation (Y-axis) 
      !3 - interpolate to the actual lattice structure (epsik,nk)
      call sort_array(redek,orderk)
      call reshuffle(rednk,orderk)
      call uniq(redek,uniq_redek,maskk)
      allocate(uniq_rednk(size(uniq_redek)))
      uniq_rednk = pack(rednk,maskk)
      call linear_spline(uniq_rednk,uniq_redek,irdnk,epsik)
    end subroutine read_nkfile

  end subroutine neq_init_run


  !+-------------------------------------------------------------------+
  !PURPOSE  : BUild the 2^nd IPT sigma functions:
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt()
    integer                               :: i,j,itau
    !Get SIgma:
    call msg("Get Sigma(t,t')")

    Sig=zero
    forall(i=0:nstep,j=0:nstep)
       Sig%less(i,j) = (U**2)*(G0%less(i,j)**2)*G0%gtr(j,i)
       Sig%gtr (i,j) = (U**2)*(G0%gtr(i,j)**2)*G0%less(j,i)
    end forall

    ! !Get impurity GF and use SPT method if required
    ! call get_impuritygf
    ! if(method=="spt")then
    !    call msg("Recalculate Sigma using SPT")
    !    forall(i=0:nstep,j=0:nstep)
    !       Sig%less(i,j) = (U**2)*(impG%less(i,j)**2)*impG%gtr(j,i)
    !       Sig%gtr (i,j) = (U**2)*(impG%gtr(i,j)**2)*impG%less(j,i)                
    !    end forall
    ! end if

    !Save data:
    if(mpiID==0)then
       call splot("Sless.data",Sig%less(0:nstep,0:nstep))
       call splot("Sgtr.data",Sig%gtr(0:nstep,0:nstep))
    endif

  end subroutine Neq_solve_ipt



  !********************************************************************
  !********************************************************************
  !********************************************************************




  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : evaluate the impurity neq Green's functions
  ! !+-------------------------------------------------------------------+
  ! subroutine get_impuritygf()
  !   integer                               :: i,j
  !   real(8)                               :: A,w
  !   complex(8),dimension(0:nstep,0:nstep) :: Uno,GammaRet,Gamma0Ret
  !   complex(8),dimension(0:nstep,0:nstep) :: dG0ret,dGret,dSret

  !   if(update_wfftw)then
  !      call get_equilibrium_impuritygf !not tested!
  !   else
  !      include "obtain_Gimp_nonequilibrium.f90"
  !   endif

  !   !Save data:
  !   if(mpiID==0)then
  !      call splot("impGless.data",impG%less(0:nstep,0:nstep))
  !      call splot("impGgtr.data",impG%gtr(0:nstep,0:nstep))
  !   endif

  ! end subroutine get_impuritygf





  !********************************************************************
  !********************************************************************
  !********************************************************************


end module IPT_NEQ
