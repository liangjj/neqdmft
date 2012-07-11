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

  !+-------------------------------------------------------------------+
  !PURPOSE  : Initialize the run guessing/reading/setting initial conditions
  !+-------------------------------------------------------------------+
  subroutine neq_init_run()
    integer                          :: i,j,ik
    integer                          :: iselect,irdL
    real(8)                          :: en,A,fmesh_
    real(8)                          :: nless,ngtr,xmu_,beta_
    complex(8)                       :: peso
    logical                          :: checkG0,checkNk,checkSm,checkEQ
    real(8),allocatable,dimension(:) :: wr_,wm_

    call system("if [ ! -d GUESS ]; then mkdir GUESS; fi")

    !Initial selection: no file exist
    iselect=0

    !Check if n(k) file exist:
    inquire(file=trim(irdnkfile),exist=checkNk)
    if(.not.checkNk)inquire(file=trim(irdNkfile)//".gz",exist=checkNk)
    inquire(file=trim(irdSiwfile),exist=checkSm)
    if(.not.checkSm)inquire(file=trim(irdSiwfile)//".gz",exist=checkSm)
    inquire(file=trim(irdG0wfile),exist=checkG0)
    if(.not.checkG0)inquire(file=trim(irdG0wfile)//".gz",exist=checkG0)
    checkEQ=checkNk.AND.checkSm.AND.checkG0
    if(checkEQ)iselect=1

    select case(iselect)        !add more solutions here (e.g. restart from read file)
    case(-1)
       write(*,*)"Iselect=",iselect
       call error("A problem in neq_init_run!")

    case default
       call msg("Temporary use of non-interacting GUESS!! ONLY if U=0",id=0)
       !Get non-interacting n(k):
       !Initial conditions:
       allocate(eq_nk(Lk))
       allocate(eq_G0w(L))
       allocate(eq_Siw(L))
       xmu_=xmu   ; if(iquench)xmu_=xmu0
       beta_=beta ; if(iquench)beta_=beta0
       do ik=1,Lk
          eq_nk(ik)=fermi0((epsik(ik)-xmu_),beta_)
       enddo
       if(mpiID==0)call splot("GUESS/guessnkVSepsk.ipt",epsik,eq_nk)
       Sgtr=zero   ; Sless=zero !to be changed when-out-half-filling
       Sgmix=zero  ; Slmix=zero
       Smat=zero
       eq_Siw=zero
       eq_G0w=zero

       G0less=zero
       G0lmix=zero
       do ik=1,Lk
          en   = epsik(ik)
          nless= fermi0(en,beta)
          ngtr = fermi0(en,beta)-1.d0
          do i=0,nstep
             do j=0,nstep
                peso=exp(-xi*en*(t(i)-t(j)))
                G0less(i,j)=G0less(i,j) + xi*nless*wt(ik)*peso
             enddo
             do j=0,Ltau
                peso=exp(-xi*en*t(i))*exp(-en*tau(j))
                G0lmix(i,j)=G0lmix(i,j) + xi*nless*wt(ik)*peso
             enddo
          enddo
       enddo
       call splot("G0less_t_t.ipt",t(0:),t(0:),G0less(0:,0:))
       call splot("G0lmix_t_tau.ipt",t(0:),tau(0:),G0lmix(0:,0:))

    case(1)
       call msg("Continuing the EQUILIBRIUM SOLUTION")
       !READ n(k) --> eq_nk(:)
       allocate(eq_nk(Lk))
       call read_nkfile(trim(irdnkfile))

       !READ \Sigma(iw) --> eq_Siw(:)
       irdL=file_length(trim(irdSiwfile))
       allocate(eq_Siw(irdL),wm_(irdL))
       call sread(trim(irdSiwfile),wm_,eq_Siw)
       deallocate(wm_)

       !READ \GG_0(w) --> eq_G0w(:)
       irdL=file_length(trim(irdG0wfile))
       allocate(eq_G0w(irdL))
       allocate(wr_(irdL))
       call sread(trim(irdG0wfile),wr_,eq_G0w)
       fmesh_=abs(wr_(2)-wr_(1)) !Get G0 mesh:
       G0lmix=zero;gf0=zero
       do ik=1,irdL
          en   = wr_(ik)
          nless= fermi0(en,beta)
          ngtr = fermi0(en,beta)-1.d0
          A    = -aimag(eq_G0w(ik))/pi*fmesh_
          do i=-nstep,nstep
             peso=exp(-xi*en*t(i))
             gf0%less%t(i)=gf0%less%t(i) + xi*nless*A*peso
             gf0%gtr%t(i) =gf0%gtr%t(i)  + xi*ngtr*A*peso
          enddo
          do i=0,nstep
             do j=0,Ltau
                peso=exp(-xi*en*t(i))*exp(-en*tau(j))
                G0lmix(i,j)=G0lmix(i,j) + xi*nless*A*peso
             enddo
          enddo
       enddo
       forall(i=0:nstep,j=0:nstep)
          G0less(i,j)=gf0%less%t(i-j)
          G0gtr(i,j) =gf0%gtr%t(i-j)
       end forall
       forall(j=0:Ltau)G0gmix(j,:)=-conjg(G0lmix(:,Ltau-j))
       deallocate(wr_)

       if(mpiID==0)then
          call splot("guessG0less.data",G0less(0:nstep,0:nstep))
          call splot("guessG0gtr.data",G0gtr(0:nstep,0:nstep))
          call splot("guessG0lmix.data",G0lmix(0:nstep,0:Ltau))
          call splot("guessG0less_t_t.ipt",t(0:),t(0:),G0less(0:,0:))
          call splot("guessG0lmix_t_tau.ipt",t(0:),tau(0:),G0lmix(0:,0:))
          call system("mv -vf *guess*.ipt GUESS/")
       endif
       call neq_solve_ipt()

    end select


  contains

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
      !work on the read arrays:
      !1 - sorting: sort the energies (X-axis), mirror on occupation (Y-axis) 
      !2 - delete duplicates energies (X-axis), mirror on occupation (Y-axis) 
      !3 - interpolate to the actual lattice structure (epsik,nk)
      call sort_array(redek,orderk)
      call reshuffle(rednk,orderk)
      call uniq(redek,uniq_redek,maskk)
      allocate(uniq_rednk(size(uniq_redek)))
      uniq_rednk = pack(rednk,maskk)
      call linear_spline(uniq_rednk,uniq_redek,eq_nk,epsik)
    end subroutine read_nkfile

  end subroutine neq_init_run




  !+-------------------------------------------------------------------+
  !PURPOSE  : BUild the 2^nd IPT sigma functions:
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt()
    integer      :: i,j,itau
    integer,save :: loop=1

    !Get SIgma:
    call msg("Get Sigma(t,t')")
    call system("if [ ! -d SIGMA ]; then mkdir SIGMA; fi")

    forall(i=0:nstep,j=0:nstep)
       Sgtr (i,j) = (U**2)*(G0gtr(i,j)**2)*G0less(j,i)
       Sless(i,j) = (U**2)*(G0less(i,j)**2)*G0gtr(j,i)
    end forall

    forall(i=0:nstep,itau=0:Ltau)Slmix(i,itau) = -U**2*(G0lmix(i,itau)**2)*G0gmix(itau,i)
    forall(j=0:Ltau)Sgmix(j,:)=-conjg(Slmix(:,Ltau-j))

    !Save data:
    if(mpiID==0)then
       call splot("Sless.data",Sless(0:nstep,0:nstep))
       call splot("Sgtr.data",Sgtr(0:nstep,0:nstep))
       call splot("Slmix.data",Slmix(0:nstep,0:Ltau))
       call splot("Slmix_t_tau.ipt",t(0:nstep),tau(0:Ltau),Slmix(0:nstep,0:Ltau))
       call splot("Sless_t_t.ipt",t(0:),t(0:),Sless(0:,0:))
       call splot("Slmix_t_tau.ipt",t(0:),tau(0:),Slmix(0:,0:))
       call system("mv -vf *Sl*_t*.ipt SIGMA/")
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
  !      call splot("impGless.data",impGless(0:nstep,0:nstep))
  !      call splot("impGgtr.data",impGgtr(0:nstep,0:nstep))
  !   endif
  ! end subroutine get_impuritygf





  !********************************************************************
  !********************************************************************
  !********************************************************************


end module IPT_NEQ
