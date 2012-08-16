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
    integer                             :: i,j,ik
    real(8)                             :: en,A,fmesh_,imt
    real(8)                             :: nless,ngtr,xmu_,beta_
    complex(8)                          :: peso
    logical                             :: checkG0
    real(8),allocatable,dimension(:)    :: wr_
    real(8),dimension(-Ltau:Ltau)       :: gmtau
    real(8),allocatable,dimension(:)    :: dummyGtau,Stau
    complex(8),allocatable,dimension(:) :: dummyGiw


    call system("if [ ! -d GUESS ]; then mkdir GUESS; fi")

    inquire(file=trim(irdG0wfile),exist=checkG0)
    if(.not.checkG0)inquire(file=trim(irdG0wfile)//".gz",exist=checkG0)

    !If irdG0wfile EXIST: read it and build the guess-->neqIPT
    if(checkG0)then
       call msg(bold("Continuing the EQUILIBRIUM SOLUTION to the KBM-Contour"))

       !READ \GG_0(w) --> eq_G0w(:)
       L=file_length(trim(irdG0wfile))
       allocate(eq_G0w(L))
       allocate(eq_Siw(L))
       allocate(eq_nk(Lk))
       allocate(wr_(L))
       call sread(trim(irdG0wfile),wr_,eq_G0w)
       fmesh_=abs(wr_(2)-wr_(1)) !Get G0 mesh:

       !Get G0(iw)
       allocate(dummyGiw(L),dummyGtau(0:Ltau))
       call get_matsubara_gf_from_DOS(wr_,eq_G0w,dummyGiw,beta)
       call fftgf_iw2tau(dummyGiw,dummyGtau,beta)

       !Get S(iw) w/ IPT:
       allocate(Stau(0:L))
       forall(i=0:L)Stau(i)=U**2*(dummyGtau(i))**2*dummyGtau(L-i)
       call fftgf_tau2iw(Stau,eq_Siw,beta)

       !Get interacting momentum-distribution n(k):
       eq_nk = square_lattice_momentum_distribution(Lk)
       !call read_nkfile(trim(irdnkfile))

       if(mpiID==0)then
          call splot("GUESS/eq_G0_iw.ipt",wm,dummyGiw)
          call splot("GUESS/eq_G_tau.ipt",tau,-dummyGtau(Ltau:0:-1))
          call splot("GUESS/eq_Sigma_iw.ipt",wm,eq_Siw)
          call splot("GUESS/eq_nkVSepsk.ipt",epsik,eq_nk)
       endif
       deallocate(dummyGiw,dummyGtau,Stau)

       G0less=zero
       G0gtr=zero
       G0lmix=zero
       do ik=1,L
          en   = wr_(ik)
          nless= fermi0(en,beta)
          ngtr = fermi0(en,beta)-1.d0
          A    = -aimag(eq_G0w(ik))/pi*fmesh_
          do i=0,nstep
             do j=0,nstep
                peso=exp(-xi*en*(t(i)-t(j)))
                G0less(i,j)=G0less(i,j) + xi*nless*A*peso
                G0gtr(i,j) =G0gtr(i,j)  + xi*ngtr*A*peso
             enddo
          enddo
          do i=0,nstep
             do j=0,Ltau
                peso=exp(-xi*en*t(i))*exp(-en*tau(j))/(1.d0+exp(beta*en))
                if(beta*en>35.d0)peso=exp(-xi*en*t(i))*exp(-(tau(j)+beta)*en)
                G0lmix(i,j)=G0lmix(i,j) + xi*A*peso
             enddo
          enddo
       enddo
       forall(j=0:Ltau)G0gmix(j,:)=-conjg(G0lmix(:,Ltau-j))

       G0mat=0.d0
       gmtau(0:Ltau)=xi*G0lmix(0,0:Ltau)
       forall(i=1:Ltau)gmtau(-i)=-gmtau(Ltau-i)
       forall(i=0:Ltau,j=0:Ltau)G0mat(i,j)=gmtau(j-i)

       deallocate(wr_)

       if(mpiID==0)then
          call splot("guessG0less.data",G0less(0:nstep,0:nstep))
          call splot("guessG0gtr.data",G0gtr(0:nstep,0:nstep))
          call splot("guessG0lmix.data",G0lmix(0:nstep,0:Ltau))
          call splot("guessG0mat.data",G0mat(0:Ltau,0:Ltau))
          call splot("guessG0less_t_t.3d",t(0:),t(0:),G0less(0:,0:))
          call splot("guessG0lmix_t_tau.3d",t(0:),tau(0:),G0lmix(0:,0:))
          call splot("guessG0mat_tau_tau.3d",tau(0:),tau(0:),G0mat(0:,0:))
          call system("mv -vf *guess*.ipt *guess*.3d GUESS/")
       endif
       call neq_solve_ipt()



    else !If irdG0wfile DOES NOT EXIST: start from the non-interacting HF solution Sigma=n-1/2=0
       call msg(bold("Starting from the non-interacting solution (U-quench!!)"))
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
       if(mpiID==0)call splot("GUESS/eq_nkVSepsk.ipt",epsik,eq_nk)

       !Get guess for Sigma from non-interacting solution (Hartree-Fock approx.):
       !Equilibrium functions:
       eq_Siw=zero
       !Sigma functions
       Sgtr=zero   ; Sless=zero !to be changed when-out-half-filling
       Sgmix=zero  ; Slmix=zero
       Smat=zero

       !These are built for comparison only.
       !Copy from the code above.
       ! !WF guess as non-interacting local GF
       ! G0less=zero ; G0gtr=zero
       ! G0lmix=zero ; G0gmix=zero
       ! G0mat=0.d0
       ! do ik=1,Lk
       !    en   = epsik(ik)
       !    nless= fermi0(en,beta)
       !    ngtr = fermi0(en,beta)-1.d0
       !    do i=0,nstep
       !       do j=0,nstep
       !          peso=exp(-xi*en*(t(i)-t(j)))
       !          G0less(i,j)=G0less(i,j) + xi*nless*wt(ik)*peso
       !          G0gtr(i,j) =G0gtr(i,j) + xi*ngtr*wt(ik)*peso
       !       enddo
       !       do j=0,Ltau
       !          peso=exp(-xi*en*t(i))*exp(-en*tau(j))/(1.d0+exp(beta*en))
       !          if(beta*en>35.d0)peso=exp(-xi*en*t(i))*exp(-(tau(j)+beta)*en)
       !          G0lmix(i,j)=G0lmix(i,j) + xi*nless*wt(ik)*peso
       !       enddo
       !    enddo
       ! enddo
       ! forall(j=0:Ltau)G0gmix(j,:)=-conjg(G0lmix(:,Ltau-j))

       ! gmtau(0:Ltau)=-dimag(G0lmix(0,0:Ltau))
       ! forall(i=1:Ltau)gmtau(-i)=-gmtau(Ltau-i)
       ! forall(i=0:Ltau,j=0:Ltau)G0mat(i,j)=-gmtau(j-i)

       ! call splot("G0less_t_t.ipt",t(0:),t(0:),G0less(0:,0:))
       ! call splot("G0lmix_t_tau.ipt",t(0:),tau(0:),G0lmix(0:,0:))
       ! call splot("G0mat_tau_tau.ipt",tau(0:),tau(0:),G0mat(0:,0:))
    endif


  contains

    function square_lattice_momentum_distribution(Lk) result(nk)
      integer            :: Lk
      integer            :: ik,i
      type(matsubara_gf) :: gm
      real(8)            :: nk(Lk)
      call allocate_gf(gm,L)
      do ik=1,Lk
         gm%iw=one/(xi*wm - epsik(ik) - eq_Siw)
         call fftgf_iw2tau(gm%iw,gm%tau,beta)
         nk(ik)=-gm%tau(L)         
      enddo
    end function square_lattice_momentum_distribution

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
    real(8),dimension(-Ltau:Ltau) :: dummyStau

    !Get SIgma:
    call msg("Get Sigma(t,t')")
    call system("if [ ! -d SIGMA ]; then mkdir SIGMA; fi")

    forall(i=0:nstep,j=0:nstep)
       Sgtr (i,j) = -(U**2)*(G0gtr(i,j)**2)*G0less(j,i)
       Sless(i,j) = -(U**2)*(G0less(i,j)**2)*G0gtr(j,i)
    end forall

    forall(i=0:nstep,itau=0:Ltau)Slmix(i,itau) = -(U**2)*(G0lmix(i,itau)**2)*G0gmix(itau,i)
    forall(j=0:Ltau)Sgmix(j,:)=-conjg(Slmix(:,Ltau-j))

    forall(i=0:Ltau,j=0:Ltau)Smat(i,j)=-(U**2)*(G0mat(i,j)**2)*G0mat(j,i)

    !Save data:
    if(mpiID==0)then
       call splot("Sless.data",Sless(0:nstep,0:nstep))
       call splot("Sgtr.data",Sgtr(0:nstep,0:nstep))
       call splot("Slmix.data",Slmix(0:nstep,0:Ltau))
       call splot("Smat.data",Smat(0:Ltau,0:Ltau))
       call splot("Sless_t_t.3d",t(0:),t(0:),Sless(0:,0:))
       call splot("Slmix_t_tau.3d",t(0:),tau(0:),Slmix(0:,0:))
       call splot("Smat_tau_tau.3d",tau(0:),tau(0:),Smat(0:,0:))
       call system("mv -vf *S*_t*.3d SIGMA/")
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
