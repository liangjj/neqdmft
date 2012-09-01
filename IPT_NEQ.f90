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

  integer :: Liw,Lw
  real(8),allocatable,dimension(:)    :: wr_,wm_

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Initialize the run guessing/reading/setting initial conditions
  !+-------------------------------------------------------------------+
  subroutine neq_init_run()
    integer                             :: i,j,ik
    real(8)                             :: en,A,fmesh_,imt
    real(8)                             :: nless,ngtr,xmu_,beta_
    complex(8)                          :: peso
    logical                             :: checkG0,checkG0w,checkG0iw

    real(8),dimension(-Ltau:Ltau)       :: tmpGtau
    real(8),allocatable,dimension(:)    :: eq_G0tau

    call create_data_dir("InitialConditions")

    inquire(file=trim(irdG0wfile),exist=checkG0w)
    if(.not.checkG0w)inquire(file=trim(irdG0wfile)//".gz",exist=checkG0w)
    inquire(file=trim(irdG0iwfile),exist=checkG0iw)
    if(.not.checkG0iw)inquire(file=trim(irdG0iwfile)//".gz",exist=checkG0iw)
    checkG0=checkG0w*checkG0iw

    !If irdG0wfile EXIST: read it and build the guess-->neqIPT
    if(checkG0)then
       call msg(bold("Continuing the EQUILIBRIUM SOLUTION to the KBM-Contour"))

       !READ \GG_0(w) --> eq_G0w(:)
       Lw=file_length(trim(irdG0wfile))
       allocate(eq_G0w(Lw),wr_(Lw))
       call sread(trim(irdG0wfile),wr_,eq_G0w)
       fmesh_=abs(wr_(2)-wr_(1)) !Get G0 mesh:

       !Get G0(iw)
       Liw=file_length(trim(irdG0iwfile))
       allocate(eq_G0iw(Liw),wm_(Liw))
       allocate(eq_G0tau(0:Ltau))
       call sread(trim(irdG0iwfile),wm_,eq_G0iw)
       ! call get_matsubara_gf_from_DOS(wr_,eq_G0w,dummyGiw,beta)
       call fftgf_iw2tau(eq_G0iw,eq_G0tau,beta)

       !Get S(iw) w/ IPT:
       allocate(eq_Stau(-Ltau:Ltau),eq_Siw(Liw))       
       forall(i=0:Ltau)eq_Stau(i)=U**2*(eq_G0tau(i))**2*eq_G0tau(Ltau-i)
       forall(i=1:Ltau)eq_Stau(-i)=-eq_Stau(Ltau-i)
       call fftgf_tau2iw(eq_Stau(0:),eq_Siw,beta)


       !Get interacting momentum-distribution n(k):
       allocate(eq_nk(Lk))
       eq_nk = square_lattice_momentum_distribution(Lk)
       !call read_nkfile(trim(irdnkfile))

       if(mpiID==0)then
          call splot("InitialConditions/ic_G0_iw.ipt",wm_,eq_G0iw)
          call splot("InitialConditions/ic_G0_tau.ipt",tau,-eq_G0tau(Ltau:0:-1))
          call splot("InitialConditions/ic_Sigma_iw.ipt",wm_,eq_Siw)
          call splot("InitialConditions/ic_Sigma_tau.ipt",tau,-eq_Stau(Ltau:0:-1))
          call splot("InitialConditions/ic_nkVSepsk.ipt",epsik,eq_nk)
       endif

       !Continue interacting solution to the KBM-Contour:
       G0=zero
       do ik=1,Lw
          en   = wr_(ik)
          nless= fermi0(en,beta)
          ngtr = fermi0(en,beta)-1.d0
          A    = -aimag(eq_G0w(ik))/pi*fmesh_
          do i=0,nstep
             do j=0,nstep
                peso=exp(-xi*en*(t(i)-t(j)))
                G0%less(i,j)=G0%less(i,j) + xi*nless*A*peso
                G0%gtr(i,j) =G0%gtr(i,j)  + xi*ngtr*A*peso
             enddo
          enddo
          do i=0,nstep
             do j=0,Ltau
                peso=exp(-xi*en*t(i))*exp(-en*tau(j))/(1.d0+exp(beta*en))
                if(beta*en>35.d0)peso=exp(-xi*en*t(i))*exp(-(tau(j)+beta)*en)
                G0%lmix(i,j)=G0%lmix(i,j) + xi*A*peso
             enddo
          enddo
       enddo
       forall(j=0:Ltau)G0%gmix(j,:)=conjg(G0%lmix(:,Ltau-j))

       tmpGtau(0:Ltau)= eq_G0tau(0:Ltau) !xi*G0lmix(0,0:Ltau)
       forall(i=1:Ltau)tmpGtau(-i)=-tmpGtau(Ltau-i)
       forall(i=0:Ltau,j=0:Ltau)G0%mats(i,j)=tmpGtau(j-i)

       deallocate(wr_,wm_)

       if(mpiID==0)then
          call write_kbm_contour_gf(G0,reg_filename(data_dir)//"/guessG0")
          if(plot3D)call plot_kbm_contour_gf(G0,t(0:),tau(0:),"PLOT/guessG0")
       endif
       call neq_solve_ipt()

    else !If irdG0wfile DOES NOT EXIST: start from the non-interacting HF solution Sigma=n-1/2=0

       call error("Can not find "//trim(irdG0wfile)//" and "//trim(irdG0iwfile)//" files")

       ! call msg(bold("Starting from the non-interacting solution (U-quench!!)"))
       ! !Get non-interacting n(k):
       ! !Initial conditions:
       ! allocate(eq_nk(Lk))
       ! allocate(eq_G0w(L))
       ! allocate(eq_Siw(L))
       ! xmu_=xmu   ; if(iquench)xmu_=xmu0
       ! beta_=beta ; if(iquench)beta_=beta0
       ! do ik=1,Lk
       !    eq_nk(ik)=fermi0((epsik(ik)-xmu_),beta_)
       ! enddo
       ! if(mpiID==0)call splot("InitialConditions/eq_nkVSepsk.ipt",epsik,eq_nk)

       ! !Get guess for Sigma from non-interacting solution (Hartree-Fock approx.):
       ! !Equilibrium functions:
       ! eq_Siw=zero
       ! !Sigma functions
       ! Sgtr=zero   ; Sless=zero !to be changed when-out-half-filling
       ! Sgmix=zero  ; Slmix=zero
       ! Smat=zero

       ! !These are built for comparison only.
       ! !Copy from the code above.
       ! ! !WF guess as non-interacting local GF
       ! ! G0less=zero ; G0gtr=zero
       ! ! G0lmix=zero ; G0gmix=zero
       ! ! G0mat=0.d0
       ! ! do ik=1,Lk
       ! !    en   = epsik(ik)
       ! !    nless= fermi0(en,beta)
       ! !    ngtr = fermi0(en,beta)-1.d0
       ! !    do i=0,nstep
       ! !       do j=0,nstep
       ! !          peso=exp(-xi*en*(t(i)-t(j)))
       ! !          G0less(i,j)=G0less(i,j) + xi*nless*wt(ik)*peso
       ! !          G0gtr(i,j) =G0gtr(i,j) + xi*ngtr*wt(ik)*peso
       ! !       enddo
       ! !       do j=0,Ltau
       ! !          peso=exp(-xi*en*t(i))*exp(-en*tau(j))/(1.d0+exp(beta*en))
       ! !          if(beta*en>35.d0)peso=exp(-xi*en*t(i))*exp(-(tau(j)+beta)*en)
       ! !          G0lmix(i,j)=G0lmix(i,j) + xi*nless*wt(ik)*peso
       ! !       enddo
       ! !    enddo
       ! ! enddo
       ! ! forall(j=0:Ltau)G0gmix(j,:)=-conjg(G0lmix(:,Ltau-j))

       ! ! gmtau(0:Ltau)=-dimag(G0lmix(0,0:Ltau))
       ! ! forall(i=1:Ltau)gmtau(-i)=-gmtau(Ltau-i)
       ! ! forall(i=0:Ltau,j=0:Ltau)G0mat(i,j)=-gmtau(j-i)

       ! ! call splot("G0less_t_t.ipt",t(0:),t(0:),G0less(0:,0:))
       ! ! call splot("G0lmix_t_tau.ipt",t(0:),tau(0:),G0lmix(0:,0:))
       ! ! call splot("G0mat_tau_tau.ipt",tau(0:),tau(0:),G0mat(0:,0:))

    endif


  contains

    function square_lattice_momentum_distribution(Lk) result(nk)
      integer            :: Lk
      integer            :: ik,i
      type(matsubara_gf) :: gm
      real(8)            :: nk(Lk)
      call allocate_gf(gm,Liw)
      do ik=1,Lk
         gm%iw=one/(xi*wm_ - epsik(ik) - eq_Siw)
         call fftgf_iw2tau(gm%iw,gm%tau,beta)
         nk(ik)=-gm%tau(Liw)         
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

    forall(i=0:nstep,j=0:nstep)
       S%gtr (i,j) = -(U**2)*(G0%gtr(i,j)**2)*G0%less(j,i)
       S%less(i,j) = -(U**2)*(G0%less(i,j)**2)*G0%gtr(j,i)
    end forall

    forall(i=0:nstep,itau=0:Ltau)S%lmix(i,itau) = -(U**2)*(G0%lmix(i,itau)**2)*G0%gmix(itau,i)
    forall(j=0:Ltau)S%gmix(j,:)=conjg(S%lmix(:,Ltau-j))

    forall(i=0:Ltau,j=0:Ltau)S%mats(i,j)=eq_Stau(j-i) 

    !Save data:
    if(mpiID==0)then
       call write_kbm_contour_gf(S,reg_filename(data_dir)//"/Sigma")
       if(plot3D)call plot_kbm_contour_gf(S,t(0:),tau(0:),"PLOT/Sigma")
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
  ! dSret=zero ; dG0ret=zero ; dGret=zero
  ! GammaRet=zero ; Gamma0Ret=zero
  ! !1 - get the Ret components of G_0 && \Sigma:
  ! forall(i=0:nstep,j=0:nstep)
  !    dG0ret(i,j)=heaviside(t(i)-t(j))*(G0gtr(i,j) - G0less(i,j))
  !    dSret(i,j) =heaviside(t(i)-t(j))*(Sgtr(i,j) - Sless(i,j))
  ! end forall
  ! !2 - get the  operator: \Gamma_0^R = \Id - \Sigma^R\circ G_0^R && invert it
  ! Uno=zero  ; forall(i=0:nstep)Uno(i,i)=One/dt
  ! Gamma0Ret(0:nstep,0:nstep) = Uno-matmul(dSret(0:nstep,0:nstep),dG0ret(0:nstep,0:nstep))*dt
  ! Gamma0Ret(0:nstep,0:nstep)=Gamma0Ret(0:nstep,0:nstep)*dt**2
  ! call mat_inversion_GJ(Gamma0Ret(0:nstep,0:nstep))
  ! !3 - get G_imp^R, G_imp^{>,<} using Dyson equations:
  ! dGret(0:nstep,0:nstep)    = matmul(dG0ret(0:nstep,0:nstep),Gamma0Ret(0:nstep,0:nstep))*dt 
  ! GammaRet(0:nstep,0:nstep) = Uno + matmul(dGret(0:nstep,0:nstep),dSret(0:nstep,0:nstep))*dt

  ! impGless(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),&
  !      matmul(G0less(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep)))))*dt**2 +&
  !      matmul(dGret(0:nstep,0:nstep),matmul(Sless(0:nstep,0:nstep),conjg(transpose(dGret(0:nstep,0:nstep)))))*dt**2

  ! impGgtr(0:nstep,0:nstep)  = matmul(GammaRet(0:nstep,0:nstep),&
  !      matmul(G0gtr(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep)))))*dt**2  +&
  !      matmul(dGret(0:nstep,0:nstep),matmul(Sgtr(0:nstep,0:nstep),conjg(transpose(dGret(0:nstep,0:nstep)))))*dt**2
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
