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

  integer                             :: Lw


contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Initialize the run guessing/reading/setting initial conditions
  !+-------------------------------------------------------------------+
  subroutine neq_init_run()
    integer                             :: i,j,ik
    real(8)                             :: en,A,fmesh_,imt
    real(8)                             :: nless,ngtr,xmu_,beta_
    complex(8)                          :: peso
    logical                             :: check(3),init
    real(8),dimension(-Ltau:Ltau)       :: tmpGtau
    real(8),allocatable,dimension(:)    :: wr_,wm_

    call create_data_dir("InitialConditions")
    init=.true.
    do i=1,3
       inquire(file=trim(irdFILE(i)),exist=check(i))
       if(.not.check(i))inquire(file=trim(irdFILE(i))//".gz",exist=check(i))
       init=init.AND.check(i)
    enddo

    !If irdG0wfile EXIST: read it and build the guess-->neqIPT
    if(init)then
       call msg(bold("Reading the EQUILIBRIUM SOLUTION"))

       !READ \GG_0(w) --> eq_G0w(:)
       Lw=file_length(trim(irdFILE(1)))
       allocate(eq_G0w(Lw),wr_(Lw))
       call sread(trim(irdFILE(1)),wr_,eq_G0w)
       fmesh_=wr_(2)-wr_(1)

       !Get S(tau):
       Ltau=file_length(trim(irdFILE(2)))-1
       allocate(eq_Stau(-Ltau:Ltau),eq_Siw(L),wm_(L))
       call sread(trim(irdFILE(2)),eq_Stau(0:))
       forall(i=1:Ltau)eq_Stau(-i)=-eq_Stau(Ltau-i)
       wm_ = pi/beta*real(2*arange(1,L)-1,8)
       call fftgf_tau2iw(eq_Stau(0:),eq_Siw,beta)

       !Get n(k):
       allocate(eq_nk(Lk))
       eq_nk = get_square_lattice_momentum_distribution(eq_Siw,Lk)

       !READ \GG_0(tau)
       if(Ltau/=file_length(trim(irdFILE(3)))-1)call error("A problem occur reading "//trim(irdFILE(3)))
       allocate(eq_G0tau(-Ltau:Ltau))
       call sread(trim(irdFILE(3)),eq_G0tau(0:))
       forall(i=1:Ltau)eq_G0tau(-i)=-eq_G0tau(Ltau-i)

       if(mpiID==0)then
          call splot("InitialConditions/ic_G0_realw.ipt",wr_,eq_G0w)
          call splot("InitialConditions/ic_Sigma_iw.ipt",wm_,eq_Siw)
          call splot("InitialConditions/ic_Sigma_tau.ipt",taureal,eq_Stau)
          call splot("InitialConditions/ic_G0_tau.ipt",taureal,eq_G0tau)
          call splot("InitialConditions/ic_nkVSepsk.ipt",epsik,eq_nk)
       endif

       !Continue interacting solution to the KBM-Contour:
       call msg(bold("Continuing the EQUILIBRIUM SOLUTION to the KBM-Contour"))
       G0=zero
       do ik=1,Lw
          en   = wr_(ik)
          nless= fermi0(en,beta)
          ngtr = fermi0(en,beta)-1.d0
          A    = -dimag(eq_G0w(ik))/pi*fmesh_
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

       ! tmpGtau(0:Ltau)= eq_G0tau(0:Ltau) !xi*G0lmix(0,0:Ltau)
       ! forall(i=1:Ltau)tmpGtau(-i)=-tmpGtau(Ltau-i)
       forall(i=0:Ltau,j=0:Ltau)G0%mats(i,j)=eq_G0tau(j-i)

       if(mpiID==0)then
          call write_kbm_contour_gf(G0,reg_filename(data_dir)//"/guessG0")
          if(plot3D)call plot_kbm_contour_gf(G0,t(0:),tau(0:),"PLOT/guessG0")
       endif

       call neq_solve_ipt()

       deallocate(wr_,wm_)

    else !If irdG0wfile DOES NOT EXIST: start from the non-interacting HF solution Sigma=n-1/2=0

       call error("Can not find "//trim(irdFILE(1))//", "//trim(irdFILE(2))//", "//trim(irdFILE(3))//" files")

    endif



  contains

    function get_square_lattice_momentum_distribution(Siw,Lk) result(nk)
      integer            :: Lk
      integer            :: ik,i
      complex(8)         :: Siw(:)
      type(matsubara_gf) :: gm
      real(8)            :: nk(Lk)
      call allocate_gf(gm,L)
      do ik=1,Lk
         gm%iw=one/(xi*wm_ - epsik(ik) - Siw)
         call fftgf_iw2tau(gm%iw,gm%tau,beta)
         nk(ik)=-gm%tau(L)         
      enddo
    end function get_square_lattice_momentum_distribution

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

    !Get SIgma:
    call msg("Get Sigma(t,t')")

    Sigma=zero
    forall(i=0:nstep,j=0:nstep)
       Sigma%gtr (i,j) = (U**2)*(G0%gtr(i,j)**2)*G0%less(j,i) !reintroducing wrong sign here, check vs. KADANOFFBAYM.f90/GFstep
       Sigma%less(i,j) = (U**2)*(G0%less(i,j)**2)*G0%gtr(j,i)
    end forall

    forall(i=0:nstep,itau=0:Ltau)Sigma%lmix(i,itau) = (U**2)*(G0%lmix(i,itau)**2)*G0%gmix(itau,i)
    forall(j=0:Ltau)Sigma%gmix(j,:)=conjg(Sigma%lmix(:,Ltau-j))

    forall(i=0:Ltau,j=0:Ltau)Sigma%mats(i,j)=eq_Stau(j-i)

    !Save data:
    if(mpiID==0)then
       call write_kbm_contour_gf(Sigma,reg_filename(data_dir)//"/Sigma")
       if(plot3D)call plot_kbm_contour_gf(Sigma,t(0:),tau(0:),"PLOT/Sigma")
    endif

  end subroutine neq_solve_ipt



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
  !      dSret=zero ; dG0ret=zero ; dGret=zero
  !      GammaRet=zero ; Gamma0Ret=zero
  !      !1 - get the Ret components of G_0 && \Sigma:
  !      forall(i=0:nstep,j=0:nstep)
  !         dG0ret(i,j)=heaviside(t(i)-t(j))*(G0gtr(i,j) - G0less(i,j))
  !         dSret(i,j) =heaviside(t(i)-t(j))*(Sgtr(i,j) - Sless(i,j))
  !      end forall
  !      !2 - get the  operator: \Gamma_0^R = \Id - \Sigma^R\circ G_0^R && invert it
  !      Uno=zero  ; forall(i=0:nstep)Uno(i,i)=One/dt
  !      Gamma0Ret(0:nstep,0:nstep) = Uno-matmul(dSret(0:nstep,0:nstep),dG0ret(0:nstep,0:nstep))*dt
  !      Gamma0Ret(0:nstep,0:nstep)=Gamma0Ret(0:nstep,0:nstep)*dt**2
  !      call mat_inversion_GJ(Gamma0Ret(0:nstep,0:nstep))
  !      !3 - get G_imp^R, G_imp^{>,<} using Dyson equations:
  !      dGret(0:nstep,0:nstep)    = matmul(dG0ret(0:nstep,0:nstep),Gamma0Ret(0:nstep,0:nstep))*dt 
  !      GammaRet(0:nstep,0:nstep) = Uno + matmul(dGret(0:nstep,0:nstep),dSret(0:nstep,0:nstep))*dt

  !      impGless(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),&
  !           matmul(G0less(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep)))))*dt**2 +&
  !           matmul(dGret(0:nstep,0:nstep),matmul(Sless(0:nstep,0:nstep),conjg(transpose(dGret(0:nstep,0:nstep)))))*dt**2

  !      impGgtr(0:nstep,0:nstep)  = matmul(GammaRet(0:nstep,0:nstep),&
  !           matmul(G0gtr(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep)))))*dt**2  +&
  !           matmul(dGret(0:nstep,0:nstep),matmul(Sgtr(0:nstep,0:nstep),conjg(transpose(dGret(0:nstep,0:nstep)))))*dt**2
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
