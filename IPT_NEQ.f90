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
  !PURPOSE  : Initialize the run guessing/reading self-energy
  !+-------------------------------------------------------------------+
  subroutine neq_init_run()
    logical                             :: init

    init = inquire_kbm_contour_gf(trim(irdFILE))
    if(init)then
       call msg(bold("Reading components of the input Self-energy"))
       call read_kbm_contour_gf(Sigma,trim(irdFILE))
    else
       call msg(bold("Start from the Hartree-Fock self-energy"))
       Sigma=zero
    endif
  end subroutine neq_init_run




  !+-------------------------------------------------------------------+
  !PURPOSE  : BUild the 2^nd IPT sigma functions:
  !comments : !reintroducing wrong sign here, check vs. KADANOFFBAYM.f90/GFstep
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt()
    integer      :: i,j,itau
    real(8),dimension(0:nstep)            :: nt             !occupation(time)
    call msg("Get Sigma(t,t')")

    forall(i=0:Ltau)eq_Stau(i)=(U**2)*(eq_G0tau(i)**2)*eq_G0tau(Ltau-i)
    forall(i=1:Ltau)eq_Stau(-i)=-eq_Stau(Ltau-i)
    call fftgf_tau2iw(eq_Stau(0:),eq_Siw,beta)                   !Get S(iw) from S(tau)

    forall(i=0:Ltau,j=0:Ltau)Sigma%mats(i,j)=eq_Stau(i-j)        !get Sigma^M(tau,tau`)

    forall(i=0:nstep,j=0:nstep)
       Sigma%less(i,j) = (U**2)*(G0%less(i,j)**2)*G0%gtr(j,i)    !get Sigma^<(t,t`)
       Sigma%gtr (i,j) = (U**2)*(G0%gtr(i,j)**2)*G0%less(j,i)    !get Sigma^>(t,t`)
    end forall

    forall(i=0:nstep,itau=0:Ltau)&
         Sigma%lmix(i,itau)=(U**2)*(G0%lmix(i,itau)**2)*G0%gmix(itau,i) !get Sigma^\lmix(t,tau`)
    !    Sigma%gmix(itau,i)=(U**2)*(G0%gmix(itau,i)**2)*G0%lmix(i,itau) !get Sigma^\gmix(t,tau`)
    ! end forall
    forall(j=0:Ltau)Sigma%gmix(j,:)=conjg(Sigma%lmix(:,Ltau-j))    !get Sigma^\gmix(tau,t`)

    !Save data:
    if(mpiID==0)then
       call write_kbm_contour_gf(Sigma,trim(data_dir)//"/Sigma")
       if(plot3D)call plot_kbm_contour_gf(Sigma,t(0:),tau(0:),trim(plot_dir)//"/Sigma")
       call splot("eq_Sigma_tau.ipt",taureal,eq_Stau,append=.true.)
       call splot("eq_Sigma_iw.ipt",wm,eq_Siw,append=.true.)
       call splot("Sigma_less_t0.ipt",t(0:),Sigma%less(0:,0))
       call splot("Sigma_lmix_tau0.ipt",t(0:),Sigma%lmix(0:,0))
       forall(i=0:nstep)nt(i)=-xi*Sigma%less(i,i)
       call splot("nsVStime.ipt",t(0:nstep),nt(0:nstep),append=TT)
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
