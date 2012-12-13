!########################################################################
!PROGRAM  : KADANOFFBAYM
!TYPE     : module
!PURPOSE  : Evolve the Green's functions G^>,< on the real time axis
!given the initial condition at t=t'=0.
!AUTHORS  : Adriano Amaricci
!########################################################################
module KADANOFBAYM
  !LOCAL:
  USE VARS_GLOBAL
  USE ELECTRIC_FIELD
  USE INTEGRATE
  implicit none
  private
  !k-dependent GF:
  type(kbm_contour_gf)                  :: Gk
  complex(8),dimension(:),allocatable   :: Gkiw
  !Vector for KB propagation solution
  complex(8),allocatable,dimension(:)   :: Ikless,Ikgtr
  complex(8),allocatable,dimension(:)   :: Ikless0,Ikgtr0
  complex(8),allocatable,dimension(:)   :: Iklmix,Iklmix0
  real(8)                               :: Ikdiag
  !Auxiliary operators:
  complex(8),allocatable,dimension(:,:) :: Udelta,Vdelta
  public                                :: neq_get_localgf

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the k-sum to get the local Green's functions (R,>,<)
  !+-------------------------------------------------------------------+
  subroutine neq_get_localgf()
    if(solve_wfftw)then
       call get_equilibrium_localgf
    else
       call kadanoff_baym_to_localgf()
    endif
    if(fchi)call get_chi
    call print_out_Gloc()
  end subroutine neq_get_localgf


  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+-------------------------------------------------------------------+
  !PURPOSE  : Build the local Green's functions G^<,>
  !along the real time axis, using a discretized verion of the Kadanoff-Baym 
  !equations. The evolution is performed by a two-pass procedure, following:
  ! H.S. K\"ohler, N.H.Kwong, Hashim A. Yousif
  !"A Fortran code for solving the Kadanoff-Baym equations for a homogeneous 
  !fermion system" Computer Physics Communications,123(1999),123-142
  !+-------------------------------------------------------------------+
  subroutine kadanoff_baym_to_localgf()
    integer              :: istep,i,j,ik
    complex(8)           :: tmpGkiw(L)
    real(8)              :: tmpnk(0:Nstep,Lk)
    type(kbm_contour_gf) :: tmpG


    !Allocate all arrays
    call allocate_funx
    call allocate_kbm_contour_gf(tmpG,Nstep,Ltau)
    tmpGkiw=zero
    tmpnk  =0.d0

    call buildUV

    call msg("Entering Kadanoff-Baym")
    call msg("Using integration method: "//bold(trim(int_method)))
    call msg("Using Ltau slices: "//bold(txtfy(Ltau)))
    !Set to Zero GF and nk:
    eq_Giw = zero
    locG   = zero
    nk     = 0.d0

    !call setup_kbe_kernel()

    !=============START K-POINTS LOOP======================
    call start_timer
    do ik=1+mpiID,Lk,mpiSIZE
       Gk  =zero
       Gkiw=zero

       !t-step loop
       call setup_initial_conditions(ik)
       do istep=0,nstep-1
          call step_kbm_contour_gf(ik,istep)
       enddo

       !sum over k-point
       tmpGkiw          = tmpGkiw          + Gkiw*wt(ik)
       tmpG%less(0:,0:) = tmpG%less(0:,0:) + Gk%less(0:,0:)*wt(ik)
       tmpG%gtr(0:,0:)  = tmpG%gtr(0:,0:)  + Gk%gtr(0:,0:)*wt(ik)
       tmpG%lmix(0:,0:) = tmpG%lmix(0:,0:) + Gk%lmix(0:,0:)*wt(ik)
       forall(istep=0:nstep)tmpnk(istep,ik)=-xi*Gk%less(istep,istep)
       call eta(ik,Lk)
    enddo
    call stop_timer
    call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
    !=============END K-POINTS LOOP======================

    !Get equilibrium-Matsubara GF:
    call MPI_ALLREDUCE(tmpGkiw,eq_Giw,L,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    if(upmflag)then
       call fftgf_iw2tau_upm(wm,eq_Giw,tau(0:),locG%mats(0:),beta)
    else
       call fftgf_iw2tau(eq_Giw,locG%mats(0:),beta)
    endif
    forall(i=1:Ltau)locG%mats(-i)=-locG%mats(Ltau-i)

    !Reduce other components of the Contour GF:
    call MPI_ALLREDUCE(tmpG%less(0:,0:),locG%less(0:,0:),(Nstep+1)**2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(tmpG%gtr(0:,0:),locG%gtr(0:,0:),(Nstep+1)**2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(tmpG%lmix(0:,0:),locG%lmix(0:,0:),(Nstep+1)*(Ltau+1),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    forall(i=0:Ltau)locG%gmix(i,0:)=conjg(locG%lmix(0:,Ltau-i))
    call MPI_ALLREDUCE(tmpnk,nk,(nstep+1)*Lk,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)

    !deallocating memory:
    call deallocate_kbm_contour_gf(tmpG)
    call deallocate_funx
  end subroutine kadanoff_baym_to_localgf



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+----------------------------------------------------------------+
  !PURPOSE  :  allocation of working array
  !+----------------------------------------------------------------+
  subroutine allocate_funx()
    call msg("Allocating KB memory:")
    !Predictor-corrector solver arrays: store the time-step
    allocate(Ikless(0:nstep),Ikgtr(0:nstep))
    allocate(Ikless0(0:nstep),Ikgtr0(0:nstep))
    allocate(Iklmix(0:Ltau),Iklmix0(0:Ltau))
    !Allocate k-dependent GF:
    call allocate_kbm_contour_gf(Gk,Nstep,Ltau)
    allocate(Gkiw(L))
    !Aux. operators
    allocate(Udelta(Lk,0:nstep),Vdelta(Lk,0:nstep))
    !allocate(Ktau(0:Ltau))
    !Chi
    if(fchi)allocate(chi_dia(2,2,0:nstep),chi_pm(2,2,0:nstep,0:nstep))
  end subroutine allocate_funx

  subroutine deallocate_funx()
    call msg("Deallocating KB memory:")
    !Predictor-corrector solver arrays: store the time-step
    deallocate(Ikless,Ikgtr)
    deallocate(Ikless0,Ikgtr0)
    deallocate(Iklmix,Iklmix0)
    !Allocate k-dependent GF:
    call deallocate_kbm_contour_gf(Gk)
    deallocate(Gkiw)
    !Aux. operators
    deallocate(Udelta,Vdelta)
    !Chi
    if(fchi)deallocate(chi_dia,chi_pm)
  end subroutine deallocate_funx



  !******************************************************************
  !******************************************************************
  !******************************************************************



  ! !+-------------------------------------------------------------------+
  ! !PURPOSE: Setup the kernel of the of the integral part for the 
  ! !KBE. In this case we have the Sigma function. The *tau integrals
  ! !are particularly bad behaving, so we use Uniform Power Mesh here.
  ! !To this end we fit or evaluate the kernel on the non-uniform mesh.
  ! !+-------------------------------------------------------------------+
  ! subroutine setup_kbe_kernel()
  !   integer :: i,ik
  !   complex(8),dimension(0:Nstep,0:Ltau) :: uSlmix,nSlmix
  !   character(len=12) :: file
  !   file='KSigma'
  !   KSigma%less(0:,0:)=Sigma%less(0:,0:)
  !   KSigma%gtr(0:,0:)=Sigma%gtr(0:,0:)
  !   KSigma%lmix(0:,0:)=Sigma%lmix(0:,0:)
  !   KSigma%gmix(0:,0:)=Sigma%gmix(0:,0:)
  !   KSigma%mats(:)=Sigma%mats(:)
  !   if(upmflag)then
  !      nSlmix(0:,0:)=Sigma%lmix(0:,0:)
  !      do i=0,Nstep
  !         call cubic_spline(nSlmix(i,0:),tau(0:),uSlmix(i,0:),ftau(0:))
  !      enddo
  !      KSigma%lmix(0:,0:)=uSlmix(0:,0:)
  !      forall(i=0:Ltau)KSigma%gmix(i,0:)=conjg(KSigma%lmix(0:,Ltau-i))
  !      call splot(trim(file)//"_less_t_t",t(0:),t(0:),KSigma%less(0:,0:))
  !      call splot(trim(file)//"_gtr_t_t",t(0:),t(0:),KSigma%gtr(0:,0:))
  !      call splot(trim(file)//"_lmix_t_tau",t(0:),ftau(0:),KSigma%lmix(0:,0:))
  !      call splot(trim(file)//"_gmix_tau_t",ftau(0:),t(0:),KSigma%gmix(0:,0:))
  !      call splot(trim(file)//"_mats_tau",ftau(0:),KSigma%lmix(0,0:))
  !      Ktau(0:)=ftau(0:)
  !   else
  !      Ktau(0:)=tau(0:)
  !   endif
  ! end subroutine setup_kbe_kernel



  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine setup_initial_conditions(ik)
    integer :: i,j,ik
    Gkiw = one/(xi*wm -epsik(ik) -eq_Siw)                  !get G_k(iw)
    if(upmflag)then
       call fftgf_iw2tau_upm(wm,Gkiw,tau(0:),Gk%mats(0:),beta) !get G_k(tau>0)
    else
       call fftgf_iw2tau(Gkiw,Gk%mats(0:),beta)
    endif
    forall(i=1:Ltau)Gk%mats(-i)=-Gk%mats(Ltau-i)               !get G_k(tau<0)          
    Gk%less(0,0) = -xi*Gk%mats(Ltau)                         !get G^<_k(0,0)=xi*G_k(0-)=-G_k(beta)
    Gk%gtr(0,0)  =  xi*Gk%mats(0)                            !get G^>_k(0,0)=xi*G_k(0+)=xi*(G_k(0-)-1.0)
    forall(i=0:Ltau)                                       !
       Gk%lmix(0,i)=-xi*Gk%mats(Ltau-i)                      !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
       Gk%gmix(i,0)= xi*Gk%mats(i)                           !get G^\gmix_k(tau,0)=xi*G_k(tau>0)
    end forall
  end subroutine setup_initial_conditions


  !+-------------------------------------------------------------------+
  !PURPOSE  : Execute the 2pass procedure to solve KB equations:
  !+-------------------------------------------------------------------+
  subroutine step_kbm_contour_gf(ik,istep)
    integer    :: ips,istep,ik
    integer    :: i,j,itau
    integer    :: it
    do ips=1,2
       select case(ips)
       case(1)
          !First step: get collision integrals up to t=T=istep
          call get_Ikcollision(istep)
          Ikless0 = Ikless ; Ikgtr0  = Ikgtr ; Iklmix0 = Iklmix 
          Ikdiag  =-2.d0*dreal(Ikless(istep))
       case(2)
          !Second step: get collision integrals up to t=T+\Delta=istep+1
          call get_Ikcollision(istep+1)
          Ikless   = (Ikless  + Ikless0)/2.d0
          Ikgtr    = (Ikgtr   + Ikgtr0)/2.d0
          Iklmix   = (Iklmix  + Iklmix0)/2.d0
          Ikdiag   = (-2.d0*dreal(Ikless(istep+1)) + Ikdiag)/2.d0 
       end select
       !
       !Evolve the solution of KB equations for all the k-points:
       forall(it=0:istep)
          Gk%less(it,istep+1)  = Gk%less(it,istep)*conjg(Udelta(ik,istep))-&
               Ikless(it)*conjg(Vdelta(ik,istep))
          Gk%gtr(istep+1,it)   = Gk%gtr(istep,it)*Udelta(ik,istep)-&
               Ikgtr(it)*Vdelta(ik,istep)
       end forall
       forall(itau=0:Ltau)Gk%lmix(istep+1,itau)=Gk%lmix(istep,itau)*Udelta(ik,istep)-&
            Iklmix(itau)*Vdelta(ik,istep)
       Gk%gtr(istep+1,istep)   =(Gk%less(istep,istep)-xi)*Udelta(ik,istep)-Ikgtr(istep)*Vdelta(ik,istep)
       Gk%less(istep+1,istep+1)= Gk%less(istep,istep)-xi*dt*Ikdiag
       Gk%gtr(istep+1,istep+1) = Gk%less(istep+1,istep+1)-xi
       !
       do i=0,istep+1
          Gk%less(istep+1,i)=-conjg(Gk%less(i,istep+1))
          Gk%gtr(i,istep+1) =-conjg(Gk%gtr(istep+1,i))
       enddo
       !
       forall(itau=0:Ltau)Gk%gmix(itau,istep+1)=conjg(Gk%lmix(istep+1,Ltau-itau))
    enddo
  end subroutine Step_kbm_contour_gf



  !******************************************************************
  !******************************************************************
  !******************************************************************






  !+-------------------------------------------------------------------+
  !COMMENTS : VER3.0_Apr2010. This version is based on the homogeneity
  ! with repects to k' of the matrices F_{k,k'}(t,t'), related to the 
  ! constant expression of Vhyb{k,k'}\== Vpd.
  ! A more general form has been tested in the past, using direct
  ! expression of the matrices in {k,k'}. See repo.
  !+-------------------------------------------------------------------+
  subroutine get_Ikcollision(Nt)
    integer,intent(in)          :: Nt
    integer                     :: i,j,itau,it,itp
    complex(8)                  :: I1,I2,Ib
    complex(8),dimension(0:Nt)  :: Vadv,Vret,Vless,Vgtr
    complex(8),dimension(0:Ltau):: Vlmix,Vgmix
    complex(8),dimension(0:Nt)  :: Fadv,Fret
    Ikless=zero; Ikgtr=zero; Iklmix=zero
    if(Nt==0)return

    do i=0,Nt
       Vless(i)= Sigma%less(i,Nt) + S0%less(i,Nt)
       Vgtr(i) = Sigma%gtr(Nt,i)  + S0%gtr(Nt,i)
       Vret(i) = SretF(Nt,i)      + S0retF(Nt,i)
       Vadv(i) = conjg(Vret(i))
    end do
    do i=0,Ltau
       Vlmix(i)= Sigma%lmix(Nt,i) + S0%lmix(Nt,i)
       Vgmix(i)= Sigma%gmix(i,Nt) + S0%gmix(i,Nt)
    end do

    select case(int_method)
    case default             !trapz
       do it=0,Nt                                  !for all t=0:T//T+\Delta
          forall(i=0:it)Fret(i) = GkretF(it,i)     !was: forall(i=0:NT)          
          I1=trapz(dt,Fret(0:it)*Vless(0:it))
          I2=trapz(dt,Gk%less(it,0:Nt)*Vadv(0:Nt))
          Ib=trapz(tau(0:Ltau),Gk%lmix(it,0:Ltau)*Vgmix(0:Ltau)) !
          Ikless(it)=I1 + I2 - Ib*xi
       enddo
       !
       do itp=0,Nt
          forall(i=0:itp)Fadv(i) = conjg(GkretF(itp,i)) !was: forall(i=0:Nt)
          I1=trapz(dt,Vret(0:Nt)*Gk%gtr(0:Nt,itp))
          I2=trapz(dt,Vgtr(0:itp)*Fadv(0:itp))
          Ib=trapz(tau(0:Ltau),Vlmix(0:Ltau)*Gk%gmix(0:Ltau,itp))
          Ikgtr(itp)=I1 + I2 - Ib*xi
       enddo
       !
       do itau=0,Ltau
          I1=trapz(dt,Vret(0:Nt)*Gk%lmix(0:Nt,itau))
          Ib=trapz(tau(0:Ltau),Vlmix(0:Ltau)*Gk%mats((0-itau):(Ltau-itau)))!*Gk%mats(0:Ltau,itau))
          Iklmix(itau)=I1 + Ib
       enddo

    case("simps")
       do it=0,Nt                                  !for all t=0:T//T+\Delta
          forall(i=0:it)Fret(i) = GkretF(it,i)     !was: forall(i=0:NT)
          I1=simps(dt,Fret(0:it)*Vless(0:it))           !
          I2=simps(dt,Gk%less(it,0:Nt)*Vadv(0:Nt))      !
          Ib=simps(tau(0:Ltau),Gk%lmix(it,0:Ltau)*Vgmix(0:Ltau)) !
          Ikless(it)=I1 + I2 - Ib*xi
       enddo
       !
       do itp=0,Nt
          forall(i=0:itp)Fadv(i) = conjg(GkretF(itp,i)) !was: forall(i=0:Nt)
          I1=simps(dt,Vret(0:Nt)*Gk%gtr(0:Nt,itp))           !
          I2=simps(dt,Vgtr(0:itp)*Fadv(0:itp))               !
          Ib=simps(tau(0:Ltau),Vlmix(0:Ltau)*Gk%gmix(0:Ltau,itp))     !
          Ikgtr(itp)=I1 + I2 - Ib*xi
       enddo
       !
       do itau=0,Ltau
          I1=simps(dt,Vret(0:Nt)*Gk%lmix(0:Nt,itau))
          Ib=simps(tau(0:Ltau),Vlmix(0:Ltau)*Gk%mats((0-itau):(Ltau-itau)))!*Gk%mats(0:Ltau,itau))
          Iklmix(itau)=I1 + Ib
       enddo

    case ("rect")               !The older version all comments kept
       !Get I^<(t=it,t`=Nt/T) it:0,...,Nt == t=0,...,T+\Delta
       !I1 =  \int_0^t  G^R*Sigma^<     = sum_i=0^t G^R(t,i)Sigma^<(i,T)
       !I2 =  \int_0^T  G^<*Sigma^A     = sum_i=0^T G^<(t,i)Sigma^A(i,T)
       !Ib =-i\int_0^\b G^\lmix*S^\gmix = -isum_i=0^\b G^\lmix(t,i)Sigma^\gmix(i,T)
       !$OMP PARALLEL SHARED(Ikless) PRIVATE(i,it,I1,I2,Fret,Fadv)
       !$OMP DO
       do it=0,Nt                                  !for all t=0:T//T+\Delta
          forall(i=0:it)Fret(i) = GkretF(it,i)     !was: forall(i=0:NT)
          I1=sum(Fret(0:it)*Vless(0:it))           !
          I2=sum(Gk%less(it,0:Nt)*Vadv(0:Nt))      !
          Ib=sum(Gk%lmix(it,0:Ltau)*Vgmix(0:Ltau)) !
          Ikless(it)=I1*dt + I2*dt - Ib*xi*dtau
       enddo
       !$OMP END DO
       !$OMP END PARALLEL
       !
       !Get I^>(t=Nt,t`=itp) itp:0,...,Nt == t` = 0,...,T+\Delta=Nt
       !I1 =  \int_0^T  S^R*G^>         = sum_i=0^T  Sigma^R(T,i)G^>(i,t`=itp)
       !I2 =  \int_0^t` S^>*G^A         = sum_i=0^t` Sigma^>(T,i)G^A(i,t`) G^A=0 for i>t`
       !Ib =-i\int_0^\b S^\lmix*G^\gmix = -isum_i=0^\b Sigma^\lmix(T,i)G^\gmix(i,t`)
       do itp=0,Nt
          forall(i=0:itp)Fadv(i) = conjg(GkretF(itp,i)) !was: forall(i=0:Nt)
          I1=sum(Vret(0:Nt)*Gk%gtr(0:Nt,itp))           !
          I2=sum(Vgtr(0:itp)*Fadv(0:itp))               !
          Ib=sum(Vlmix(0:Ltau)*Gk%gmix(0:Ltau,itp))     !
          Ikgtr(itp)=I1*dt + I2*dt - Ib*xi*dtau
       enddo
       !
       !Get I^\lmix(Nt,itau) itau=0:Ltau == -tau=0:beta
       !I1 = \int_0^T    S^Ret*G^\rceil = sum_i=0^L Sigma^R(T,i)G^\lmix(i,tau)
       !Ib = \int_0^\b S^\lceil* G^\M   = sum_i=0^L Sigma^\lmix(T,i)(xi*G(i,tau))
       do itau=0,Ltau
          I1=sum(Vret(0:Nt)*Gk%lmix(0:Nt,itau))
          Ib=sum(Vlmix(0:Ltau)*Gk%mats((0-itau):(Ltau-itau)))!*Gk%mats(0:Ltau,itau))
          Iklmix(itau)=I1*dt + Ib*dtau
       enddo
    end select
  end subroutine get_Ikcollision


  !******************************************************************
  !******************************************************************
  !******************************************************************


  pure function GkretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: GkretF
    GkretF = heaviside(t(i)-t(j))*(Gk%gtr(i,j)-Gk%less(i,j))
  end function GkretF

  pure function S0retF(i,j)
    integer,intent(in) :: i,j
    complex(8)         :: S0retF
    S0retF = heaviside(t(i)-t(j))*(S0%gtr(i,j)-S0%less(i,j))
  end function S0retF

  pure function SretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: SretF
    SretF = heaviside(t(i)-t(j))*(Sigma%gtr(i,j)-Sigma%less(i,j))
  end function SretF

  ! pure function KSretF(i,j)      
  !   integer,intent(in) :: i,j
  !   complex(8)         :: KSretF
  !   KSretF = heaviside(t(i)-t(j))*(KSigma%gtr(i,j)-KSigma%less(i,j))
  ! end function KSretF


  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine buildUV
    integer :: ik,i
    do  ik=1,Lk
       do i=0,nstep
          Udelta(ik,i)=UdeltaF(ik,i)
          Vdelta(ik,i)=VdeltaF(ik,i)
       enddo
    enddo
  end subroutine buildUV

  function UdeltaF(ik,istep) 
    integer,intent(in) :: ik,istep
    complex(8)         :: UdeltaF
    real(8)            :: arg
    arg=Hbar(ik,istep)
    UdeltaF=exp(-xi*arg*dt)
  end function UdeltaF

  function VdeltaF(ik,istep)
    integer,intent(in) :: ik,istep
    complex(8)         :: VdeltaF
    real(8)            :: arg
    arg=Hbar(ik,istep)
    VdeltaF=exp(-xi*arg*dt)
    if(abs(arg*dt) <= 1.d-5)then
       VdeltaF=xi*dt
    else
       VdeltaF=(1.d0-VdeltaF)/arg
    endif
  end function VdeltaF



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : This is the time-dependent hamiltonian
  ! a more general code would require this quantity to be definied 
  ! externally as an array for every k-point (ik) and time (istep) but
  ! I am lazy...
  !+-------------------------------------------------------------------+
  function Hbar(ik,istep)
    integer,intent(in) :: ik,istep  
    integer      :: i,j
    real(8)      :: Hbar
    real(8)      :: tbar
    type(vect2D) :: kt,Ak
    tbar=t(istep) + dt/2.d0
    i=ik2ix(ik)
    j=ik2iy(ik)
    Ak=Afield(tbar,Ek)
    kt=kgrid(i,j) - Ak
    Hbar=square_lattice_dispersion(kt)
  end function Hbar




  !******************************************************************
  !******************************************************************
  !******************************************************************





  subroutine get_chi
    integer :: i
    if(mpiID==0)then
       call msg("Get Chi:")
       call get_chi_pm
       call get_chi_dia
       chi(:,:,0:nstep,0:nstep)=chi_pm(:,:,0:nstep,0:nstep)
       do i=0,nstep
          chi(:,:,i,i)=chi(:,:,i,i)+chi_dia(:,:,i)
       enddo
       call MPI_BCAST(chi(:,:,0:nstep,0:nstep),2*2*(nstep+1)**2,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
       call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
    endif
  end subroutine get_chi


  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : evaluate the susceptibility (para-magnetic contribution)
  !+-------------------------------------------------------------------+
  subroutine get_chi_pm
    integer                            :: i,j,ik,ix,iy
    type(vect2D)                       :: Ak,kt,vel
    chi_pm=0.d0
    do ik=1,Lk
       ix=ik2ix(ik)
       iy=ik2iy(ik)
       do i=0,nstep
          Ak = Afield(t(i),Ek)
          kt = kgrid(ix,iy)-Ak
          vel= square_lattice_velocity(kt)
          do j=0,nstep
             chi_pm(1,1,i,j)=chi_pm(1,1,i,j)-2.d0*vel%x*vel%x*wt(ik)*aimag(GkretF(i,j)*Gk%less(j,i))
             chi_pm(1,2,i,j)=chi_pm(1,1,i,j)-2.d0*vel%x*vel%y*wt(ik)*aimag(GkretF(i,j)*Gk%less(j,i))
             chi_pm(2,1,i,j)=chi_pm(1,1,i,j)-2.d0*vel%y*vel%x*wt(ik)*aimag(GkretF(i,j)*Gk%less(j,i))
             chi_pm(2,2,i,j)=chi_pm(1,1,i,j)-2.d0*vel%y*vel%y*wt(ik)*aimag(GkretF(i,j)*Gk%less(j,i))
          enddo
       enddo
    enddo
  end subroutine get_chi_pm



  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : evaluate the susceptibility (dia-magnetic contribution)
  !+-------------------------------------------------------------------+
  subroutine get_chi_dia
    integer       :: i,j,ik,ix,iy
    type(vect2D)  :: vel,kt,Ak
    real(8)       :: eab(2)
    eab=0.d0
    chi_dia=0.d0
    do ik=1,Lk
       ix=ik2ix(ik)
       iy=ik2iy(ik)
       do i=0,nstep
          Ak = Afield(t(i),Ek)
          kt = kgrid(ix,iy)-Ak
          eab(1)=2.d0*ts*cos(kt%x)
          eab(2)=2.d0*ts*cos(kt%y)
          chi_dia(1,1,i)=chi_dia(1,1,i)+2.d0*wt(ik)*eab(1)*xi*Gk%less(i,i)
          chi_dia(2,2,i)=chi_dia(2,2,i)+2.d0*wt(ik)*eab(2)*xi*Gk%less(i,i)
       enddo
    enddo
  end subroutine get_chi_dia



  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : Solve the equilibrium case
  !+-------------------------------------------------------------------+
  subroutine get_equilibrium_localgf()
    integer    :: i,j,ik
    complex(8) :: A,zetan
    real(8)    :: w,n
    complex(8) :: funcM(L),sigmaM(L)
    real(8)    :: funcT(0:L) 
    if(mpiID==0)then
       !Get Sret(w) = FFT(Sret(t-t'))
       forall(i=0:nstep,j=0:nstep) sf%ret%t(i-j)=heaviside(t(i-j))*(Sigma%gtr(i,j)-Sigma%less(i,j))
       sf%ret%t=exa*sf%ret%t ; call fftgf_rt2rw(sf%ret%t,sf%ret%w,nstep) ; sf%ret%w=dt*sf%ret%w

       !Get locGret(w)
       gf%ret%w=zero
       do i=1,2*nstep
          w=wr(i)
          zetan=cmplx(w,eps,8)-sf%ret%w(i) !-eqsbfret(i)
          do ik=1,Lk
             gf%ret%w(i)=gf%ret%w(i)+wt(ik)/(zetan-epsik(ik))
          enddo
       enddo

       !Get locG<(w/t),locG>(w/t)
       gf%less%w=less_component_w(gf%ret%w,wr,beta)
       gf%gtr%w=gtr_component_w(gf%ret%w,wr,beta)
       call fftgf_rw2rt(gf%less%w,gf%less%t,nstep)  ; gf%less%t=exa*fmesh/pi2*gf%less%t
       call fftgf_rw2rt(gf%gtr%w,gf%gtr%t,nstep)    ; gf%gtr%t=exa*fmesh/pi2*gf%gtr%t


       forall(i=0:nstep,j=0:nstep)
          locG%less(i,j) = gf%less%t(i-j)
          locG%gtr(i,j)  = gf%gtr%t(i-j)
          gf%ret%t(i-j) = heaviside(t(i-j))*(locG%gtr(i,j)-locG%less(i,j))
       end forall

       !This is just to get n(k)
       call get_matsubara_gf_from_dos(wr,sf%ret%w,sigmaM,beta)
       do ik=1,Lk
          funcM=zero
          do i=1,L
             w=pi/beta*dble(2*i-1) ; zetan=cmplx(0.d0,w,8) - sigmaM(i)
             funcM(i)=one/(zetan - epsik(ik))
          enddo
          call fftgf_iw2tau(funcM,funcT,beta)
          n=-funcT(L)
          nk(:,ik)=n
       enddo
    endif
    call MPI_BCAST(locG%less,(nstep+1)**2,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BCAST(locG%gtr,(nstep+1)**2,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BCAST(nk,(nstep+1)*Lk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
    call splot('nkVSepsk.ipt',epsik,nk(nstep/2,:),append=TT)
    call splot('locSM_iw.ipt',wm,sigmaM,append=TT)
    call splot("eqG_w.ipt",wr,gf%ret%w,append=TT)
    call splot("eqSigma_w.ipt",wr,sf%ret%w,append=TT)
    return
  end subroutine get_equilibrium_localgf



  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : print out useful information
  !+-------------------------------------------------------------------+
  subroutine print_out_Gloc()
    integer                         :: i,j,ix,iy,ik
    type(vect2D)                    :: Jk,Ak
    type(vect2D),dimension(0:nstep) :: Jloc                   !local Current 
    real(8),dimension(0:nstep)      :: nt,modJloc             !occupation(time)
    if(mpiID==0)then
       call write_kbm_contour_gf(locG,trim(data_dir)//"/locG")
       call splot(trim(data_dir)//"/nk.data",nk(0:,:))
       if(plot3D)call plot_kbm_contour_gf(locG,t(0:),tau(0:),trim(plot_dir)//"/locG")
       call splot("locGlesst0.ipt",t(0:),locG%less(0:,0))
       call splot("locGlmixtau0.ipt",t(0:),locG%lmix(0:,0))
       call splot("locGlmix_t0_tau.ipt",tau(0:),locG%lmix(0,0:))
       call splot("eq_G_tau.ipt",tau(0:),locG%mats(0:),append=.true.)
       call splot("eq_G_iw.ipt",wm,eq_Giw,append=.true.)
       call splot("eq_nk.ipt",epsik,nk(0,:),append=.true.)

       forall(i=0:nstep,j=0:nstep)
          gf%less%t(i-j) = locG%less(i,j)
          gf%gtr%t(i-j)  = locG%gtr(i,j)
          gf%ret%t(i-j)  = heaviside(t(i-j))*(locG%gtr(i,j)-locG%less(i,j))
       end forall
       if(heaviside(0.d0)==1.d0)gf%ret%t(0)=gf%ret%t(0)/2.d0
       call fftgf_rt2rw(gf%ret%t,gf%ret%w,nstep) ;  gf%ret%w=gf%ret%w*dt ; call swap_fftrt2rw(gf%ret%w)
       call splot("locGless_t.ipt",t,gf%less%t,append=TT)
       call splot("locGgtr_t.ipt",t,gf%gtr%t,append=TT)
       call splot("locGret_t.ipt",t,gf%ret%t,append=TT)
       call splot("locGret_realw.ipt",wr,gf%ret%w,append=TT)
       call splot("locDOS.ipt",wr,-aimag(gf%ret%w)/pi,append=TT)

       if(fchi)then
          call splot(trim(data_dir)//"/locChi_11.data",chi(1,1,0:,0:))
          call splot(trim(data_dir)//"/locChi_12.data",chi(1,2,0:,0:))
          call splot(trim(data_dir)//"/locChi_21.data",chi(2,1,0:,0:))
          call splot(trim(data_dir)//"/locChi_22.data",chi(2,2,0:,0:))
          if(plot3D)then
             call splot(trim(plot_dir)//"/locChi_11",t(0:),t(0:),chi(1,1,0:,0:))
             call splot(trim(plot_dir)//"/locChi_12",t(0:),t(0:),chi(1,2,0:,0:))
             call splot(trim(plot_dir)//"/locChi_21",t(0:),t(0:),chi(2,1,0:,0:))
             call splot(trim(plot_dir)//"/locChi_22",t(0:),t(0:),chi(2,2,0:,0:))
          endif
       endif

       forall(i=0:nstep)nt(i)=-xi*locG%less(i,i)

       Jloc=Vzero    
       do ik=1,Lk
          ix=ik2ix(ik);iy=ik2iy(ik)
          do i=0,nstep
             Ak= Afield(t(i),Ek)
             Jk= nk(i,ik)*square_lattice_velocity(kgrid(ix,iy) - Ak)
             Jloc(i) = Jloc(i) +  wt(ik)*Jk
          enddo
       enddo
       call splot("nVStime.ipt",t(0:nstep),2.d0*nt(0:nstep),append=TT)
       modJloc(0:nstep)=modulo(Jloc(0:nstep))
       if(Efield/=0.d0)then
          call splot("JlocVStime.ipt",t(0:nstep),Jloc(0:nstep)%x,Jloc(0:nstep)%y,append=TT)
          call splot("modJlocVStime.ipt",t(0:nstep),modJloc(0:nstep),append=TT)
       endif

    endif
  end subroutine print_out_Gloc



  !******************************************************************
  !******************************************************************
  !******************************************************************

end module KADANOFBAYM
