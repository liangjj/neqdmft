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
  complex(8),allocatable,dimension(:,:,:)            :: Gk
  !Vector for KB propagation solution
  complex(8),allocatable,dimension(:,:)              :: Ikless,Ikgtr
  complex(8),allocatable,dimension(:,:)              :: Ikless0,Ikgtr0
  real(8),allocatable,dimension(:)                   :: Ikdiag

  public                                             :: neq_get_localgf

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the k-sum to get the local Green's functions (R,>,<)
  !+-------------------------------------------------------------------+
  subroutine neq_get_localgf()
    call kadanoff_baym_to_localgf()
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
    integer :: i,j,ik,istep
    call msg("Entering Kadanoff-Baym")
    call allocate_funx
    locG   = zero
    nk     = 0.d0

    !=============START T-STEP LOOP======================
    call start_timer
    call setup_initial_conditions
    do istep=0,nstep-1

       call step_keldysh_contour_gf(istep)

       call eta(istep,nstep-1)
    enddo
    call stop_timer
    !=============END T-STEP LOOP======================    

    forall(i=0:Nstep,j=0:Nstep,i>j)locG%less(i,j)=-conjg(locG%less(j,i))
    forall(i=0:Nstep,j=0:Nstep,i<j)locG%gtr(i,j) =-conjg(locG%gtr(j,i))
    forall(i=0:Nstep)locG%gtr(i,i)=locG%less(i,i)-xi
    if(fchi)call get_chi
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
    !Allocate k-dependent GF:
    allocate(Gk(Lk,0:Nstep,0:Nstep));Gk=zero
    !Predictor-corrector solver arrays: store the time-step
    allocate(Ikless(Lk,0:nstep),Ikgtr(Lk,0:nstep))
    allocate(Ikless0(Lk,0:nstep),Ikgtr0(Lk,0:nstep))
    allocate(Ikdiag(Lk))
    !Chi
    if(fchi)allocate(chi_dia(2,2,0:nstep),chi_pm(2,2,0:nstep,0:nstep))
  end subroutine allocate_funx

  subroutine deallocate_funx()
    call msg("Deallocating KB memory:")
    deallocate(Gk)
    deallocate(Ikless,Ikgtr)
    deallocate(Ikless0,Ikgtr0)
    deallocate(Ikdiag)
    if(fchi)deallocate(chi_dia,chi_pm)
  end subroutine deallocate_funx


  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : Build up the initial conditions at T=0
  !+-------------------------------------------------------------------+
  subroutine setup_initial_conditions
    integer :: ik
    do ik=1,Lk
       Gk(ik,0,0)= xi*eq_nk(ik)
       locG%less(0,0) = locG%less(0,0) + Gk(ik,0,0)*wt(ik)
    enddo
    locG%gtr(0,0)  = locG%less(0,0)-xi
  end subroutine setup_initial_conditions



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : Execute the 2pass procedure to solve KB equations:
  !+-------------------------------------------------------------------+
  subroutine step_keldysh_contour_gf(Nt)
    integer :: Nt,ik
    integer :: i,j,itau
    integer :: it
    real(8) :: nloc
    nloc=-xi*locG%less(Nt,Nt)
    SigmaHF(Nt) = U*(nloc-0.5d0)
    !Start building the Gk(T+\D,T+\D) with the two-step (predictor-corrector) method
    !First step: get collision integrals up to t=T=Nt
    nloc=0.d0
    do ik=1,Lk
       call get_Ikcollision(ik,Nt)
       Ikless0(ik,:) = Ikless(ik,:)
       Ikgtr0(ik,:)  = Ikgtr(ik,:)
       Ikdiag(ik)    =-2.d0*dreal(Ikless(ik,Nt))
       !
       call evolve_gk(ik,Nt)
       !
       !locG%less(Nt+1,Nt+1) = locG%less(Nt+1,Nt+1) + Gk(ik,Nt+1,Nt+1)*wt(ik)
       nloc = nloc - xi*Gk(ik,Nt+1,Nt+1)*wt(ik)
    enddo

    SigmaHF(Nt+1) = (U*(nloc-0.5d0) + SigmaHF(Nt))/2.d0
    !Second step: get collision integrals up to t=T+\Delta=Nt+1
    do ik=1,Lk
       call get_Ikcollision(ik,Nt+1)
       Ikless(ik,:)  = (Ikless(ik,:)  + Ikless0(ik,:))/2.d0
       Ikgtr(ik,:)   = (Ikgtr(ik,:)   + Ikgtr0(ik,:))/2.d0
       Ikdiag(ik)    = (-2.d0*dreal(Ikless(ik,Nt+1)) + Ikdiag(ik))/2.d0
       !
       call evolve_gk(ik,Nt)
       !
       !sum over k-point
       locG%less(0:Nt+1,Nt+1) = locG%less(0:Nt+1,Nt+1) + Gk(ik,0:Nt+1,Nt+1)*wt(ik) !i<=j
       locG%gtr(Nt+1,0:Nt)    = locG%gtr(Nt+1,0:Nt)    + Gk(ik,Nt+1,0:Nt)*wt(ik)   !i>j
       nk(Nt+1,ik)            = -xi*Gk(ik,Nt+1,Nt+1)
    enddo

  end subroutine step_keldysh_contour_gf


  subroutine evolve_gk(ik,Nt)
    integer :: it,ik,Nt
    !Evolve the solution of KB equations G_k(T)--> G_k(T+\D):
    !the upper triangular terms evolve G_k^<(t,T+\D)
    !the lower triangular terms evolve G_k^>(T+\D,t)
    do it=0,Nt
       Gk(ik,it,Nt+1) = Gk(ik,it,Nt)*conjg(Udelta(ik,Nt))-Ikless(ik,it)*conjg(Vdelta(ik,Nt))
       Gk(ik,Nt+1,it) = Gk(ik,Nt,it)*Udelta(ik,Nt)-Ikgtr(ik,it)*Vdelta(ik,Nt)
    end do
    !the evolution of the last term for G_k^>(T+\D,T) needs to start from G_k(T,T):
    Gk(ik,Nt+1,Nt)=(Gk(ik,Nt,Nt)-xi)*Udelta(ik,Nt)-Ikgtr(ik,Nt)*Vdelta(ik,Nt)
    !the evolution along the diagonal is only for the G_k^<(T+\D,T+\D)
    Gk(ik,Nt+1,Nt+1)= Gk(ik,Nt,Nt)-xi*dt*Ikdiag(ik) !diagonal part is Gk^<(t,t)
    !
    !this is not needed anymore
    ! Gk%gtr(Nt+1,Nt+1) = Gk%less(Nt+1,Nt+1)-xi
    ! do i=0,Nt+1
    !    Gk%less(Nt+1,i)=-conjg(Gk%less(i,Nt+1))
    !    Gk%gtr(i,Nt+1) =-conjg(Gk%gtr(Nt+1,i))
    ! enddo
  end subroutine evolve_gk

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
  subroutine get_Ikcollision(ik,Nt)
    integer,intent(in)              :: ik,Nt
    integer                         :: i,j,itau,it,itp
    complex(8)                      :: I1,I2,Ib
    complex(8),dimension(0:Nt)      :: Vadv,Vret,Vless,Vgtr
    complex(8),dimension(0:Nt)      :: Gadv,Gret
    complex(8),dimension(0:Nt,0:Nt) :: Gless,Ggtr

    Ikless=zero; Ikgtr=zero
    if(Nt==0)return

    do i=0,Nt
       Vless(i)= Sigma%less(i,Nt)    + S0%less(i,Nt)
       Vgtr(i) = Sigma%gtr(Nt,i)     + S0%gtr(Nt,i)
       Vret(i) = Sret(Nt,i)         + S0ret(Nt,i)
       Vadv(i) = conjg(Vret(i))
    end do

    do i=0,Nt
       do j=0,Nt
          if(i<j)then           !upper triangle: G^<\=Gk, G^>=-Gk^+
             Gless(i,j)=       Gk(ik,i,j)
             Ggtr(i,j) =-conjg(Gk(ik,j,i))
          elseif(i==j)then
             Gless(i,i)=       Gk(ik,i,i)
             Ggtr(i,i) =       Gk(ik,i,i)-xi
          else
             Gless(i,j)=-conjg(Gk(ik,j,i))
             Ggtr(i,j) =       Gk(ik,i,j)
          endif
       enddo
    enddo


    select case(int_method)
    case default             !trapz
       !Get I^<(t=it,t'=Nt) it:0,...,Nt == t=0,...,T+\Delta
       !I1 = \int_0^{t=it}  G^R*Sigma^<
       !I2 = \int_0^{t'=Nt} G^<*Sigma^A
       do it=0,Nt
          forall(i=0:it)Gret(i) = Fret(it,i)
          I1=trapz(dt,Gret(0:it)*Vless(0:it))
          I2=trapz(dt,Gless(it,0:Nt)*Vadv(0:Nt))
          Ikless(ik,it)=I1 + I2
       enddo
       !
       !Get I^>(t=Nt,t'=itp) itp:0,...,Nt == t' = 0,...,T+\Delta=Nt
       !I1 = \int_0^{t=Nt}   S^R*G^>
       !I2 = \int_0^{t'=itp} S^>*G^A
       do itp=0,Nt
          forall(i=0:itp)Gadv(i) = conjg(Fret(itp,i))
          I1=trapz(dt,Vret(0:Nt)*Ggtr(0:Nt,itp))
          I2=trapz(dt,Vgtr(0:itp)*Gadv(0:itp))
          Ikgtr(ik,itp)=I1 + I2
       enddo

    case ("simps")
       do it=0,Nt
          forall(i=0:it)Gret(i) = Fret(it,i)
          I1=simps(dt,Gret(0:it)*Vless(0:it))
          I2=simps(dt,Gless(it,0:Nt)*Vadv(0:Nt))
          Ikless(ik,it)=I1 + I2
       enddo
       !
       do itp=0,Nt
          forall(i=0:itp)Gadv(i) = conjg(Fret(itp,i))
          I1=simps(dt,Vret(0:Nt)*Ggtr(0:Nt,itp))
          I2=simps(dt,Vgtr(0:itp)*Gadv(0:itp))
          Ikgtr(ik,itp)=I1 + I2
       enddo

    case ("rect")
       do it=0,Nt
          forall(i=0:it)Gret(i) = Fret(it,i)
          I1=sum(Gret(0:it)*Vless(0:it))*dt
          I2=sum(Gless(it,0:Nt)*Vadv(0:Nt))*dt
          Ikless(ik,it)=I1 + I2
       enddo
       !
       do itp=0,Nt
          forall(i=0:itp)Gadv(i) = conjg(Fret(itp,i))
          I1=sum(Vret(0:Nt)*Ggtr(0:Nt,itp))*dt
          I2=sum(Vgtr(0:itp)*Gadv(0:itp))*dt
          Ikgtr(ik,itp)=I1 + I2
       enddo
    end select

  contains
    pure function Fret(i,j)      
      integer,intent(in) :: i,j
      complex(8)         :: Fret
      Fret = heaviside(t(i)-t(j))*(Ggtr(i,j)-Gless(i,j))
    end function Fret

    pure function S0ret(i,j)
      integer,intent(in) :: i,j
      complex(8)         :: S0ret
      S0ret = heaviside(t(i)-t(j))*(S0%gtr(i,j)-S0%less(i,j))
    end function S0ret

    pure function Sret(i,j)      
      integer,intent(in) :: i,j
      complex(8)         :: Sret
      Sret = heaviside(t(i)-t(j))*(Sigma%gtr(i,j)-Sigma%less(i,j))
    end function Sret
  end subroutine get_Ikcollision


  !******************************************************************
  !******************************************************************
  !******************************************************************


  function Udelta(ik,Nt)  result(UdeltaF)
    integer,intent(in)    :: ik,Nt
    complex(8) :: UdeltaF
    real(8) :: arg
    arg=Hbar(ik,Nt)
    UdeltaF=exp(-xi*arg*dt)
  end function Udelta

  function Vdelta(ik,Nt) result(VdeltaF)
    integer,intent(in)    :: ik,Nt
    complex(8) :: VdeltaF
    real(8) :: arg
    arg=Hbar(ik,Nt)
    VdeltaF=exp(-xi*arg*dt)
    if(abs(arg*dt) <= 1.d-9)then
       VdeltaF=xi*dt
    else
       VdeltaF=(1.d0-VdeltaF)/arg
    endif
  end function Vdelta



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : This is the time-dependent hamiltonian
  ! a more general code would require this quantity to be definied 
  ! externally as an array for every k-point (ik) and time (Nt) but
  ! I am lazy...
  !+-------------------------------------------------------------------+
  function Hbar(ik,Nt)
    integer,intent(in) :: ik,Nt  
    integer      :: i,j
    real(8)      :: Hbar
    real(8)      :: tbar
    type(vect2D) :: kt,Ak
    tbar=t(Nt) + dt/2.d0
    i=ik2ix(ik)
    j=ik2iy(ik)
    Ak=Afield(tbar,Ek)
    kt=kgrid(i,j) - Ak
    Hbar=square_lattice_dispersion(kt)+SigmaHF(Nt)-xmu
  end function Hbar


  !******************************************************************
  !******************************************************************
  !******************************************************************





  subroutine get_chi
    integer :: i
    if(mpiID==0)then
       call msg("Get Chi:")
       ! call get_chi_pm
       ! call get_chi_dia
       ! chi(:,:,0:nstep,0:nstep)=chi_pm(:,:,0:nstep,0:nstep)
       ! do i=0,nstep
       !    chi(:,:,i,i)=chi(:,:,i,i)+chi_dia(:,:,i)
       ! enddo
       call msg("Sorry no Chi...")
    endif
  end subroutine get_chi


  !******************************************************************
  !******************************************************************
  !******************************************************************




  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : evaluate the susceptibility (para-magnetic contribution)
  ! !+-------------------------------------------------------------------+
  ! subroutine get_chi_pm
  !   integer                               :: i,j,ik,ix,iy
  !   type(vect2D)                          :: Ak,kt,vel


  !   chi_pm=0.d0
  !   do ik=1,Lk
  !      ix=ik2ix(ik)
  !      iy=ik2iy(ik)
  !      do i=0,nstep
  !         Ak = Afield(t(i),Ek)
  !         kt = kgrid(ix,iy)-Ak
  !         vel= square_lattice_velocity(kt)
  !         do j=0,nstep
  !            chi_pm(1,1,i,j)=chi_pm(1,1,i,j)-2.d0*vel%x*vel%x*wt(ik)*aimag(GkretF(i,j)*Gk%less(j,i))
  !            chi_pm(1,2,i,j)=chi_pm(1,1,i,j)-2.d0*vel%x*vel%y*wt(ik)*aimag(GkretF(i,j)*Gk%less(j,i))
  !            chi_pm(2,1,i,j)=chi_pm(1,1,i,j)-2.d0*vel%y*vel%x*wt(ik)*aimag(GkretF(i,j)*Gk%less(j,i))
  !            chi_pm(2,2,i,j)=chi_pm(1,1,i,j)-2.d0*vel%y*vel%y*wt(ik)*aimag(GkretF(i,j)*Gk%less(j,i))
  !         enddo
  !      enddo
  !   enddo
  ! end subroutine get_chi_pm



  ! !******************************************************************
  ! !******************************************************************
  ! !******************************************************************




  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : evaluate the susceptibility (dia-magnetic contribution)
  ! !+-------------------------------------------------------------------+
  ! subroutine get_chi_dia
  !   integer       :: i,j,ik,ix,iy
  !   type(vect2D)  :: vel,kt,Ak
  !   real(8)       :: eab(2)
  !   eab=0.d0
  !   chi_dia=0.d0
  !   do ik=1,Lk
  !      ix=ik2ix(ik)
  !      iy=ik2iy(ik)
  !      do i=0,nstep
  !         Ak = Afield(t(i),Ek)
  !         kt = kgrid(ix,iy)-Ak
  !         eab(1)=2.d0*ts*cos(kt%x)
  !         eab(2)=2.d0*ts*cos(kt%y)
  !         chi_dia(1,1,i)=chi_dia(1,1,i)+2.d0*wt(ik)*eab(1)*xi*Gk%less(i,i)
  !         chi_dia(2,2,i)=chi_dia(2,2,i)+2.d0*wt(ik)*eab(2)*xi*Gk%less(i,i)
  !      enddo
  !   enddo
  ! end subroutine get_chi_dia



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

       call write_keldysh_contour_gf(locG,trim(data_dir)//"/locG")
       call splot(trim(data_dir)//"/nk.data",nk(0:,:))

       if(plot3D)then
          call plot_keldysh_contour_gf(locG,t(0:),trim(plot_dir)//"/locG")
       end if

       forall(i=0:nstep,j=0:nstep)
          gf%less%t(i-j) = locG%less(i,j)
          gf%gtr%t(i-j)  = locG%gtr(i,j)
          gf%ret%t(i-j)  = heaviside(t(i-j))*(locG%gtr(i,j)-locG%less(i,j))
       end forall
       !if(heaviside(0.d0)==1.d0)gf%ret%t(0)=gf%ret%t(0)/2.d0
       call fftgf_rt2rw(gf%ret%t,gf%ret%w,nstep) ;  gf%ret%w=gf%ret%w*dt ; call swap_fftrt2rw(gf%ret%w)
       call splot("locGless_t.ipt",t,gf%less%t,append=TT)
       call splot("locGgtr_t.ipt",t,gf%gtr%t,append=TT)
       call splot("locGret_t.ipt",t,gf%ret%t,append=TT)
       call splot("locGret_realw.ipt",wr,gf%ret%w,append=TT)
       call splot("locDOS.ipt",wr,-aimag(gf%ret%w)/pi,append=TT)
       call splot("SigmaHF_t.ipt",t(0:),SigmaHF(0:),append=TT)
       ! if(fchi)then
       !    call splot(trim(data_dir)//"/locChi_11.data",chi(1,1,0:,0:))
       !    call splot(trim(data_dir)//"/locChi_12.data",chi(1,2,0:,0:))
       !    call splot(trim(data_dir)//"/locChi_21.data",chi(2,1,0:,0:))
       !    call splot(trim(data_dir)//"/locChi_22.data",chi(2,2,0:,0:))
       !    if(plot3D)then
       !       call splot(trim(plot_dir)//"/locChi_11",t(0:),t(0:),chi(1,1,0:,0:))
       !       call splot(trim(plot_dir)//"/locChi_12",t(0:),t(0:),chi(1,2,0:,0:))
       !       call splot(trim(plot_dir)//"/locChi_21",t(0:),t(0:),chi(2,1,0:,0:))
       !       call splot(trim(plot_dir)//"/locChi_22",t(0:),t(0:),chi(2,2,0:,0:))
       !    endif
       ! endif

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
