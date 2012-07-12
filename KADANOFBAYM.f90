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
  USE EQUILIBRIUM
  USE FUNX_NEQ
  implicit none
  private
  !Initial conditions arrays:
  complex(8),allocatable,dimension(:)     :: icGkless 
  real(8),allocatable,dimension(:,:)      :: icGktau
  !k-dependent GF:
  complex(8),allocatable,dimension(:,:)   :: Gkless,Gkgtr,Gklmix,Gkgmix
  real(8),allocatable,dimension(:)        :: Gktau
  !Vector for KB propagation solution: predictor-corrector method
  complex(8),allocatable,dimension(:)     :: Ikless,Ikgtr
  complex(8),allocatable,dimension(:)     :: Ikless0,Ikgtr0
  complex(8),allocatable,dimension(:)     :: Iklmix,Iklmix0
  real(8)                                 :: Ikdiag
  !Auxiliary operators:
  complex(8),allocatable,dimension(:,:)   :: Udelta,Vdelta

  public    :: neq_get_localgf

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the k-sum to get the local Green's functions (R,>,<)
  !+-------------------------------------------------------------------+
  subroutine neq_get_localgf()
    if(solve_wfftw)then
       call get_equilibrium_localgf
    else
       call kadanoff_baym_localgf()
    endif
    if(fchi)call get_chi
    call print_out_Gloc()
    stop
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
  subroutine kadanoff_baym_localgf()
    integer      :: istep,i,j,ik
    integer,save :: loop=1
    complex(8)   :: tmpGless(0:nstep,0:nstep),tmpGgtr(0:nstep,0:nstep)
    complex(8)   :: tmpGlmix(0:nstep,0:Ltau)
    real(8)      :: tmpnk(0:nstep,Lk),tmpGmat(-Ltau:Ltau),Dummy(-Ltau:Ltau)

    !First loop ever, build the Initial Condition
    if(loop==1)then
       call allocate_funx
       call build_ic
       call buildUV
    end if

    call msg("Entering Kadanoff-Baym")

    !Set to Zero loc GF:
    locGless=zero; locGgtr=zero ; locGlmix=zero 
    tmpGless=zero; tmpGgtr =zero; tmpGlmix=zero
    tmpnk=0.d0 ; nk=0.d0
    !=============START K-POINTS LOOP======================
    call start_timer
    do ik=1+mpiID,Lk,mpiSIZE
       Gkless=zero;   Gkgtr=zero

       !Recover Initial Condition for Gk^{<,>}
       call read_ic(ik)
       !======T-STEP LOOP=====================
       do istep=0,nstep-1
          call GFstep(1,ik,istep) !1st-pass
          call GFstep(2,ik,istep) !2nd-pass
       enddo

       tmpGless(0:nstep,0:nstep) =tmpGless(0:nstep,0:nstep) + Gkless(0:nstep,0:nstep)*wt(ik)
       tmpGgtr(0:nstep,0:nstep)  =tmpGgtr(0:nstep,0:nstep)  + Gkgtr(0:nstep,0:nstep)*wt(ik)
       tmpGlmix(0:nstep,0:Ltau)  =tmpGlmix(0:nstep,0:Ltau)  + Gklmix(0:nstep,0:Ltau)*wt(ik)      
       forall(istep=0:nstep)tmpnk(istep,ik)=-xi*Gkless(istep,istep)
       call eta(ik,Lk,unit=6)
    enddo
    call stop_timer
    call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
    !=============END K-POINTS LOOP======================


    call MPI_REDUCE(tmpGless(0:,0:),locGless(0:,0:),(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_REDUCE(tmpGgtr(0:,0:),locGgtr(0:,0:),(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_REDUCE(tmpGlmix(0:,0:),locGlmix(0:,0:),(nstep+1)*(Ltau+1),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)

    !Gloc^>(t',t)= - Gloc^>(t,t')^T
    if(mpiID==0)then
       forall(i=0:nstep,j=0:nstep,i>j)locGless(i,j)=-conjg(locGless(j,i))
       forall(i=0:nstep,j=0:nstep,i<j)locGgtr(i,j) =-conjg(locGgtr(j,i))
       forall(i=0:Ltau)locGgmix(i,:)=-conjg(locGlmix(:,Ltau-i))
    endif
    call MPI_BCAST(locGless(0:,0:),(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BCAST(locGgtr(0:,0:),(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BCAST(locGlmix(0:,0:),(nstep+1)*(Ltau+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)


    call MPI_REDUCE(tmpnk,nk,(nstep+1)*Lk,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(nk,(nstep+1)*Lk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)

    loop=loop+1
  end subroutine kadanoff_baym_localgf



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+----------------------------------------------------------------+
  !PURPOSE  :  allocation of working array
  !+----------------------------------------------------------------+
  subroutine allocate_funx()
    call msg("Allocating KB memory:")
    !Initial conditions for the KBE solution:
    allocate(icGkless(Lk),icGktau(Lk,-Ltau:Ltau))
    !k-dependent Green's functions:
    allocate(Gkless(0:nstep,0:nstep),Gkgtr(0:nstep,0:nstep))
    allocate(Gklmix(0:nstep,0:Ltau),Gkgmix(0:Ltau,0:nstep))
    allocate(Gktau(-Ltau:Ltau))
    !Predictor-corrector solver arrays: store the time-step
    allocate(Ikless(0:nstep),Ikgtr(0:nstep))
    allocate(Ikless0(0:nstep),Ikgtr0(0:nstep))
    allocate(Iklmix(0:Ltau),Iklmix0(0:Ltau))
    !Aux. operators
    allocate(Udelta(Lk,0:nstep),Vdelta(Lk,0:nstep))
    if(fchi)allocate(chi_dia(2,2,0:nstep),chi_pm(2,2,0:nstep,0:nstep))
  end subroutine allocate_funx



  !******************************************************************
  !******************************************************************
  !******************************************************************


  ! !+----------------------------------------------------------------+
  ! !PURPOSE  :  deallocation of working array
  ! !+----------------------------------------------------------------+
  ! subroutine deallocate_funx()
  !   call msg("Deallocating KB memory:")
  !   deallocate(icGkless)
  !   deallocate(Gkless,Gkgtr)      
  !   deallocate(Ikless,Ikgtr)
  !   deallocate(Ikless0,Ikgtr0)
  !   deallocate(Udelta,Vdelta)
  !   if(fchi)deallocate(chi_dia,chi_pm)
  ! end subroutine deallocate_funx



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : Execute the 2pass procedure to solve KB equations:
  !+-------------------------------------------------------------------+
  subroutine GFstep(ips,ik,istep)
    integer :: ips,istep,ik
    integer :: i,j,itau
    integer :: it

    if(ips == 1)then
       !First Pass: get collision integrals up to t=T=istep
       Ikless0 = zero
       Ikgtr0  = zero
       Iklmix0 = zero
       Ikdiag  = 0.d0
       call get_Ikcollision(istep)
       Ikless0  = Ikless
       Ikgtr0   = Ikgtr
       Iklmix0  = iklmix
       Ikdiag   = real(Ikgtr(istep))-real(Ikless(istep))

    elseif(ips==2)then
       !Second Pass: get collision integrals up to t=T+\Delta=istep+1
       call get_Ikcollision(istep+1)
       Ikless   = (Ikless  + Ikless0)/2.d0
       Ikgtr    = (Ikgtr   + Ikgtr0)/2.d0
       Iklmix   = (Iklmix  + Iklmix0)/2.d0
       Ikdiag   = (real(Ikgtr(istep+1))-real(Ikless(istep+1)) + Ikdiag)/2.d0
    endif

    !Evolve the solution of KB equations for all the k-points:
    forall(it=0:istep)
       Gkless(it,istep+1) = Gkless(it,istep)*conjg(Udelta(ik,istep))+Ikless(it)*conjg(Vdelta(ik,istep))
       Gkgtr(istep+1,it)  = Gkgtr(istep,it)*Udelta(ik,istep)+Ikgtr(it)*Vdelta(ik,istep)
    end forall
    Gkgtr(istep+1,istep)=(Gkless(istep,istep)-xi)*Udelta(ik,istep)+Ikgtr(istep)*Vdelta(ik,istep)
    Gkless(istep+1,istep+1)= Gkless(istep,istep)-xi*dt*Ikdiag
    Gkgtr(istep+1,istep+1) = Gkless(istep+1,istep+1)-xi

    forall(itau=0:Ltau)Gklmix(istep+1,itau)=Gklmix(istep,itau)*Udelta(ik,istep)+Iklmix(itau)*Vdelta(ik,istep)
    forall(itau=0:Ltau)Gkgmix(itau,istep+1)=-conjg(Gklmix(istep+1,Ltau-itau))

    !$OMP PARALLEL PRIVATE(i,j)
    !$OMP DO
    do i=0,istep+1!nstep
       do j=0,istep+1!nstep
          if(i>j)Gkless(i,j)=-conjg(Gkless(j,i))
          if(i<j)Gkgtr(i,j)=-conjg(Gkgtr(j,i)) 
       enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    return
  end subroutine GFstep



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
    integer                     :: i,itau,it,itp,count=0
    complex(8)                  :: I1,I2,Ib
    complex(8),dimension(0:Nt)  :: Vadv,Vret,Vless,Vgtr
    complex(8),dimension(0:Ltau):: Vlmix,Vgmix
    complex(8),dimension(0:Nt)  :: Fadv,Fret

    Ikless=zero; Ikgtr=zero; Iklmix=zero
    if(Nt==0)return

    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=0,Nt
       Vless(i)= Sless(i,Nt)    + S0less(i-Nt)
       Vgtr(i) = Sgtr(Nt,i)     + S0gtr(Nt-i)
       Vret(i) = SretF(Nt,i)    + S0retF(Nt-i)
       Vadv(i) = conjg(Vret(i)) 
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=0,Ltau
       Vlmix(i)= Slmix(Nt,i) + S0lmix(Nt,i)
       Vgmix(i)= Sgmix(i,Nt) + S0gmix(i,Nt)
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !Get I^<(t=it,t'=Nt) it:0,...,Nt == t=0,...,T+\Delta
    !I1 = \int_0^{t=it}  G^R*Sigma^<
    !I2 = \int_0^{t'=Nt} G^<*Sigma^A
    !Ib = \int_0^\beta   G^\lmix*S^\gmix
    !===============================================================
    !$OMP PARALLEL SHARED(Ikless) PRIVATE(i,it,I1,I2,Fret,Fadv)
    !$OMP DO
    do it=0,Nt
       forall(i=0:Nt)Fret(i) = GkretF(it,i)
       I1=sum(Fret(0:it)*Vless(0:it))*dt
       I2=sum(Gkless(it,0:Nt)*Vadv(0:Nt))*dt
       Ib=sum(Gklmix(it,0:Ltau)*Vgmix(0:Ltau))*dtau
       Ikless(it)=I1 + I2 - xi*Ib
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    !Get I^>(t=Nt,t'=itp) itp:0,...,Nt == t' = 0,...,T+\Delta=Nt
    !I1 = \int_0^{t=Nt}   S^R*G^>
    !I2 = \int_0^{t'=itp} S^>*G^A
    !Ib = \int_0^\beta    S^\lmix*G^\gmix
    !===============================================================
    !$OMP PARALLEL SHARED(Ikgtr) PRIVATE(i,itp,I1,I2,Fadv)
    !$OMP DO
    do itp=0,Nt
       forall(i=0:Nt)Fadv(i) = conjg(GkretF(itp,i))
       I1=sum(Vret(0:Nt)*Gkgtr(0:Nt,itp))*dt
       I2=sum(Vgtr(0:itp)*Fadv(0:itp))*dt
       Ib=sum(Vlmix(0:Ltau)*Gkgmix(0:Ltau,itp))*dtau
       Ikgtr(itp)=I1 + I2 - xi*Ib
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    !Get I^\lmix(Nt,itau) itau:0:Ltau == tau=-beta:0
    !I1 = \int_0^Nt    S^Ret*G^\rceil
    !Ib = \int_0^\beta S^\lceil* G^\Mats
    !===============================================================
    !$OMP PARALLEL SHARED(Iklmix) PRIVATE(i,itau,I1,Ib)
    !$OMP DO
    do itau=0,Ltau
       I1=sum(Vret(0:Nt)*Gklmix(0:Nt,itau))*dt
       Ib=zero
       do i=0,Ltau
          Ib = Ib+Vlmix(i)*Gktau(i-itau)*dtau
       enddo
       Iklmix(itau)=I1+Ib
       print*,I1,Ib
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine get_Ikcollision



  !******************************************************************
  !******************************************************************
  !******************************************************************


  !-------------------------------------------------------!
  pure function GkretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: GkretF
    GkretF = heaviside(t(i-j))*(Gkgtr(i,j)-Gkless(i,j))
  end function GkretF
  !-------------------------------------------------------!

  !-------------------------------------------------------!
  pure function S0retF(i)      
    integer,intent(in) :: i
    complex(8)         :: S0retF
    S0retF = heaviside(t(i))*(S0gtr(i)-S0less(i))
  end function S0retF
  !-------------------------------------------------------!

  !-------------------------------------------------------!
  pure function SretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: SretF
    SretF = heaviside(t(i-j))*(Sgtr(i,j)-Sless(i,j))
  end function SretF
  !-------------------------------------------------------!





  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : Build the initial condition for the solution of KB equations
  ! G_k^<(0,0)        = xi*G_k^M(tau=0-)
  ! G_k^>(0,0)        = xi*(G_k^<(0,0) - 1.0)
  ! G_k^\lceil(0,tau) = xi*G_k^M(tau<0) = -xi*G_k^M(beta-tau>0)
  ! G_k^\rceil(0,tau) = xi*G_k^M(tau>0)  
  !+-------------------------------------------------------------------+
  subroutine build_ic
    integer    :: ik,i,j
    complex(8) :: funcM(L)
    real(8)    :: n,funcT(0:Ltau),gmtau(-Ltau:Ltau)

    call msg("Building initial conditions:")
    call system("if [ ! -d InitialConditions ]; then mkdir InitialConditions; fi")
    icGkless = xi*eq_Nk
    icGktau  = 0.d0
    gmtau    = 0.d0
    do ik=1,Lk
       funcM=one/(xi*wm - epsik(ik) - eq_Siw)
       call fftgf_iw2tau(funcM,funcT,beta)
       n=-funcT(Ltau) !; nk(:,ik)=n
       icGktau(ik,0:Ltau)=funcT(0:Ltau)
       forall(i=1:Ltau)icGktau(ik,-i)=-funcT(Ltau-i)
       call splot("InitialConditions/icGktau.ipt",tau(0:),icGktau(ik,0:),append=TT)
       gmtau = gmtau + icGktau(ik,:)*wt(ik)
    enddo
    call splot("InitialConditions/icGklessVSepsik.ipt",epsik(1:Lk),dimag(icGkless(1:Lk)))
    forall(i=0:Ltau,j=0:Ltau)locGmat(i,j)=gmtau(j-i)
  end subroutine build_ic



  !******************************************************************
  !******************************************************************
  !******************************************************************




  subroutine read_ic(ik)
    integer :: ik,i
    Gkless(0,0) = icGkless(ik)
    Gkgtr(0,0)  = icGkless(ik)-xi
    forall(i=0:Ltau)
       Gklmix(0,i)=-xi*icGktau(ik,Ltau-i)
       Gkgmix(i,0)= xi*icGktau(ik,i)
    end forall
    Gktau(-Ltau:Ltau)=icGktau(ik,-Ltau:Ltau)
    return
  end subroutine read_ic



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
    complex(8)   :: Hbar
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




  function UdeltaF(ik,istep) 
    integer,intent(in)    :: ik,istep
    complex(8) :: UdeltaF
    complex(8) :: arg
    arg=Hbar(ik,istep)
    UdeltaF=exp(-xi*arg*dt)
  end function UdeltaF



  !******************************************************************
  !******************************************************************
  !******************************************************************



  function VdeltaF(ik,istep)
    integer,intent(in)    :: ik,istep
    complex(8) :: VdeltaF
    complex(8) :: arg
    arg=Hbar(ik,istep)
    VdeltaF=exp(-xi*arg*dt)
    if(abs(arg*dt) <= 1.d-5)then
       VdeltaF=-xi*dt
    else
       VdeltaF=(VdeltaF-1.d0)/arg
    endif
  end function VdeltaF



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
             chi_pm(1,1,i,j)=chi_pm(1,1,i,j)-2.d0*vel%x*vel%x*wt(ik)*aimag(GkretF(i,j)*Gkless(j,i))
             chi_pm(1,2,i,j)=chi_pm(1,1,i,j)-2.d0*vel%x*vel%y*wt(ik)*aimag(GkretF(i,j)*Gkless(j,i))
             chi_pm(2,1,i,j)=chi_pm(1,1,i,j)-2.d0*vel%y*vel%x*wt(ik)*aimag(GkretF(i,j)*Gkless(j,i))
             chi_pm(2,2,i,j)=chi_pm(1,1,i,j)-2.d0*vel%y*vel%y*wt(ik)*aimag(GkretF(i,j)*Gkless(j,i))
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
          eab(1)=2*ts*cos(kt%x)
          eab(2)=2*ts*cos(kt%y)
          chi_dia(1,1,i)=chi_dia(1,1,i)+2.d0*wt(ik)*eab(1)*xi*Gkless(i,i)
          chi_dia(2,2,i)=chi_dia(2,2,i)+2.d0*wt(ik)*eab(2)*xi*Gkless(i,i)
       enddo
    enddo
  end subroutine get_chi_dia



  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+-------------------------------------------------------------------+
  !PURPOSE  : print out useful information
  !+-------------------------------------------------------------------+
  subroutine print_out_Gloc()
    integer :: i,j
    if(mpiID==0)then
       call system("if [ ! -d GLOC ]; then mkdir GLOC; fi")
       call splot("locGless.data",locGless(0:,0:))
       call splot("locGgtr.data",locGgtr(0:,0:))
       call splot("locGlmix.data",locGlmix(0:,0:))
       call splot("locGmat.data",locGmat(0:,0:))
       call splot("nk.data",nk(0:,:))

       forall(i=0:nstep,j=0:nstep)
          gf%less%t(i-j) = locGless(i,j)
          gf%gtr%t(i-j)  = locGgtr(i,j)
          gf%ret%t(i-j)  = heaviside(t(i-j))*(locGgtr(i,j)-locGless(i,j))
       end forall
       if(heaviside(0.d0)==1.d0)gf%ret%t(0)=gf%ret%t(0)/2.d0
       call fftgf_rt2rw(gf%ret%t,gf%ret%w,nstep) ;  gf%ret%w=gf%ret%w*dt ; call swap_fftrt2rw(gf%ret%w)
       call splot("locGless_t.ipt",t,gf%less%t,append=TT)
       call splot("locGgtr_t.ipt",t,gf%gtr%t,append=TT)
       call splot("locGret_t.ipt",t,gf%ret%t,append=TT)
       call splot("locGret_realw.ipt",wr,gf%ret%w,append=TT)
       call splot("locDOS.ipt",wr,-aimag(gf%ret%w)/pi,append=TT)
       if(fchi)then
          call splot("locChi_11.data",chi(1,1,0:nstep,0:nstep))
          call splot("locChi_12.data",chi(1,2,0:nstep,0:nstep))
          call splot("locChi_21.data",chi(2,1,0:nstep,0:nstep))
          call splot("locChi_22.data",chi(2,2,0:nstep,0:nstep))
       endif
       call splot("locGless_t_t.3d",t(0:),t(0:),locGless(0:,0:))
       call splot("locGlmix_t_tau.3d",t(0:),tau(0:),locGlmix(0:,0:))
       call splot("locGmat_tau_tau.3d",tau(0:),tau(0:),locGmat(0:,0:))
       call system("mv -vf *locG*.3d GLOC/")
    endif



  end subroutine print_out_Gloc



  !******************************************************************
  !******************************************************************
  !******************************************************************

end module KADANOFBAYM
