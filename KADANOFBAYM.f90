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
  !k-dependent GF:
  type(keldysh_contour_gf)                :: Gk
  !Vector for KB propagation solution
  complex(8),allocatable,dimension(:)     :: Ikless,Ikgtr
  complex(8),allocatable,dimension(:)     :: Ikless0,Ikgtr0
  real(8)                                 :: Ikdiag
  !Auxiliary operators:
  complex(8),allocatable,dimension(:,:)   :: Udelta,Vdelta

  public                                  :: neq_get_localgf

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
    integer              :: istep,i,j,ik
    integer,save         :: loop=1
    real(8),allocatable  :: tmpnk(:,:)
    type(keldysh_contour_gf) :: tmpG,tmpG1,tmpG2

    !First loop ever, build the Initial Condition
    if(loop==1)then
       call allocate_funx
       call build_ic
       call buildUV
    end if
    loop=loop+1

    call msg("Entering Kadanoff-Baym")

    !Set to Zero loc GF:
    locG=zero ; locG1=zero ; locG2=zero

    !Tmp array for MPI storage, set to zero
    call allocate_keldysh_contour_gf(tmpG,Nstep)
    tmpG=zero
    if(volterra)then
       call allocate_keldysh_contour_gf(tmpG1,Nstep)
       call allocate_keldysh_contour_gf(tmpG2,Nstep)
       tmpG1=zero ; tmpG2=zero
    endif

    !set to zero n(k,t)
    allocate(tmpnk(0:nstep,Lk))
    tmpnk=0.d0 ; nk=0.d0


    !=============START K-POINTS LOOP======================
    call start_timer
    call allocate_keldysh_contour_gf(Gk,Nstep)
    do ik=1+mpiID,Lk,mpiSIZE
       Gk=zero

       !Recover Initial Conditions
       call read_ic(ik)

       !======T-STEP LOOP=====================
       do istep=0,nstep-1
          call GFstep(1,ik,istep) !1st-pass
          call GFstep(2,ik,istep) !2nd-pass
       enddo

       !sum over k-point
       call keldysh_contour_gf_sum(tmpG,Gk,wt(ik))

       if(volterra)then
          !get G1=sum_k h(k,t)G(t,t')
          do i=0,nstep
             tmpG1%less(i,0:) = tmpG1%less(i,0:) + Hkt(ik,i)*Gk%less(i,0:)*wt(ik)
             tmpG1%gtr(i,0:)  = tmpG1%gtr(i,0:)  + Hkt(ik,i)*Gk%gtr(i,0:)*wt(ik)
          end do

          !get G2=sum_k h(k,t)G(t,t')h(k,t')
          do i=0,nstep
             do j=0,nstep
                tmpG2%less(i,j) = tmpG2%less(i,j) + Hkt(ik,i)*Gk%less(i,j)*Hkt(ik,j)*wt(ik)
                tmpG2%gtr(i,j)  = tmpG2%gtr(i,j)  + Hkt(ik,i)*Gk%gtr(i,j)*Hkt(ik,j)*wt(ik)
             enddo
          enddo
       endif

       forall(istep=0:nstep)tmpnk(istep,ik)=-xi*Gk%less(istep,istep)
       call eta(ik,Lk,unit=6)
    enddo
    call stop_timer
    call deallocate_keldysh_contour_gf(Gk)
    call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
    !=============END K-POINTS LOOP======================


    !Reduce Contour GF to master:
    call MPI_REDUCE_keldysh_contour_gf(tmpG,locG)
    call deallocate_keldysh_contour_gf(tmpG)
    if(volterra)then
       call MPI_REDUCE_keldysh_contour_gf(tmpG1,locG1)
       call MPI_REDUCE_keldysh_contour_gf(tmpG2,locG2)
       call deallocate_keldysh_contour_gf(tmpG1)
       call deallocate_keldysh_contour_gf(tmpG2)
    endif

    !Bcast the local contour GF to every node
    call MPI_BCAST_keldysh_contour_gf(locG)
    if(volterra)then
       call MPI_BCAST_keldysh_contour_gf(locG1)
       call MPI_BCAST_keldysh_contour_gf(locG2)
    endif

    call MPI_REDUCE(tmpnk,nk,(nstep+1)*Lk,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(nk,(nstep+1)*Lk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
    deallocate(tmpnk)

    call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
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
    allocate(icGkless(Lk))

    !Predictor-corrector solver arrays: store the time-step
    allocate(Ikless(0:nstep),Ikgtr(0:nstep))
    allocate(Ikless0(0:nstep),Ikgtr0(0:nstep))

    !Aux. operators
    allocate(Udelta(Lk,0:nstep),Vdelta(Lk,0:nstep))

    !Chi
    if(fchi)allocate(chi_dia(2,2,0:nstep),chi_pm(2,2,0:nstep,0:nstep))
  end subroutine allocate_funx



  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+----------------------------------------------------------------+
  !PURPOSE  :  deallocation of working array
  !+----------------------------------------------------------------+
  subroutine deallocate_funx()
    call msg("Deallocating KB memory:")
    deallocate(icGkless)
    deallocate(Ikless,Ikgtr)
    deallocate(Ikless0,Ikgtr0)
    deallocate(Udelta,Vdelta)
    if(fchi)deallocate(chi_dia,chi_pm)
  end subroutine deallocate_funx



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

    select case(ips)
    case(1)
       !First Pass: get collision integrals up to t=T=istep
       Ikless0 = zero
       Ikgtr0  = zero
       Ikdiag  = 0.d0
       call get_Ikcollision(istep)
       Ikless0  = Ikless
       Ikgtr0   = Ikgtr
       Ikdiag   = real(Ikgtr(istep))-real(Ikless(istep))

    case(2)
       !Second Pass: get collision integrals up to t=T+\Delta=istep+1
       call get_Ikcollision(istep+1)
       Ikless   = (Ikless  + Ikless0)/2.d0
       Ikgtr    = (Ikgtr   + Ikgtr0)/2.d0
       Ikdiag   = (real(Ikgtr(istep+1))-real(Ikless(istep+1)) + Ikdiag)/2.d0

    end select

    !Evolve the solution of KB equations for all the k-points:
    forall(it=0:istep)
       Gk%less(it,istep+1) = Gk%less(it,istep)*conjg(Udelta(ik,istep))+Ikless(it)*conjg(Vdelta(ik,istep))
       Gk%gtr(istep+1,it)  = Gk%gtr(istep,it)*Udelta(ik,istep)+Ikgtr(it)*Vdelta(ik,istep)
    end forall
    Gk%gtr(istep+1,istep)=(Gk%less(istep,istep)-xi)*Udelta(ik,istep)+Ikgtr(istep)*Vdelta(ik,istep)
    Gk%less(istep+1,istep+1)= Gk%less(istep,istep)-xi*dt*Ikdiag
    Gk%gtr(istep+1,istep+1) = Gk%less(istep+1,istep+1)-xi


    !$OMP PARALLEL PRIVATE(i,j)
    !$OMP DO
    do i=0,istep+1
       do j=0,istep+1
          if(i>j)Gk%less(i,j)=-conjg(Gk%less(j,i))
          if(i<j)Gk%gtr(i,j)=-conjg(Gk%gtr(j,i)) 
       enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

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
    integer                     :: i,itau,it,itp
    complex(8)                  :: I1,I2,Ib
    complex(8),dimension(0:Nt)  :: Vadv,Vret,Vless,Vgtr
    complex(8),dimension(0:Nt)  :: Fadv,Fret

    Ikless=zero; Ikgtr=zero
    if(Nt==0)return

    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=0,Nt
       Vless(i)= Sigma%less(i,Nt)    + S0%less(i,Nt)
       Vgtr(i) = Sigma%gtr(Nt,i)     + S0%gtr(Nt,i)
       Vret(i) = SretF(Nt,i)         + S0retF(Nt,i)
       Vadv(i) = conjg(Vret(i))
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !Get I^<(t=it,t'=Nt) it:0,...,Nt == t=0,...,T+\Delta
    !I1 = \int_0^{t=it}  G^R*Sigma^<
    !I2 = \int_0^{t'=Nt} G^<*Sigma^A
    !Ib = \int_0^\beta   G^\rceil*S^\rceil
    !===============================================================
    !$OMP PARALLEL SHARED(Ikless) PRIVATE(i,it,I1,I2,Fret,Fadv)
    !$OMP DO
    do it=0,Nt
       forall(i=0:Nt)Fret(i) = GkretF(it,i)
       I1=sum(Fret(0:it)*Vless(0:it))*dt
       I2=sum(Gk%less(it,0:Nt)*Vadv(0:Nt))*dt
       Ikless(it)=I1 + I2
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    !Get I^>(t=Nt,t'=itp) itp:0,...,Nt == t' = 0,...,T+\Delta=Nt
    !I1 = \int_0^{t=Nt}   S^R*G^>
    !I2 = \int_0^{t'=itp} S^>*G^A
    !Ib = \int_0^\beta    S^\lceil*G^\rceil
    !===============================================================
    !$OMP PARALLEL SHARED(Ikgtr) PRIVATE(i,itp,I1,I2,Fadv)
    !$OMP DO
    do itp=0,Nt
       forall(i=0:Nt)Fadv(i) = conjg(GkretF(itp,i))
       I1=sum(Vret(0:Nt)*Gk%gtr(0:Nt,itp))*dt
       I2=sum(Vgtr(0:itp)*Fadv(0:itp))*dt
       Ikgtr(itp)=I1 + I2
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
    GkretF = heaviside(t(i)-t(j))*(Gk%gtr(i,j)-Gk%less(i,j))
  end function GkretF
  !-------------------------------------------------------!

  !-------------------------------------------------------!
  pure function S0retF(i,j)   
    integer,intent(in) :: i,j
    complex(8)         :: S0retF
    S0retF = heaviside(t(i)-t(j))*(S0%gtr(i,j)-S0%less(i,j))
  end function S0retF
  !-------------------------------------------------------!

  !-------------------------------------------------------!
  pure function SretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: SretF
    SretF = heaviside(t(i)-t(j))*(Sigma%gtr(i,j)-Sigma%less(i,j))
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
    integer :: ik,i
    call msg("Building initial conditions:")
    call system("if [ ! -d InitialConditions ]; then mkdir InitialConditions; fi")
    icGkless = xi*eq_Nk
    call splot("InitialConditions/icGklessVSepsik.ipt",epsik(1:Lk),dimag(icGkless(1:Lk)))
  end subroutine build_ic



  !******************************************************************
  !******************************************************************
  !******************************************************************




  subroutine read_ic(ik)
    integer :: ik
    Gk%less(0,0)=icGkless(ik)
    Gk%gtr(0,0) =icGkless(ik)-xi
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
  function Hkt(ik,istep)
    integer,intent(in) :: ik,istep  
    integer      :: i,j
    complex(8)   :: Hkt
    type(vect2D) :: kt,Ak
    i=ik2ix(ik)
    j=ik2iy(ik)
    Ak=Afield(t(istep),Ek)
    kt=kgrid(i,j) - Ak
    Hkt=square_lattice_dispersion(kt)
  end function Hkt


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
          eab(1)=2*ts*cos(kt%x)
          eab(2)=2*ts*cos(kt%y)
          chi_dia(1,1,i)=chi_dia(1,1,i)+2.d0*wt(ik)*eab(1)*xi*Gk%less(i,i)
          chi_dia(2,2,i)=chi_dia(2,2,i)+2.d0*wt(ik)*eab(2)*xi*Gk%less(i,i)
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

       call write_keldysh_contour_gf(locG,reg_filename(data_dir)//"/locG")
       if(volterra)then
          call write_keldysh_contour_gf(locG1,reg_filename(data_dir)//"/locG1")
          call write_keldysh_contour_gf(locG2,reg_filename(data_dir)//"/locG2")
       endif
       call splot(reg_filename(data_dir)//"/nk.data",nk(0:,:))

       if(plot3D)then
          call plot_keldysh_contour_gf(locG,t(0:),"PLOT/locG")
          if(volterra)then
             call plot_keldysh_contour_gf(locG1,t(0:),"PLOT/locG1")
             call plot_keldysh_contour_gf(locG2,t(0:),"PLOT/locG2")
          endif
       end if

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
          call splot(reg_filename(data_dir)//"/locChi_11.data",chi(1,1,0:,0:))
          call splot(reg_filename(data_dir)//"/locChi_12.data",chi(1,2,0:,0:))
          call splot(reg_filename(data_dir)//"/locChi_21.data",chi(2,1,0:,0:))
          call splot(reg_filename(data_dir)//"/locChi_22.data",chi(2,2,0:,0:))
          if(plot3D)then
             call splot("PLOT/locChi_11",t(0:),t(0:),chi(1,1,0:,0:))
             call splot("PLOT/locChi_12",t(0:),t(0:),chi(1,2,0:,0:))
             call splot("PLOT/locChi_21",t(0:),t(0:),chi(2,1,0:,0:))
             call splot("PLOT/locChi_22",t(0:),t(0:),chi(2,2,0:,0:))
          endif
       endif

    endif
  end subroutine print_out_Gloc



  !******************************************************************
  !******************************************************************
  !******************************************************************

end module KADANOFBAYM
