!########################################################
!     Program  : HMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Hubbard model using DMFT-IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
module EQUILIBRIUM
  USE VARS_GLOBAL
  USE IPT_SOPT
  USE IPT_MATS
  implicit none
  private

  !Equilibrium Function
  complex(8),allocatable,dimension(:) :: sigma_,fg_,fg0_
  complex(8),allocatable,dimension(:) :: fiw_,siw_,f0iw_
  real(8),allocatable,dimension(:)    :: ftau_,stau_,f0tau_
  real(8),allocatable,dimension(:)    :: wr_,wm_

  public  :: solve_equilibrium_ipt
  public  :: update_equilibrium_weiss_field
  public  :: get_equilibrium_localgf
  !public :: get_equilibrium_impuritygf

contains

  subroutine solve_equilibrium_ipt
    integer    :: ik,i
    real(8)    :: n
    complex(8) :: zeta
    logical    :: init,check(2)
    if(mpiID==0)then
       ! init=.true.
       ! do i=1,3
       !    inquire(file=trim(irdFILE(i)),exist=check(i))
       !    if(.not.check(i))inquire(file=trim(irdFILE(i))//".gz",exist=check(i))
       !    init=init.AND.check(i)
       ! enddo
       ! if(init)then                !equilibrium solution exists:
       !    call msg("Equilibrium solution exist:")
       ! else

       call create_data_dir("Equilibrium")
       call solve_equilibrium_ipt_matsubara
       call solve_equilibrium_ipt_realaxis

       ! !Write the Weiss-fields G0(w):
       ! call splot(trim(irdFILE(1)),wr_,fg0_)

       ! !Write the Sigma(iw)
       ! call splot(trim(irdFILE(2)),stau_(0:))

       ! !Write the Weiss-fields G0(iw):
       ! call splot(trim(irdFILE(3)),f0tau_(0:))

       deallocate(sigma_,fg_,fg0_)
       deallocate(siw_,fiw_,f0iw_)
       deallocate(stau_,ftau_,f0tau_)
       deallocate(wr_,wm_)
       ! endif
    endif
  end subroutine solve_equilibrium_ipt



  !******************************************************************
  !******************************************************************
  !******************************************************************



  subroutine solve_equilibrium_ipt_matsubara()
    integer    :: i,j,loop
    logical    :: converged
    complex(8) :: zeta
    real(8)    :: n,z,wmax_
    complex(8),allocatable,dimension(:) :: sold
    call msg("Solving  Equilibrium problem in Matsubara:")
    call create_data_dir("Equilibrium/Matsubara")
    allocate(fiw_(L),siw_(L),f0iw_(L))
    allocate(ftau_(0:Ltau),stau_(0:Ltau),f0tau_(0:Ltau))
    allocate(sold(L))
    allocate(wm_(L))
    wm_ = pi/beta*real(2*arange(1,L)-1,8)
    siw_=zero ; sold=siw_
    loop=0 ; converged=.false.
    do while (.not.converged)
       loop=loop+1
       write(*,"(A,i5)",advance="no")"DMFT-loop",loop
       do i=1,L
          zeta    = xi*wm_(i) - siw_(i)
          fiw_(i) = sum_overk_zeta(zeta,epsik,wt)
       enddo
       call  fftgf_iw2tau(fiw_,ftau_(0:),beta)
       n     =-real(ftau_(Ltau))
       f0iw_= one/(one/fiw_ + siw_)
       siw_ = solve_ipt_matsubara(f0iw_)
       siw_ = weight*siw_ + (1.d0-weight)*sold ; sold=siw_
       converged=check_convergence(siw_,eps_error,Nsuccess,eqNloop)
       z=1.d0 - dimag(siw_(1))/wm_(1);z=1.d0/z
       call splot("Equilibrium/Matsubara/observables_all.ipt",dble(loop),beta,u,n,z,append=TT)
    enddo
    !Use IPT to get S(tau).
    call fftgf_iw2tau(f0iw_,f0tau_(0:),beta)
    forall(i=0:Ltau)stau_(i)=U**2*(f0tau_(i))**2*f0tau_(Ltau-i) 
    call close_file("Equilibrium/Matsubara/observables_all.ipt")
    call splot("Equilibrium/Matsubara/G_iw.ipt",wm_,fiw_)
    call splot("Equilibrium/Matsubara/G0_iw.ipt",wm_,f0iw_)
    call splot("Equilibrium/Matsubara/Sigma_iw.ipt",wm_,siw_)
    call splot("Equilibrium/Matsubara/G0_tau.ipt",taureal(0:),f0tau_(0:))
    call splot("Equilibrium/Matsubara/Sigma_tau.ipt",taureal(0:),stau_(0:))
  end subroutine solve_equilibrium_ipt_matsubara


  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine solve_equilibrium_ipt_realaxis()
    integer    :: i,j,loop
    logical    :: converged
    complex(8) :: zeta
    real(8)    :: n,z,wmax_
    complex(8),allocatable,dimension(:) :: sold
    call msg("Solving  Equilibrium problem on Real-axis:")
    call create_data_dir("Equilibrium/Real")
    allocate(fg_(L))
    allocate(sigma_(L))
    allocate(fg0_(L))
    allocate(sold(L))
    allocate(wr_(L))
    wmax_= min(20.d0,wmax)
    wr_ = linspace(-wmax_,wmax_,L)
    sigma_=zero ; sold=sigma_
    loop=0 ; converged=.false.             
    do while (.not.converged)
       loop=loop+1
       write(*,"(A,i5)",advance="no")"DMFT-loop",loop
       do i=1,L
          zeta  = cmplx(wr_(i),eps) - sigma_(i)
          fg_(i) = sum_overk_zeta(zeta,epsik,wt)
       enddo
       n      = sum(aimag(fg_)*fermi(wr_,beta))/sum(aimag(fg_))
       fg0_   = one/(one/fg_ + sigma_)
       sigma_= solve_ipt_sopt(fg0_,wr_)
       sigma_= weight*sigma_ + (1.d0-weight)*sold
       sold  = sigma_
       converged=check_convergence(sigma_,eps_error,nsuccess,eqNloop)
       call splot("Equilibrium/Real/nVSiloop.ipt",loop,n,append=TT)
    enddo
    call close_file("Equilibrium/Real/nVSiloop.ipt")
    call splot("Equilibrium/Real/DOS.ipt",wr_,-aimag(fg_)/pi)
    call splot("Equilibrium/Real/G_realw.ipt",wr_,fg_)
    call splot("Equilibrium/Real/G0_realw.ipt",wr_,fg0_)
    call splot("Equilibrium/Real/Sigma_realw.ipt",wr_,sigma_)
  end subroutine solve_equilibrium_ipt_realaxis



  !********************************************************************
  !********************************************************************
  !********************************************************************



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



  subroutine update_equilibrium_weiss_field
    integer :: M,i,j,k,itau,jtau,NN
    real(8) :: R,deg
    real(8) :: w,A,An
    forall(i=0:nstep,j=0:nstep)
       gf%ret%t(i-j) = heaviside(t(i-j))*(locG%gtr(i,j)-locG%less(i,j))
       sf%ret%t(i-j) = heaviside(t(i-j))*(Sigma%gtr(i,j)-Sigma%less(i,j))
    end forall
    if(heaviside(0.d0)==1.d0)gf%ret%t(0)=gf%ret%t(0)/2.d0
    if(heaviside(0.d0)==1.d0)sf%ret%t(0)=sf%ret%t(0)/2.d0

    call fftgf_rt2rw(gf%ret%t,gf%ret%w,nstep) ; gf%ret%w=gf%ret%w*dt ; call swap_fftrt2rw(gf%ret%w)
    call fftgf_rt2rw(sf%ret%t,sf%ret%w,nstep) ; sf%ret%w=sf%ret%w*dt ; call swap_fftrt2rw(sf%ret%w)
    gf0%ret%w  = one/(one/gf%ret%w + sf%ret%w)
    gf0%less%w = less_component_w(gf0%ret%w,wr,beta)
    gf0%gtr%w  = gtr_component_w(gf0%ret%w,wr,beta)
    ! call splot("updateG0ret_w.ipt",wr,gf0%ret%w,append=TT)
    ! call splot("updateG0less_w.ipt",wr,gf0%less%w,append=TT)
    ! call splot("updateG0gtr_w.ipt",wr,gf0%gtr%w,append=TT)

    call fftgf_rw2rt(gf0%less%w,gf0%less%t,nstep) ; gf0%less%t=exa*fmesh/pi2*gf0%less%t
    call fftgf_rw2rt(gf0%gtr%w, gf0%gtr%t,nstep)  ; gf0%gtr%t =exa*fmesh/pi2*gf0%gtr%t
    call fftgf_rw2rt(gf0%ret%w, gf0%ret%t,nstep)  ; gf0%ret%t =exa*fmesh/pi2*gf0%ret%t
    forall(i=0:nstep,j=0:nstep)
       G0%less(i,j)= gf0%less%t(i-j)
       G0%gtr(i,j) = gf0%gtr%t(i-j)
    end forall
    ! call splot("updateG0ret_t.ipt",t,gf0%ret%t,append=TT)
    ! call splot("G0less3D",t(0:nstep)/dt,t(0:nstep)/dt,G0less(0:nstep,0:nstep))
    ! call splot("G0gtr3D",t(0:nstep)/dt,t(0:nstep)/dt,G0gtr(0:nstep,0:nstep))

    ! !PLus this:
    ! forall(i=0:nstep,j=0:nstep)
    !    G0ret(i,j)=heaviside(t(i-j))*(G0gtr(i,j) - G0less(i,j))
    !    gf0%ret%t(i-j)=G0ret(i,j)
    ! end forall
    ! call fftgf_rt2rw(gf0%ret%t,gf0%less%w,nstep) ; gf0%less%w=gf0%less%w*dt ; call swap_fftrt2rw(gf0%less%w)
  end subroutine update_equilibrium_weiss_field




  !********************************************************************
  !********************************************************************
  !********************************************************************



  ! subroutine get_equilibrium_impuritygf
  !   integer :: i,j,itau
  !   real(8) :: A,w

  !   forall(i=0:nstep,j=0:nstep)
  !      gf0%ret%t(i-j)=heaviside(t(i-j))*(G0gtr(i,j) - G0less(i,j))
  !      sf%ret%t(i-j)=heaviside(t(i-j))*(Sgtr(i,j) - Sless(i,j))
  !   end forall
  !   if(heaviside(0.d0)==1.d0)gf0%ret%t(0)=gf0%ret%t(0)/2.d0
  !   if(heaviside(0.d0)==1.d0)sf%ret%t(0)=sf%ret%t(0)/2.d0
  !   call fftgf_rt2rw(gf0%ret%t,gf0%ret%w,nstep) ;  gf0%ret%w=gf0%ret%w*dt ; call swap_fftrt2rw(gf0%ret%w) !swap because F(t) are not oscillating in this formalism:
  !   call fftgf_rt2rw(sf%ret%t,sf%ret%w,nstep)   ;  sf%ret%w=dt*sf%ret%w   ; call swap_fftrt2rw(sf%ret%w)   !swap because F(t) are not oscillating in this formalism:
  !   gf%ret%w = one/(one/gf0%ret%w - sf%ret%w)
  !   do i=1,2*nstep
  !      w = wr(i)
  !      A=-aimag(gf%ret%w(i))/pi
  !      gf%less%w(i)= pi2*xi*fermi(w,beta)*A
  !      gf%gtr%w(i) = pi2*xi*(fermi(w,beta)-1.d0)*A
  !   enddo
  !   call fftgf_rw2rt(gf%less%w,gf%less%t,nstep)  ; gf%less%t=fmesh/pi2*gf%less%t ;  gf%less%t=gf%less%t*exa 
  !   call fftgf_rw2rt(gf%gtr%w,gf%gtr%t,nstep)   ; gf%gtr%t =fmesh/pi2*gf%gtr%t  ;  gf%gtr%t=gf%gtr%t*exa
  !   forall(i=0:nstep,j=0:nstep)
  !      impGless(i,j)= gf%less%t(i-j)
  !      impGgtr(i,j) = gf%gtr%t(i-j)
  !   end forall

  ! end subroutine get_equilibrium_impuritygf


end module EQUILIBRIUM
