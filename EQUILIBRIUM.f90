!########################################################
!     Program  : HMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Hubbard model using DMFT-IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
module EQUILIBRIUM
  USE VARS_GLOBAL
  USE INTEGRATE
  implicit none
  private

  !Equilibrium Function
  complex(8),allocatable,dimension(:) :: sigma_,fg_,fg0_
  complex(8),allocatable,dimension(:) :: fiw_,siw_,f0iw_
  real(8),allocatable,dimension(:)    :: ftau_,stau_,f0tau_
  real(8),allocatable,dimension(:)    :: wr_,wm_

  public  :: solve_equilibrium_ipt
  !public :: get_equilibrium_impuritygf

contains

  subroutine solve_equilibrium_ipt
    integer    :: ik,i
    real(8)    :: n
    complex(8) :: zeta
    logical    :: init,check(2)
    if(mpiID==0)then
       call create_data_dir("Equilibrium")
       call solve_equilibrium_ipt_matsubara
       call solve_equilibrium_ipt_realaxis

       deallocate(sigma_,fg_,fg0_)
       deallocate(siw_,fiw_,f0iw_)
       deallocate(stau_,ftau_,f0tau_)
       deallocate(wr_,wm_)
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
    call splot("Equilibrium/Matsubara/G0_tau.ipt",tau(0:),f0tau_(0:))
    call splot("Equilibrium/Matsubara/Sigma_tau.ipt",tau(0:),stau_(0:))
  contains
    function solve_ipt_matsubara(fg0_) result(sigma_)
      complex(8),dimension(L)  :: fg0_
      complex(8),dimension(L)  :: sigma_
      real(8),dimension(0:L)   :: fg0tau_,sigmatau_
      call fftgf_iw2tau(fg0_,fg0tau_,beta)
      forall(i=0:L)sigmatau_(i)=U**2*(fg0tau_(i))**2*fg0tau_(L-i)
      call fftgf_tau2iw(sigmatau_,sigma_,beta)
    end function solve_ipt_matsubara
  end subroutine solve_equilibrium_ipt_matsubara


  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine solve_equilibrium_ipt_realaxis()
    integer                             :: i,j,loop
    logical                             :: converged
    complex(8)                          :: zeta
    real(8)                             :: n,z,wmax_,mesh
    complex(8),allocatable,dimension(:) :: sold,sigt
    integer,allocatable,dimension(:,:)  :: iy_m_ix
    real(8),dimension(:),allocatable    :: A0m,A0p,P1,P2
    call msg("Solving  Equilibrium problem on Real-axis:")
    call create_data_dir("Equilibrium/Real")
    allocate(fg_(L))
    allocate(sigma_(L))
    allocate(fg0_(L))
    allocate(sold(L))
    allocate(wr_(L))
    call get_frequency_index
    wmax_= min(20.d0,wmax)
    wr_ = linspace(-wmax_,wmax_,L,mesh=mesh)
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
       sold  = sigma_
       call solve_ipt_sopt()
       sigma_= weight*sigma_ + (1.d0-weight)*sold
       converged=check_convergence(sigma_,eps_error,nsuccess,eqNloop)
       call splot("Equilibrium/Real/nVSiloop.ipt",loop,n,append=TT)
    enddo
    call close_file("Equilibrium/Real/nVSiloop.ipt")
    call splot("Equilibrium/Real/DOS.ipt",wr_,-aimag(fg_)/pi)
    call splot("Equilibrium/Real/G_realw.ipt",wr_,fg_)
    call splot("Equilibrium/Real/G0_realw.ipt",wr_,fg0_)
    call splot("Equilibrium/Real/Sigma_realw.ipt",wr_,sigma_)
  contains
    subroutine solve_ipt_sopt()
      call getAs
      call getPolarization
      call Sopt
    end subroutine solve_ipt_sopt
    !
    subroutine get_frequency_index()
      integer :: ix,iy,iz
      if(.not.allocated(iy_m_ix))allocate(iy_m_ix(L,L))
      iy_m_ix=0
      do ix=1,L
         do iy=1,L
            iz = iy - ix + L/2 
            if(iz<1 .OR. iz>L) iz=-1 !out of range-> if(iz>-L)
            iy_m_ix(iy,ix)=iz
         enddo
      enddo
      if(.not.allocated(A0m))allocate(A0m(L))
      if(.not.allocated(A0p))allocate(A0p(L))
      if(.not.allocated(P1)) allocate(P1(L))
      if(.not.allocated(P2)) allocate(P2(L))
    end subroutine get_frequency_index
    !
    subroutine getAs
      real(8) :: dos(L)
      dos(:) =-aimag(fg0_(:))/pi
      A0p(:) = dos(:)*fermi(wr_(:),beta)
      A0m(:) = dos(:)*(1.d0-fermi(wr_(:),beta))
    end subroutine getAs
    !
    subroutine getPolarization
      integer :: ix,iy,iz    
      P1=zero
      P2=zero
      do ix=1,L
         do iy=1,L
            iz= iy_m_ix(iy,ix)
            if(iz>0)then
               P1(ix)=P1(ix) + A0p(iy)*A0m(iz)*mesh
               P2(ix)=P2(ix) + A0m(iy)*A0p(iz)*mesh
            endif
         enddo
      enddo
    end subroutine getPolarization
    !
    subroutine Sopt
      integer              :: ix,iy,iz
      real(8)              :: sum1,sum2
      real(8),dimension(L) :: reS,imS
      do ix=1,L
         sum1=zero
         sum2=zero
         do iy=1,L
            iz= iy_m_ix(iy,ix)
            if(iz>0)then
               sum1=sum1+A0p(L-iz+1)*P1(iy)*mesh
               sum2=sum2+A0m(L-iz+1)*P2(iy)*mesh
            end if
         enddo
         imS(ix)=-(U**2)*(sum1+sum2)*pi
      enddo
      reS = kronig(imS,wr_,size(ImS))
      sigma_ = reS + xi*imS
    end subroutine Sopt
  end subroutine solve_equilibrium_ipt_realaxis



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
  !      w = wr_(i)
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
