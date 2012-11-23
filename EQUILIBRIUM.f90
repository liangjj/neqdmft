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

  public  :: solve_equilibrium_ipt

contains

  subroutine solve_equilibrium_ipt
    if(mpiID==0)then
       call create_data_dir("Equilibrium")
       call solve_equilibrium_ipt_matsubara
       call solve_equilibrium_ipt_realaxis
    endif
  end subroutine solve_equilibrium_ipt



  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine solve_equilibrium_ipt_matsubara()
    integer                 :: i,j,loop
    logical                 :: converged
    complex(8)              :: zeta
    real(8)                 :: n,z,wmax_
    complex(8),dimension(L) :: fiw_,siw_,f0iw_,sold
    real(8),dimension(0:L)  :: ftau_,stau_,f0tau_
    real(8),dimension(L)    :: wm_
    real(8),dimension(0:L)  :: tau_
    call msg("Solving  Equilibrium problem in Matsubara:")
    call create_data_dir("Equilibrium/Matsubara")
    wm_      = pi/beta*real(2*arange(1,L)-1,8)
    tau_(0:) = linspace(0.d0,beta,L+1)
    siw_     = zero
    loop     = 0 ; converged=.false.
    do while (.not.converged)
       loop=loop+1
       sold=f0iw_
       write(*,"(A,i5)",advance="no")"DMFT-loop",loop
       do i=1,L
          zeta    = xi*wm_(i) - siw_(i)
          fiw_(i) = sum_overk_zeta(zeta,epsik,wt)
       enddo
       call  fftgf_iw2tau(fiw_,ftau_(0:),beta)
       n     =-real(ftau_(L))
       f0iw_= one/(one/fiw_ + siw_)
       !
       call fftgf_iw2tau(f0iw_,f0tau_,beta)
       forall(i=0:L)stau_(i)=U**2*(f0tau_(i))**2*f0tau_(L-i) 
       call fftgf_tau2iw(stau_,siw_,beta)
       !
       f0iw_ = weight*f0iw_ + (1.d0-weight)*sold
       converged=check_convergence(f0iw_,eps_error,Nsuccess,eqNloop)
       z=1.d0 - dimag(siw_(1))/wm_(1);z=1.d0/z
       call splot("Equilibrium/Matsubara/observables_all.ipt",&
            dfloat(loop),beta,u,n,z,append=TT)
    enddo
    call close_file("Equilibrium/Matsubara/observables_all.ipt")
    call splot("Equilibrium/Matsubara/G_iw.ipt",wm_,fiw_)
    call splot("Equilibrium/Matsubara/G0_iw.ipt",wm_,f0iw_)
    call splot("Equilibrium/Matsubara/Sigma_iw.ipt",wm_,siw_)
    call splot("Equilibrium/Matsubara/G0_tau.ipt",tau_(0:),f0tau_(0:))
    call splot("Equilibrium/Matsubara/Sigma_tau.ipt",tau_(0:),stau_(0:))
  end subroutine solve_equilibrium_ipt_matsubara


  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine solve_equilibrium_ipt_realaxis()
    integer                 :: i,j,loop
    logical                 :: converged
    complex(8)              :: zeta
    real(8)                 :: n,z,wmax_,mesh
    complex(8),dimension(L) :: sigma_,fg_,fg0_,sold
    real(8),dimension(L)    :: wr_
    integer,dimension(L,L)  :: iy_m_ix
    real(8),dimension(L)    :: A0m,A0p,P1,P2
    call msg("Solving  Equilibrium problem on Real-axis:")
    call create_data_dir("Equilibrium/Real")
    call get_frequency_index
    wmax_ = min(20.d0,wmax)
    wr_   = linspace(-wmax_,wmax_,L,mesh=mesh)
    sigma_= zero ; sold=sigma_
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
       call getAs
       call getPolarization
       call sopt
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
    subroutine get_frequency_index()
      integer :: ix,iy,iz
      iy_m_ix=0
      do ix=1,L
         do iy=1,L
            iz = iy - ix + L/2 
            if(iz<1 .OR. iz>L) iz=-1 !out of range-> if(iz>-L)
            iy_m_ix(iy,ix)=iz
         enddo
      enddo
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
      reS = kronig(imS,wr_,L)
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
