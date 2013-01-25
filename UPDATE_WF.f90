!###############################################################
!     PURPOSE  : Constructs some functions used in other places. 
!     AUTHORS  : Adriano Amaricci
!###############################################################
module UPDATE_WF
  USE MATRIX
  USE VARS_GLOBAL
  implicit none
  private

  public  :: neq_update_weiss_field

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Update the Weiss Fields G0^{>,<,Ret} using Dyson equation
  !+-------------------------------------------------------------------+
  subroutine neq_update_weiss_field
    integer :: M,i,j,k,itau,jtau,NN
    real(8) :: R,deg
    real(8) :: w,A,An
    complex(8),dimension(0:nstep,0:nstep) :: locGret,Sret
    complex(8),dimension(0:nstep,0:nstep) :: Uno,GammaRet
    complex(8),dimension(0:nstep,0:nstep) :: G0ret,G0adv
    !
    complex(8),dimension(:,:),allocatable :: mat_Delta
    complex(8),dimension(:,:),allocatable :: mat_Gamma
    complex(8),dimension(:,:),allocatable :: mat_G0,mat_Sigma,mat_locG
    complex(8),dimension(:,:),allocatable :: mat_SigmaHF
    !
    type(keldysh_contour_gf),save             :: G0_old

    if(G0_old%status.EQV..false.)call allocate_keldysh_contour_gf(G0_old,Nstep)
    G0_old=G0    

    call msg("Update WF")
    include "update_G0_nonequilibrium.f90"
    G0%less = weight*G0%less + (1.d0-weight)*G0_old%less
    G0%gtr  = weight*G0%gtr  + (1.d0-weight)*G0_old%gtr

    !Save data:
    if(mpiID==0)then
       call write_keldysh_contour_gf(G0,trim(data_dir)//"/G0")
       if(plot3D)call plot_keldysh_contour_gf(G0,t(0:),trim(plot_dir)//"/G0")
    end if
  end subroutine neq_update_weiss_field



  !********************************************************************
  !********************************************************************
  !********************************************************************



  function build_keldysh_matrix_gf(G,N) result(matG)
    type(keldysh_contour_gf)              :: G
    complex(8),dimension(0:2*N+1,0:2*N+1) :: matG
    integer                               :: i,j,N
    forall(i=0:N,j=0:N)
       matG(i,j)         = step(t(i)-t(j))*G%gtr(i,j) + step(t(j)-t(i))*G%less(i,j)
       matG(i,N+1+j)     =-G%less(i,j)
       matG(N+1+i,j)     = G%gtr(i,j)
       matG(N+1+i,N+1+j) =-(step(t(i)-t(j))*G%less(i,j)+ step(t(j)-t(i))*G%gtr(i,j))
    end forall
  end function build_keldysh_matrix_gf



  !********************************************************************
  !********************************************************************
  !********************************************************************





  ! subroutine update_equilibrium_weiss_field
  !   integer :: M,i,j,k,itau,jtau,NN
  !   real(8) :: R,deg
  !   real(8) :: w,A,An
  !   forall(i=0:nstep,j=0:nstep)
  !      gf%ret%t(i-j) = heaviside(t(i-j))*(locG%gtr(i,j)-locG%less(i,j))
  !      sf%ret%t(i-j) = heaviside(t(i-j))*(Sigma%gtr(i,j)-Sigma%less(i,j))
  !   end forall
  !   if(heaviside(0.d0)==1.d0)gf%ret%t(0)=gf%ret%t(0)/2.d0
  !   if(heaviside(0.d0)==1.d0)sf%ret%t(0)=sf%ret%t(0)/2.d0

  !   call fftgf_rt2rw(gf%ret%t,gf%ret%w,nstep) ; gf%ret%w=gf%ret%w*dt ; call swap_fftrt2rw(gf%ret%w)
  !   call fftgf_rt2rw(sf%ret%t,sf%ret%w,nstep) ; sf%ret%w=sf%ret%w*dt ; call swap_fftrt2rw(sf%ret%w)
  !   gf0%ret%w  = one/(one/gf%ret%w + sf%ret%w)
  !   gf0%less%w = less_component_w(gf0%ret%w,wr,beta)
  !   gf0%gtr%w  = gtr_component_w(gf0%ret%w,wr,beta)

  !   call fftgf_rw2rt(gf0%less%w,gf0%less%t,nstep) ; gf0%less%t=exa*fmesh/pi2*gf0%less%t
  !   call fftgf_rw2rt(gf0%gtr%w, gf0%gtr%t,nstep)  ; gf0%gtr%t =exa*fmesh/pi2*gf0%gtr%t
  !   call fftgf_rw2rt(gf0%ret%w, gf0%ret%t,nstep)  ; gf0%ret%t =exa*fmesh/pi2*gf0%ret%t
  !   forall(i=0:nstep,j=0:nstep)
  !      G0%less(i,j)= gf0%less%t(i-j)
  !      G0%gtr(i,j) = gf0%gtr%t(i-j)
  !   end forall
  !   ! !PLus this:
  !   ! forall(i=0:nstep,j=0:nstep)
  !   !    G0ret(i,j)=heaviside(t(i-j))*(G0gtr(i,j) - G0less(i,j))
  !   !    gf0%ret%t(i-j)=G0ret(i,j)
  !   ! end forall
  !   ! call fftgf_rt2rw(gf0%ret%t,gf0%less%w,nstep) ; gf0%less%w=gf0%less%w*dt ; call swap_fftrt2rw(gf0%less%w)
  ! end subroutine update_equilibrium_weiss_field




end module UPDATE_WF
