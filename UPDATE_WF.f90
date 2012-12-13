!###############################################################
!     PURPOSE  : Constructs some functions used in other places. 
!     AUTHORS  : Adriano Amaricci
!###############################################################
module UPDATE_WF
  USE VARS_GLOBAL
  USE MATRIX
  implicit none
  private

  public  :: neq_update_weiss_field

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Update the Weiss Fields G0^{>,<,Ret} using Dyson equation
  !+-------------------------------------------------------------------+
  subroutine neq_update_weiss_field
    integer                               :: i,j,k,itau,jtau
    real(8),dimension(0:nstep)            :: nt             !occupation(time)
    complex(8),dimension(:,:),allocatable :: mat_calG,mat_Sigma,mat_locG
    type(kbm_contour_gf),save             :: G0_old
    integer                               :: t1min,t1max,t2min,t2max,t3min,t3max
    complex(8),allocatable,dimension(:)   :: dtloc


    if(G0_old%status.EQV..false.)call allocate_kbm_contour_gf(G0_old,Nstep,Ltau)
    G0_old=G0    


    call msg("Update WF: Dyson (no mixing)")
    if(update_wfftw)then
       call update_equilibrium_weiss_field
    else
       !==== Contour setup =====================================================
       !Contour indices
       t1min = 0         ; t1max=nstep          !size(nstep+1)
       t2min = nstep+1   ; t2max=nstep+1+nstep  !size(nstep+1)
       t3min = 2*nstep+2 ; t3max=2*nstep+2+Ltau !size(Ltau+1)
       !Contour differential
       allocate(dtloc(t1min:t3max))
       dtloc(t1min:t1max)=  dt
       dtloc(t2min:t2max)= -dt
       if(upmflag)then
          dtloc(t3min)=tau(1)-tau(0)
          do i=1,Ltau
             dtloc(t3min+i)=tau(i)-tau(i-1)
          enddo
          dtloc(t3min:t3max)= -xi*dtloc(t3min:t3max)
       else
          dtloc(t3min:t3max)= -xi*dtau
       endif



       !====direct inversion of the KBM Contour gf (in G0 *kbm_contour_gf)======
       allocate(mat_locG(t1min:t3max,t1min:t3max))
       allocate(mat_Sigma(t1min:t3max,t1min:t3max))
       allocate(mat_calG(t1min:t3max,t1min:t3max))

       mat_locG(0:,0:) = build_kbm_matrix_gf(locG,Nstep,Ltau)
       mat_Sigma(0:,0:)= build_kbm_matrix_gf(Sigma,Nstep,Ltau)

       !====inversion of the KBM Contour local GF ==============================
       forall(i=t1min:t3max,j=t1min:t3max)mat_locG(i,j)=dtloc(i)*mat_locG(i,j)*dtloc(j)
       call matrix_inverse(mat_locG(0:,0:))

       !====build calG_0^-1 = G_loc^-1 + Sigma =================================
       mat_calG(0:,0:) = mat_locG(0:,0:) + mat_Sigma(0:,0:)

       !====inversion of the KBM Contour calG_0 ================================
       forall(i=t1min:t3max,j=t1min:t3max)mat_calG(i,j)=dtloc(i)*mat_calG(i,j)*dtloc(j)
       call matrix_inverse(mat_calG(0:,0:))

       !====scatter of the compontents into G0 *kbm_contour_gf)=================
       call scatter_kbm_matrix_gf(mat_calG(0:,0:),Nstep,Ltau,G0)

       !====Update of the Matubara component=========
       eq_G0iw = one/(one/eq_Giw + eq_Siw)
       ! if(upmflag)then
       !    call fftgf_iw2tau_upm(wm,eq_G0iw,tau(0:),G0%mats(0:),beta)
       ! else
       !    call fftgf_iw2tau(eq_G0iw,G0%mats(0:),beta)
       ! endif
       ! forall(i=1:Ltau)G0%mats(-i)=-G0%mats(Ltau-i)
    endif

    !Save data:
    if(mpiID==0)then
       call write_kbm_contour_gf(G0,reg_filename(data_dir)//"/G0")
       if(plot3D)call plot_kbm_contour_gf(G0,t(0:),tau(0:),trim(plot_dir)//"/G0")
       call splot("eq_G0_tau.ipt",tau,G0%mats,append=.true.)
       call splot("eq_G0_iw.ipt",wm,eq_G0iw,append=.true.)
       call splot("G0_less_t0.ipt",t(0:),G0%less(0:,0))
       call splot("G0_lmix_tau0.ipt",t(0:),G0%lmix(0:,0))
       call splot("G0_lmix_t0_tau.ipt",tau(0:),G0%lmix(0,0:))
       forall(i=0:nstep)nt(i)=-xi*G0%less(i,i)
       call splot("n0VStime.ipt",t(0:nstep),2.d0*nt(0:nstep),append=TT)
    end if


    G0%less(0:,0:) = weight*G0%less(0:,0:) + (1.d0-weight)*G0_old%less(0:,0:)
    G0%gtr(0:,0:)  = weight*G0%gtr(0:,0:)  + (1.d0-weight)*G0_old%gtr(0:,0:)
    G0%lmix(0:,0:) = weight*G0%lmix(0:,0:) + (1.d0-weight)*G0_old%lmix(0:,0:)
    G0%gmix(0:,0:) = weight*G0%gmix(0:,0:) + (1.d0-weight)*G0_old%gmix(0:,0:)
    G0%mats(:)     = weight*G0%mats(:)     + (1.d0-weight)*G0_old%mats(:)

  end subroutine neq_update_weiss_field

  !********************************************************************
  !********************************************************************
  !********************************************************************

  function build_kbm_matrix_gf(G,N,L) result(matG)
    type(kbm_contour_gf)                      :: G
    integer                                   :: i,j,N,L
    complex(8),dimension(0:2*N+L+2,0:2*N+L+2) :: matG
    matG=zero
    forall(i=0:N,j=0:N)
       matG(i,        j)   = step(t(i)-t(j))*G%gtr(i,j)  + step(t(j)-t(i))*G%less(i,j)
       matG(i,    N+1+j)   = G%less(i,j)
       matG(N+1+i,    j)   = G%gtr(i,j)
       matG(N+1+i,N+1+j)   = (step(t(i)-t(j))*G%less(i,j)+ step(t(j)-t(i))*G%gtr(i,j))
    end forall
    forall(i=0:N,j=0:L)
       matG(i    ,2*N+2+j) =  G%lmix(i,j)
       matG(N+1+i,2*N+2+j) =  G%lmix(i,j)
    end forall
    forall(i=0:L,j=0:N)
       matG(2*N+2+i,    j) =  conjg(G%lmix(j,L-i))
       matG(2*N+2+i,N+1+j) =  conjg(G%lmix(j,L-i))
    end forall
    forall(i=0:L,j=0:L)matG(2*N+2+i,2*N+2+j) = xi*G%mats(i-j)
  end function build_kbm_matrix_gf

  !********************************************************************
  !********************************************************************
  !********************************************************************

  subroutine scatter_kbm_matrix_gf(matG,N,L,G)
    integer              :: i,j,N,L
    complex(8)           :: matG(0:2*N+L+2,0:2*N+L+2)
    type(kbm_contour_gf) :: G
    if(.not.G%status)call error("Error 1")
    if(G%N/=N)call error("Error 2: N")
    if(G%L/=L)call error("Error 3: L")
    forall(i=0:N,j=0:N)
       G%less(i,j) =  matG(i,N+1+j) !matG12
       G%gtr(i,j)  =  matG(N+1+i,j) !matG21
    end forall
    G%lmix(0:N,0:L) =  matG(0:N,2*N+2:2*N+2+L) !matG13/matG23
    forall(i=0:L)G%gmix(i,:)=conjg(G%lmix(:,L-i))    
    G%mats(0:L) = dimag(matG(2*N+2:2*N+2+L,2*N+2))!,2*N+2:2*N+2+L))
    forall(i=1:L)G%mats(-i)=-G%mats(L-i)
  end subroutine scatter_kbm_matrix_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

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


end module UPDATE_WF
