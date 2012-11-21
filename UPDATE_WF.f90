!###############################################################
!     PURPOSE  : Constructs some functions used in other places. 
!     AUTHORS  : Adriano Amaricci
!###############################################################
module UPDATE_WF
  USE VARS_GLOBAL
  USE ELECTRIC_FIELD
  USE BATH
  !USE EQUILIBRIUM
  USE MATRIX
  implicit none
  private

  public                           :: neq_update_weiss_field
  !public                           :: convergence_check

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Update the Weiss Fields G0^{>,<,Ret} using Dyson equation
  !+-------------------------------------------------------------------+
  subroutine neq_update_weiss_field
    integer                               :: i,j,k,itau,jtau
    real(8),dimension(0:nstep)            :: nt             !occupation(time)
    real(8),dimension(-Ltau:Ltau)         :: tmpG0tau
    complex(8),dimension(:,:),allocatable :: mat_calG,mat_Sigma,mat_locG
    type(kbm_contour_gf),save             :: G0_old
    integer                               :: t1min,t1max,t2min,t2max,t3min,t3max
    complex(8),allocatable,dimension(:)   :: dtloc


    ! if(G0_old%status.EQV..false.)call allocate_kbm_contour_gf(G0_old,Nstep,Ltau)
    ! G0_old=G0    


    call msg("Update WF: Dyson (no mixing)")
    if(update_wfftw)then
       call update_equilibrium_weiss_field
    else
       !==== Contour setup =====================================================
       !Contour indices
       t1min = 0         ; t1max=nstep          !size(nstep+1)
       t2min = nstep+1   ; t2max=2*nstep+1      !size(nstep+1)
       t3min = 2*nstep+2 ; t3max=2*nstep+Ltau+2 !size(Ltau+1)
       !Contour differential
       allocate(dtloc(t1min:t3max))
       dtloc(t1min:t1max)=  dt
       dtloc(t2min:t2max)= -dt
       dtloc(t3min:t3max)= -xi*dtau

       !====Update of the Matubara component first (stored in eq_G0tau)=========
       eq_G0iw = one/(one/eq_Giw + eq_Siw)
       call fftgf_iw2tau(eq_G0iw,eq_G0tau(0:),beta)
       forall(i=1:Ltau)eq_G0tau(-i)=-eq_G0tau(Ltau-i)

       !====direct inversion of the KBM Contour gf (in G0 *kbm_contour_gf)======
       allocate(mat_locG(t1min:t3max,t1min:t3max))
       mat_locG(0:,0:) = build_kbm_matrix_gf(locG,Nstep,Ltau)

       allocate(mat_Sigma(t1min:t3max,t1min:t3max))
       mat_Sigma(0:,0:) = build_kbm_matrix_gf(Sigma,Nstep,Ltau)

       forall(i=t1min:t3max,j=t1min:t3max)mat_locG(i,j)=dtloc(i)*mat_locG(i,j)*dtloc(j)
       call matrix_inverse(mat_locG(0:,0:))

       allocate(mat_calG(t1min:t3max,t1min:t3max))
       mat_calG(0:,0:) = mat_locG(0:,0:) + mat_Sigma(0:,0:)

       forall(i=t1min:t3max,j=t1min:t3max)mat_calG(i,j)=dtloc(i)*mat_calG(i,j)*dtloc(j)
       call matrix_inverse(mat_calG(0:,0:))

       !====scatter of the compontents into G0 *kbm_contour_gf)=================
       call scatter_kbm_matrix_gf(mat_calG(0:,0:),Nstep,Ltau,G0)

       forall(j=0:Ltau)G0%gmix(j,:)=conjg(G0%lmix(:,Ltau-j))  !enforce symmetry

    endif

    ! G0%less(0:,0:) = weight*G0%less(0:,0:) + (1.d0-weight)*G0_old%less(0:,0:)
    ! G0%gtr(0:,0:)  = weight*G0%gtr(0:,0:)  + (1.d0-weight)*G0_old%gtr(0:,0:)
    ! G0%lmix(0:,0:) = weight*G0%lmix(0:,0:) + (1.d0-weight)*G0_old%lmix(0:,0:)
    ! G0%gmix(0:,0:) = weight*G0%gmix(0:,0:) + (1.d0-weight)*G0_old%gmix(0:,0:)
    ! G0%mats(0:,0:) = weight*G0%mats(0:,0:) + (1.d0-weight)*G0_old%mats(0:,0:)

    !Save data:
    if(mpiID==0)then
       call write_kbm_contour_gf(G0,reg_filename(data_dir)//"/G0")
       if(plot3D)call plot_kbm_contour_gf(G0,t(0:),tau(0:),trim(plot_dir)//"/G0")
       call splot("G0_less_t0.ipt",t(0:),G0%less(0:,0))
       call splot("G0_lmix_tau0.ipt",t(0:),G0%lmix(0:,0))
       call splot("eq_G0_iw.ipt",wm,eq_G0iw,append=.true.)
       forall(i=0:nstep)nt(i)=-xi*G0%less(i,i)
       call splot("n0VStime.ipt",t(0:nstep),2.d0*nt(0:nstep),append=TT)
       forall(i=0:Ltau,j=0:Ltau)tmpG0tau(i-j)=G0%mats(i,j)
       call splot("eq_G0_tau.ipt",tau,eq_G0tau,tmpG0tau,append=.true.)
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

  function build_kbm_matrix_gf(G,N,L) result(matG)
    type(kbm_contour_gf)                      :: G
    integer                                   :: i,j,N,L
    complex(8),dimension(0:2*N+L+2,0:2*N+L+2) :: matG
    !
    real(8)                                   :: ftau(-L:L),fdtau
    complex(8)                                :: Glmix(0:N,0:L)
    real(8)                                   :: Gmats(0:L,0:L),eqGtau(-L:L),upmGtau(0:L)
    !insert a layer to interpolate the UPM functions to LinM in tau:
    ftau = linspace(-beta,beta,2*L+1,mesh=fdtau)               !linear mesh in tau [-beta,beta]
    forall(i=0:L)upmGtau(i) = G%mats(i,0)                      !get the eq. G(tau) on the UPM grid:
    call cubic_spline(upmGtau(0:),tau(0:),eqGtau(0:),ftau(0:)) !interpolate G^M(tau) in (0:beta)
    forall(i=1:L)eqGtau(-i)=-eqGtau(L-i)                       !get the (-beta:0) part
    forall(i=0:L,j=0:L)Gmats(i,j)=eqGtau(i-j)                  !build G^M(tau,tau`)
    do i=0,N                                                   !interpolate G^lmix(t,tau`)
       call cubic_spline(G%lmix(i,0:),tau(0:),Glmix(i,0:),ftau(0:))
    enddo
    !
    matG=zero
    forall(i=0:N,j=0:N)
       matG(i,        j)   = step(t(i)-t(j))*G%gtr(i,j)  + step(t(j)-t(i))*G%less(i,j)
       matG(i,    N+1+j)   = G%less(i,j)
       matG(N+1+i,    j)   = G%gtr(i,j)
       matG(N+1+i,N+1+j)   = (step(t(i)-t(j))*G%less(i,j)+ step(t(j)-t(i))*G%gtr(i,j))
    end forall

    forall(i=0:N,j=0:L)
       matG(i    ,2*N+2+j) =  Glmix(i,j)!G%lmix(i,j)
       matG(N+1+i,2*N+2+j) =  Glmix(i,j)!G%lmix(i,j)
    end forall

    forall(i=0:L,j=0:N)
       matG(2*N+2+i,    j) =  conjg(Glmix(j,L-i))!G%lmix(j,L-i))  !=G%gmix
       matG(2*N+2+i,N+1+j) =  conjg(Glmix(j,L-i))!G%lmix(j,L-i))  !=G%gmix
    end forall

    forall(i=0:L,j=0:L)matG(2*N+2+i,2*N+2+j) = xi*Gmats(i,j)!G%mats(i,j)

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
    forall(i=0:L)G%gmix(i,:)=conjg(G%lmix(:,Ltau-i))
    G%mats(0:L,0:L) = aimag(matG(2*N+2:2*N+2+L,2*N+2:2*N+2+L))
  end subroutine scatter_kbm_matrix_gf

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


end module UPDATE_WF
