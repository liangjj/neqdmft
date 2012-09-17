!###############################################################
!     PURPOSE  : Constructs some functions used in other places. 
!     AUTHORS  : Adriano Amaricci
!###############################################################
module FUNX_NEQ
  USE MATRIX
  USE VARS_GLOBAL
  USE ELECTRIC_FIELD
  USE BATH
  USE EQUILIBRIUM
  implicit none
  private

  public                           :: neq_update_weiss_field
  public                           :: print_observables
  public                           :: convergence_check
contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Update the Weiss Fields G0^{>,<,Ret} using Dyson equation
  !+-------------------------------------------------------------------+
  subroutine neq_update_weiss_field
    integer :: M,i,j,k,itau,jtau,NN
    real(8) :: R,deg
    real(8) :: w,A,An
    !
    complex(8),dimension(:,:),allocatable :: mat_calG,mat_Sigma,mat_locG
    !
    type(kbm_contour_gf),save             :: G0_old

    if(G0_old%status.EQV..false.)call allocate_kbm_contour_gf(G0_old,Nstep,Ltau)
    G0_old=G0    

    call msg("Update WF: Dyson")
    if(update_wfftw)then
       call update_equilibrium_weiss_field
    else
       include "update_G0_nonequilibrium.f90"
    endif

    G0%less(0:,0:) = weight*G0%less(0:,0:) + (1.d0-weight)*G0_old%less(0:,0:)
    G0%gtr(0:,0:)  = weight*G0%gtr(0:,0:)  + (1.d0-weight)*G0_old%gtr(0:,0:)
    G0%lmix(0:,0:) = weight*G0%lmix(0:,0:) + (1.d0-weight)*G0_old%lmix(0:,0:)
    G0%gmix(0:,0:) = weight*G0%gmix(0:,0:) + (1.d0-weight)*G0_old%gmix(0:,0:)
    G0%mats(0:,0:) = weight*G0%mats(0:,0:) + (1.d0-weight)*G0_old%mats(0:,0:)

    !Save data:
    if(mpiID==0)then
       call write_kbm_contour_gf(G0,reg_filename(data_dir)//"/G0")
       if(plot3D)call plot_kbm_contour_gf(G0,t(0:),tau(0:),"PLOT/G0")
       call splot("PLOT/G0_less_t0.ipt",t(0:),G0%less(0:,0))
       call splot("PLOT/G0_lmix_tau0.ipt",t(0:),G0%lmix(0:,0))
    end if
  end subroutine neq_update_weiss_field





  !********************************************************************
  !********************************************************************
  !********************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : evaluate and print some observables to check the calculation
  !+-------------------------------------------------------------------+
  subroutine print_observables
    integer                          :: i,ik,ix,iy
    type(vect2D)                     :: Jk,Ak
    type(vect2D),dimension(0:nstep)  :: Jloc                   !local Current 
    real(8),dimension(0:nstep)       :: nt,modJloc             !occupation(time)
    if(mpiID==0)then
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
  end subroutine print_observables



  !********************************************************************
  !********************************************************************
  !********************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : check convergence of the calculation:
  !+-------------------------------------------------------------------+
  function convergence_check() result(converged)
    logical                         :: converged
    integer                         :: i,ik,ix,iy
    type(vect2D)                    :: Jk,Ak
    type(vect2D),dimension(0:nstep) :: Jloc                   !local Current 
    real(8),dimension(0:nstep)      :: test_func
    integer                         :: selector

    if(mpiID==0)then
       if(solve_wfftw)then
          forall(i=0:nstep)test_func(i)=-xi*locG%less(i,i)
       elseif(Efield/=0.d0)then
          Jloc=Vzero
          do ik=1,Lk
             ix=ik2ix(ik);iy=ik2iy(ik)
             do i=0,nstep
                Ak= Afield(t(i),Ek)
                Jk= nk(i,ik)*square_lattice_velocity(kgrid(ix,iy) - Ak)
                Jloc(i) = Jloc(i) +  wt(ik)*Jk
             enddo
          enddo
          test_func(0:nstep)=modulo(Jloc(0:nstep))
       else
          forall(i=0:nstep)test_func(i)=-xi*locG%less(i,i)
       endif
       converged=check_convergence(test_func(0:nstep),eps_error,Nsuccess,nloop,id=0)
    endif
  end function convergence_check





end module FUNX_NEQ
