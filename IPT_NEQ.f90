!###############################################################
!     PURPOSE  : A non-equilibrium IPT solver module. 
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_NEQ
  USE VARS_GLOBAL
  USE EQUILIBRIUM
  USE MATRIX
  implicit none
  private

  public  :: neq_init_sigma
  public  :: neq_solve_ipt


contains


  subroutine neq_init_sigma()
    logical :: check,check1,check2
    inquire(file="Sless.data",exist=check1)
    if(.not.check1)inquire(file="Sless.data.gz",exist=check1)
    inquire(file="Sgtr.data",exist=check2)
    if(.not.check2)inquire(file="Sgtr.data.gz",exist=check2)
    check=check1.AND.check2
    if(check)then
       write(*,*)"Reading Sigma in input:"
       call sread("Sless.data",Sless(0:nstep,0:nstep))
       call sread("Sgtr.data",Sgtr(0:nstep,0:nstep))
    else
       print*,"Using Hartree-Fock self-energy"
       print*,"===================================="
       Sless=zero ; Sgtr=zero
    endif
  end subroutine neq_init_sigma


  !+-------------------------------------------------------------------+
  !PURPOSE  : BUild the 2^nd IPT sigma functions:
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt()
    integer                               :: i,j,itau
    !Get SIgma:
    call msg("Get Sigma(t,t')")

    forall(i=0:nstep,j=0:nstep)
       Sless(i,j) = (U**2)*(G0less(i,j)**2)*G0gtr(j,i)
       Sgtr (i,j) = (U**2)*(G0gtr(i,j)**2)*G0less(j,i)
    end forall

    !Get impurity GF and use SPT method if required
    call get_impuritygf
    if(method=="spt")then
       call msg("Recalculate Sigma using SPT")
       forall(i=0:nstep,j=0:nstep)
          Sless(i,j) = (U**2)*(impGless(i,j)**2)*impGgtr(j,i)
          Sgtr (i,j) = (U**2)*(impGgtr(i,j)**2)*impGless(j,i)                
       end forall
    end if

    !Save data:
    if(mpiID==0)then
       call splot("Sless.data",Sless(0:nstep,0:nstep))
       call splot("Sgtr.data",Sgtr(0:nstep,0:nstep))
    endif

  end subroutine Neq_solve_ipt



  !********************************************************************
  !********************************************************************
  !********************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : evaluate the impurity neq Green's functions
  !+-------------------------------------------------------------------+
  subroutine get_impuritygf()
    integer                               :: i,j
    real(8)                               :: A,w
    complex(8),dimension(0:nstep,0:nstep) :: Uno,GammaRet,Gamma0Ret
    complex(8),dimension(0:nstep,0:nstep) :: dG0ret,dGret,dSret

    if(update_wfftw)then
       call get_equilibrium_impuritygf !not tested!
    else
       include "obtain_Gimp_nonequilibrium.f90"
    endif

    !Save data:
    if(mpiID==0)then
       call splot("impGless.data",impGless(0:nstep,0:nstep))
       call splot("impGgtr.data",impGgtr(0:nstep,0:nstep))
    endif

  end subroutine get_impuritygf





  !********************************************************************
  !********************************************************************
  !********************************************************************


end module IPT_NEQ
