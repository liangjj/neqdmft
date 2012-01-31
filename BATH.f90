!###############################################################
!     PROGRAM  : BATH
!     TYPE     : Module
!     PURPOSE  : Constructs some functions used in other places. 
!     AUTHORS  : Adriano Amaricci
!     LAST UPDATE: 07/2009
!###############################################################
module BATH
  USE VARS_GLOBAL
  implicit none
  private
  integer,parameter                :: Lw=1024 !# of frequencies
  real(8),allocatable,dimension(:) :: bath_dens,bath_epsik,wfreq
  public                           :: get_Bath

contains
  !+-------------------------------------------------------------------+
  !PURPOSE  : Build the Bath part of the system using exclusively time
  !dependent formalism. The bath is not interacting so it is 
  ! time translation invariant.
  !+-------------------------------------------------------------------+
  subroutine get_Bath()
    integer          :: iw,i
    real(8)          :: en,w,ebini,de,dw,wfin,wini
    complex(8)       :: peso
    real(8)          :: ngtr,nless,arg

    print '(A)',"Get Bath:"
    call system("if [ ! -d BATH ]; then mkdir BATH; fi")
    call msg("Using "//trim(adjustl(trim(bath_type)))//" dissipative bath")

    allocate(bath_epsik(Lmu),bath_dens(-Lw:Lw),wfreq(-Lw:Lw))

    !fake linear dispersion
    ebini=-Wbath/2.d0 ; de=Wbath/dble(Lmu)
    forall(i=1:Lmu)bath_epsik(i)=ebini + real(i,8)*de

    select case(bath_type)
    case default
       wfin=Wbath ; wini=-wfin
       wfreq(-Lw:Lw)= linspace(wini,wfin,2*Lw+1,mesh=dw)
       call get_bath_constant_dos()

    case("bethe")
       wfin=Wbath ; wini=-wfin
       wfreq(-Lw:Lw)= linspace(wini,wfin,2*Lw+1,mesh=dw)
       call get_bath_bethe_dos()

    case("gaussian")
       wfin=4.d0*Wbath ; wini=-wfin
       wfreq(-Lw:Lw)= linspace(wini,wfin,2*Lw+1,mesh=dw)
       call get_bath_gaussian_dos()
    end select


    S0less=zero ; S0gtr=zero
    do iw=-Lw,Lw
       en   = wfreq(iw)
       nless= fermi0(en,beta)
       ngtr = fermi0(en,beta)-1.d0
       do i=-nstep,nstep
          peso=exp(-xi*t(i)*en)
          S0less(i)=S0less(i)+ xi*Vpd**2*nless*peso*bath_dens(iw)*dw
          S0gtr(i) =S0gtr(i) + xi*Vpd**2*ngtr*peso*bath_dens(iw)*dw
       enddo
    enddo
    call splot("BATH/S0less_t.ipt",t,S0less)
    call splot("BATH/S0gtr_t.ipt",t,S0gtr)
    call splot("DOSbath.lattice",wfreq,bath_dens)
    call system("mv *.lattice BATH/ 2>/dev/null")
  end subroutine get_Bath



  !********************************************************************
  !********************************************************************
  !********************************************************************



  !+-----------------------------------------------------------------+
  !PURPOSE  : Build BATH dispersion arrays \epsilon_bath(\ka) = bath_epski(i)
  !+-----------------------------------------------------------------+
  subroutine get_bath_constant_dos()
    integer    :: i,ik
    real(8)    :: w
    complex(8) :: gf,iw,zeta

    do i=-Lw,Lw
       w=wfreq(i) ; zeta=cmplx(w,eps,8)       
       gf=zero
       do ik=1,Lmu
          gf=gf+one/(zeta-bath_epsik(ik))/dble(Lmu)
       enddo
       bath_dens(i)=-aimag(gf)/pi
    enddo

  end subroutine get_bath_constant_dos



  !********************************************************************
  !********************************************************************
  !********************************************************************





  !+-----------------------------------------------------------------+
  !PURPOSE  : Build BATH dispersion arrays \epsilon_bath(\ka) = bath_epski(i)
  !+-----------------------------------------------------------------+
  subroutine get_bath_gaussian_dos()
    integer    :: i,ik
    real(8)    :: w,sig,alpha
    complex(8) :: gf,zeta

    bath_dens = exp(-0.5d0*(wfreq/Wbath)**2)/(sqrt(pi2)*Wbath) !standard Gaussian
    !bath_dens = exp(-(wfreq/Wbath)**2)/(sqrt(pi)*Wbath) !Camille's choice

    !    !!w/ erf in frquency space: coded from Abramowitz-Stegun
    ! do i=-Lw,Lw
    !    w=wfreq(i)
    !    !zeta=cmplx(w,eps,8)
    !    !sig=aimag(zeta)/abs(dimag(zeta))
    !    !gf=-sig*xi*sqrt(pi)*wfun(zeta/Wbath)/Wbath
    !    !bath_dens(i)=-aimag(gf)/pi
    ! enddo
  end subroutine get_bath_gaussian_dos




  !********************************************************************
  !********************************************************************
  !********************************************************************




  !+-----------------------------------------------------------------+
  !PURPOSE  : Build BATH dispersion arrays \epsilon_bath(\ka) = bath_epski(i)
  !+-----------------------------------------------------------------+
  subroutine get_bath_bethe_dos()
    integer    :: i,ik
    real(8)    :: w,sig,alpha
    complex(8) :: gf,zeta

    do i=-Lw,Lw
       w=wfreq(i)
       zeta=cmplx(w,eps,8)
       gf=gfbether(w,zeta,wbath/2.d0)
       bath_dens(i)=-aimag(gf)/pi
    enddo

  end subroutine get_bath_bethe_dos
  !********************************************************************
  !********************************************************************
  !********************************************************************









  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : Build the Bath part of the system using frequency domain
  ! !formalism. 
  ! !COMMENT  : Construct the Bath-Sigma functions from freq. domain
  ! !(that is 'cause you know the analytic expression, same
  ! !procedure as in IPTkeldysh @equilibrium)
  ! !This turns out to be the same provided one takes 
  ! !the limit \eta\goto 0.d0
  ! !+-------------------------------------------------------------------+
  ! subroutine get_Bath_w
  !   integer :: i,im
  !   real(8) :: A,An,w
  !   complex(8) :: iw,fg
  !   complex(8),dimension(-nstep:nstep) :: gb0fless,gb0fgtr
  !   complex(8),dimension(-nstep:nstep) :: gb0tgtr,gb0tless
  !   print '(A)',"Get Bath:"

  !   do i=-nstep,nstep
  !      w = wr(i) ;iw= cmplx(w,eps,8)
  !      fg=zero
  !      do im=0,Lmu
  !         fg=fg + de*bath_dens(im)/(iw-bath_epsik(im))!/dble(Lmu)
  !      enddo
  !      A = -aimag(fg)/pi
  !      An= A*fermi0(w,beta)
  !      gb0fless(i)= pi2*xi*An
  !      gb0fgtr(i) = pi2*xi*(An-A)
  !   enddo
  !   call fftgf_rw2rt(gb0fless,gb0tless,nstep)  ; gb0tless=fmesh/pi2*gb0tless
  !   call fftgf_rw2rt(gb0fgtr, gb0tgtr,nstep)   ; gb0tgtr =fmesh/pi2*gb0tgtr
  !   gb0tgtr =exa*gb0tgtr
  !   gb0tless=exa*gb0tless
  !   !Get the Self-energy contribution of the bath:
  !   !=================================================
  !   S0gtr=zero
  !   S0less=zero
  !   S0less=Vpd**2*gb0tless
  !   S0gtr =Vpd**2*gb0tgtr
  !   call system("mkdir BATH")   
  !   call splot("BATH/bathG0less_t.ipt",t(-nstep:nstep),gb0tless)
  !   call splot("BATH/bathG0gtr_t.ipt",t(-nstep:nstep),gb0tgtr)
  !   call splot("BATH/S0less_t.ipt",t(-nstep:nstep),S0less)
  !   call splot("BATH/S0gtr_t.ipt",t(-nstep:nstep),S0gtr)
  !   print '(A)',"" 
  !   return
  ! end subroutine Get_Bath_w



  !********************************************************************
  !********************************************************************
  !********************************************************************




end module BATH

