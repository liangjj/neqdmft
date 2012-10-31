!#####################################################################
!     Program  : ELECTRIC_FIELD
!     TYPE     : Module
!     PURPOSE  : Function for the Electric field vector potential
!     AUTHORS  : Adriano Amaricci
!#####################################################################
MODULE ELECTRIC_FIELD
  USE VARS_GLOBAL
  implicit none
  private

  public :: Afield
  public :: set_efield_vector
  public :: print_Afield_form


contains

  !+------------------------------------------------------------+
  !PURPOSE: set the normalized electric field versors using given direction
  !+------------------------------------------------------------+
  function set_efield_vector(Ex,Ey) result(E)
    type(vect2D)  :: E
    real(8)       :: Ex,Ey
    real(8)       :: modulo
    integer       :: i
    !Normalize the Electric Field components
    !Keep unaltered the Electric Field Strenght Efield=E0
    modulo=sqrt(Ex**2+Ey**2)
    Ex=Ex/modulo
    Ey=Ey/modulo
    E%x=Ex;E%y=Ey
    call msg("|E|=E0="//trim(txtfy(Efield/modulo)),id=0)
    if(alat==0)call error("a_lat=0! EXIT")
    select case(field_profile)
    case default
       stop "ELECTRIC_FIELD/Afield: wrong field_profile. set:dc,ac,acdc,pulse,ramp"
    case ("dc")                !DC ELECTRIC FIELD:
       return
    case("ac")                  !AC ELECTRIC FIELD
       return
    case("acdc")                !AC+DC ELECTRIC FIELD (super-bloch)
       return
    case("pulse")               !LIGHT PULSE (for Pump&Probe)
       return
    case("ramp")                !RAMP TO CONSTANT DC-FIELD:
       return
       !!add more here:
    end select
  end function set_efield_vector


  !+------------------------------------------------------------+
  !PURPOSE : 
  !+------------------------------------------------------------+
  pure function Afield(t,E)
    type(vect2D),intent(in) :: E
    real(8),intent(in)      :: t
    real(8)                 :: ftime,tau1
    type(vect2D)            :: Afield
    complex(8)              :: zp,zm

    select case(field_profile)
    case ("dc")                !DC ELECTRIC FIELD:
       ftime=-(step(t-t0)*(t-t0 + (t1-t)*step(t-t1) - (t1-t0)*step(t0-t1)))
       Afield=E*Efield*ftime       !A(t) = E0*F(t)*(e_x + e_y)

    case("ac")                  !AC ELECTRIC FIELD
       ftime=-sin(Omega0*(t-t0))/Omega0
       Afield=E*Efield*ftime       !A(t) = E0*F(t)*(e_x + e_y)

    case("acdc")                !AC+DC ELECTRIC FIELD (super-bloch)
       !ftime=-(t+sin(Omega0*(t-t0))/Omega0)
       ftime =-sin(Omega0*(t-t0))/Omega0
       Afield=E*(Efield*ftime - E1*t)       !A(t) = E0*F(t)*(e_x + e_y)

    case("pulse")               !LIGHT PULSE (for Pump&Probe)
       tau1=tau0/pi2
       zp=cmplx(t-t0,tau1**2*w0,8)/(sqrt(2.d0)*tau1)
       zm=cmplx(t-t0,-tau1**2*w0,8)/(sqrt(2.d0)*tau1)
       ftime =-real(sqrt(pi/2.d0)/2.d0*tau1*exp(-(tau1*w0)**2/2.d0)*(zerf(zm)+zerf(zp)),8)
       Afield=E*Efield*ftime       !A(t) = E0*F(t)*(e_x + e_y)

    case("ramp")                !RAMP TO CONSTANT DC-FIELD:
       ftime=-(24.d0*pi*(t+(t-t0)*step(t-t0)+2.d0*(t1-t)*step(t-t0)*step(t-t1)-&
            2.d0*(t0-t1)*step(t-t0)*step(t0-t1))+                              &
            27.d0*t0*(step(t-t0)-1.d0)*Sin(pi*t/t0) - &
            t0*(step(t-t0)-1.d0)*Sin(3.d0*pi*t/t0))/48.d0/pi
       Afield=E*Efield*ftime       !A(t) = E0*F(t)*(e_x + e_y)

       !!add more here:
    end select
    !-----------------------------



  end function Afield



  !***************************************************************
  !***************************************************************
  !***************************************************************





  subroutine print_Afield_form(t)
    real(8),dimension(:) :: t
    integer              :: i
    type(vect2D)         :: A
    real(8),dimension(size(t)) :: Ax,Ay,Ex,Ey
    do i=1,size(t)
       A=Afield(t(i),Ek)
       Ax(i)=A%x
       Ay(i)=A%y
    enddo
    Ex = deriv(Ax,dt)
    Ey = deriv(Ay,dt)
    open(10,file="Avector_shape.ipt")
    open(11,file="Efield_shape.ipt")
    do i=1,size(t)
       write(10,*)t(i),Afield(t(i),Ek)
       write(11,*)t(i),Ex(i),Ey(i)
    enddo
    close(10)
    close(11)
    if(field_profile=="ac")call msg("Root condition: "//trim(txtfy(bessel_j0(Efield/Omega0))))
  end subroutine print_Afield_form


  !***************************************************************
  !***************************************************************
  !***************************************************************


end MODULE ELECTRIC_FIELD
