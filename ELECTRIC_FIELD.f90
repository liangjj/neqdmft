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
  !PURPOSE: set the electric field vector using given direction
  !+------------------------------------------------------------+
  function set_efield_vector(Ex,Ey) result(E)
    type(vect2D)  :: E
    real(8)       :: Ex,Ey
    real(8)       :: modulo
    integer       :: i

    modulo=sqrt(Ex**2+Ey**2)
    if(modulo/=1)then
       Ex=Ex/modulo
       Ey=Ey/modulo
    endif
    E%x=Ex;E%y=Ey
    if(alat==0)then
       print*, "a_lat=0! EXIT"
       stop
    endif
    E=Efield*E                  !Normalized Electric field vector    
  end function set_efield_vector


  !+------------------------------------------------------------+
  !PURPOSE : 
  !+------------------------------------------------------------+
  pure function Afield(t,E)
    type(vect2D),intent(in) :: E
    real(8),intent(in)      :: t
    real(8)                 :: field,tau1
    type(vect2D)            :: Afield
    complex(8)              :: zp,zm

    select case(field_profile)
    case default                !!CONSTANT ELECTRIC FIELD:
       field=-(step(t-t0)*(t-t0 + (t1-t)*step(t-t1) - (t1-t0)*step(t0-t1)))

    case("pulse")               !!LIGHT PULSE (for Pump&Probe)
       tau1=tau0/pi2
       zp=cmplx(t-t0,tau1**2*w0,8)/(sqrt(2.d0)*tau1)
       zm=cmplx(t-t0,-tau1**2*w0,8)/(sqrt(2.d0)*tau1)
       field =-real(sqrt(pi/2.d0)/2.d0*tau1*exp(-(tau1*w0)**2/2.d0)*(zerf(zm)+zerf(zp)),8)


    case("ac")
       field=-sin(Omega0*(t-t0))/Omega0

    case("ac1")
       field=-(t+sin(Omega0*(t-t0))/Omega0)

    case("ramp")                !!RAMP CONSTANT FIELD:
       field=-(24.d0*pi*(t+(t-t0)*step(t-t0)+2.d0*(t1-t)*step(t-t0)*step(t-t1)-&
            2.d0*(t0-t1)*step(t-t0)*step(t0-t1))+                              &
            27.d0*t0*(step(t-t0)-1.d0)*Sin(pi*t/t0) - t0*(step(t-t0)-1.d0)*Sin(3.d0*pi*t/t0))/48.d0/pi

       !!add more here:
    end select

    !-----------------------------
    Afield=E*field
  end function Afield
  !***************************************************************
  !***************************************************************
  !***************************************************************





  subroutine print_Afield_form(t)
    real(8),dimension(:) :: t
    integer              :: i
    open(10,file="Avector_shape.ipt")
    do i=1,size(t)
       write(10,*)t(i),Afield(t(i),Ek)
    enddo
    close(10)
  end subroutine print_Afield_form


  !***************************************************************
  !***************************************************************
  !***************************************************************


end MODULE ELECTRIC_FIELD
