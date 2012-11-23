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
  integer,parameter                :: Lw=2048 !# of frequencies
  real(8),allocatable,dimension(:) :: bath_dens,wfreq

  public                           :: get_thermostat_bath

contains
  !+-------------------------------------------------------------------+
  !PURPOSE  : Build the Bath part of the system using exclusively time
  !dependent formalism. The bath is not interacting so it is 
  ! time translation invariant.
  !+-------------------------------------------------------------------+
  subroutine get_thermostat_bath()
    integer          :: iw,i,j,itau
    real(8)          :: en,w,dw,wfin,wini
    complex(8)       :: peso
    real(8)          :: ngtr,nless,arg

    call msg("Get Bath. Type: "//bold_green(trim(adjustl(trim(bath_type))))//" dissipative bath",id=0)
    call create_data_dir("Bath")
    allocate(bath_dens(Lw),wfreq(Lw))

    select case(trim(adjustl(trim(bath_type))))
    case("bethe")
       wfin  = 2.d0*Wbath ; wini=-wfin
       wfreq = linspace(wini,wfin,Lw,mesh=dw)
       call get_bath_bethe_dos()

    case("gaussian")
       wfin = 4.d0*Wbath ; wini=-wfin
       wfreq= linspace(wini,wfin,Lw,mesh=dw)
       call get_bath_gaussian_dos()

    case ("constant")
       wfin  = 2.d0*Wbath ; wini=-wfin
       wfreq = linspace(wini,wfin,Lw,mesh=dw)
       call get_bath_constant_dos()

    case default
       call abort("Bath type:"//trim(adjustl(trim(bath_type)))//" not supported. Accepted values are: constant,gaussian,bethe.")

    end select


    ! S0less=zero ; S0gtr=zero
    ! S0lmix=zero; S0gmix=zero
    S0 = zero
    if(Vbath/=0.d0)then
       do iw=1,Lw
          en   = wfreq(iw)
          nless= fermi0(en,beta)
          ngtr = fermi0(en,beta)-1.d0 !it absorbs the minus sign of the greater functions
          do i=0,nstep
             do j=0,nstep
                peso=exp(-xi*(t(i)-t(j))*en)
                S0%less(i,j)=S0%less(i,j)+ xi*Vbath**2*nless*peso*bath_dens(iw)*dw
                S0%gtr(i,j) =S0%gtr(i,j) + xi*Vbath**2*ngtr*peso*bath_dens(iw)*dw
             enddo

             if(en>=0.d0)then
                do itau=0,Ltau
                   peso=exp(-xi*en*t(i))*exp(-en*(beta-tau(itau)))/(1.d0+exp(-en*beta))
                   if(beta*en>20.d0)peso=exp(-xi*en*t(i))*exp(-en*(beta-tau(itau)))
                   S0%lmix(i,itau) = S0%lmix(i,itau) + xi*Vbath**2*peso*bath_dens(iw)*dw
                enddo
             else
                do itau=0,Ltau
                   peso=exp(-xi*en*t(i))*exp(-en*(beta-tau(itau)))/(1.d0+exp(-en*beta))
                   if(beta*en<-20.d0)peso=exp(-xi*en*t(i))*exp(en*tau(itau))
                   S0%lmix(i,itau) = S0%lmix(i,itau) + xi*Vbath**2*peso*bath_dens(iw)*dw
                enddo
             endif
          enddo
       enddo
       forall(i=0:Ltau)S0%gmix(i,:)=conjg(S0%lmix(:,Ltau-i))
    endif

    if(mpiID==0)then
       call splot("Bath/DOSbath.lattice",wfreq,bath_dens)
       if(Vbath/=0.d0 .AND. plot3D)call plot_kbm_contour_gf(S0,t(0:),tau(0:),"Bath/S0")          
       ! call splot("Bath/S0less_t.ipt",t,S0less)
       ! call splot("Bath/S0gtr_t.ipt",t,S0gtr)
       ! call splot("Bath/S0lmix_t_tau",t(0:nstep),tau(0:Ltau),S0lmix(0:nstep,0:Ltau))
    endif

  end subroutine get_thermostat_bath



  !********************************************************************
  !********************************************************************
  !********************************************************************



  !+-----------------------------------------------------------------+
  !PURPOSE  : Build constant BATH 
  !+-----------------------------------------------------------------+
  subroutine get_bath_constant_dos()
    integer    :: i
    real(8)    :: w

    do i=1,Lw
       w=wfreq(i)
       bath_dens(i)= heaviside(Wbath-abs(w))/(2.d0*Wbath)
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
    !bath_dens = exp(-((wfreq)/Wbath)**2)/(sqrt(pi)*Wbath) !Camille's choice

    !    !!w/ erf in frquency space: coded from Abramowitz-Stegun
    ! do i=-Lw,Lw
    !    !w=wfreq(i)
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

    do i=1,Lw
       w=wfreq(i)
       zeta=cmplx(w,eps,8)
       gf=gfbether(w,zeta,wbath/2.d0)
       bath_dens(i)=-aimag(gf)/pi
    enddo

  end subroutine get_bath_bethe_dos
  !********************************************************************
  !********************************************************************
  !********************************************************************



end module BATH

