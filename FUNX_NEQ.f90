!###############################################################
!     PURPOSE  : Constructs some functions used in other places. 
!     AUTHORS  : Adriano Amaricci
!###############################################################
module FUNX_NEQ
  USE MATRIX
  USE VARS_GLOBAL
  USE ELECTRIC_FIELD
  USE BATH
  implicit none
  private
  !
  integer                          :: irdL
  real(8)                          :: irdfmesh
  real(8),allocatable,dimension(:) :: irdwr
  !
  public                           :: guess_g0_sigma
  public                           :: update_weiss_field
  public                           :: get_sigma
  public                           :: obtain_gimp
  public                           :: print_observables,convergence_check

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Build a guess for the initial Weiss Fields G0^{<,>} 
  !as non-interacting GFs. If required read construct it starting 
  !from a read seed (cf. routine for seed reading)
  !+-------------------------------------------------------------------+
  subroutine read_init_seed()
    logical :: control
    real(8) :: w1,w2
    integer :: ik,redLk
    real(8),allocatable :: rednk(:),redek(:)
    integer,allocatable :: orderk(:)
    real(8),allocatable :: uniq_rednk(:),uniq_redek(:)
    logical,allocatable :: maskk(:)

    !GOret:
    !Check if G0file exist
    inquire(file=trim(irdG0file),exist=control)
    if(.not.control)call abort("Can not find irdG0file")
    !Get file length to allocate arrays
    irdL=file_length(trim(irdG0file))
    !Read the function: the WF
    allocate(irdG0w(irdL))
    allocate(irdwr(irdL))
    call sread(trim(irdG0file),irdwr,irdG0w)
    !Get G0 mesh:
    irdfmesh=abs(irdwr(2)-irdwr(1))

    !n(k):
    inquire(file=trim(irdnkfile),exist=control)
    if(.not.control)call abort("Can not find irdnkfile")
    !Get file length:
    redLk=file_length(trim(irdnkfile))
    allocate(rednk(redLk),redek(redLk),orderk(redLk))
    call sread(trim(irdnkfile),redek,rednk)
    !work on the read arrays:
    !1 - sorting: sort the energies (X-axis), mirror on occupation (Y-axis) 
    !2 - delete duplicates energies (X-axis), mirror on occupation (Y-axis) 
    !3 - interpolate to the actual lattice structure (epsik,nk)
    call sort_array(redek,orderk)
    call reshuffle(rednk,orderk)
    call uniq(redek,uniq_redek,maskk)
    allocate(uniq_rednk(size(uniq_redek)))
    uniq_rednk = pack(rednk,maskk)
    allocate(irdnk(Lk))
    call linear_spline(uniq_rednk,uniq_redek,irdnk,epsik)

    call system("if [ ! -d InitialConditions ]; then mkdir InitialConditions; fi")
    call splot("InitialConditions/read_G0_realw.ipt",irdwr,irdG0w)
    call splot("InitialConditions/read_nkVSek.ipt",epsik,irdnk)
  end subroutine read_init_seed




  !********************************************************************
  !********************************************************************
  !********************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : Build a guess for the initial Weiss Fields G0^{<,>} 
  !as non-interacting GFs. If required read construct it starting 
  !from a read seed (cf. routine for seed reading)
  !+-------------------------------------------------------------------+
  subroutine guess_g0_sigma
    integer    :: i,j,ik,redLk
    real(8)    :: en,intE,A
    complex(8) :: peso
    real(8)    :: nless,ngtr

    call msg("Get G0guess(t,t')",id=0)
    gf0=zero ; G0gtr=zero ; G0less=zero
    if(mpiID==0)then

       if(irdeq)then            !Read from equilibrium solution
          call read_init_seed()
          do ik=1,irdL             !2*L
             en   = irdwr(ik)
             nless= fermi0(en,beta)
             ngtr = fermi0(en,beta)-1.d0
             A    = -aimag(irdG0w(ik))/pi*irdfmesh
             do i=-nstep,nstep
                peso=exp(-xi*en*t(i))
                gf0%less%t(i)=gf0%less%t(i) + xi*nless*A*peso
                gf0%gtr%t(i) =gf0%gtr%t(i)  + xi*ngtr*A*peso
             enddo
          enddo
          forall(i=0:nstep,j=0:nstep)
             G0less(i,j)=gf0%less%t(i-j)
             G0gtr(i,j) =gf0%gtr%t(i-j)
          end forall

       else

          if(equench)then
             do ik=1,Lk
                en   = epsik(ik)
                nless= fermi0(en,beta)
                ngtr = fermi0(en,beta)-1.d0
                do j=0,nstep
                   do i=0,nstep
                      intE=int_Ht(ik,i,j)
                      peso=exp(-xi*intE)
                      G0less(i,j)= G0less(i,j) + xi*nless*peso*wt(ik)
                      G0gtr(i,j) = G0gtr(i,j)  + xi*ngtr*peso*wt(ik)
                   enddo
                enddo
             enddo

          else

             do ik=1,Lk
                en   = epsik(ik)
                nless= fermi0(en,beta)
                ngtr = fermi0(en,beta)-1.d0
                A    = wt(ik)
                do i=-nstep,nstep
                   peso=exp(-xi*en*t(i))
                   gf0%less%t(i)=gf0%less%t(i) + xi*nless*A*peso
                   gf0%gtr%t(i) =gf0%gtr%t(i)  + xi*ngtr*A*peso
                enddo
             enddo
             forall(i=0:nstep,j=0:nstep)
                G0less(i,j)=gf0%less%t(i-j)
                G0gtr(i,j) =gf0%gtr%t(i-j)
             end forall

          endif
       endif

       call splot("guessG0less.data",G0less(0:nstep,0:nstep))
       call splot("guessG0gtr.data",G0gtr(0:nstep,0:nstep))
    endif
    call MPI_BCAST(G0less,(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BCAST(G0gtr,(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)

    call get_sigma 

  contains

    function int_Ht(ik,it,jt)
      real(8)      :: int_Ht
      integer      :: i,j,ii,ik,it,jt,sgn
      type(vect2D) :: kt,Ak
      int_Ht=0.d0 ; if(it==jt)return
      sgn=1 ; if(jt > it)sgn=-1
      i=ik2ix(ik); j=ik2iy(ik)
      do ii=jt,it,sgn
         Ak=Afield(t(ii),Ek)
         kt=kgrid(i,j) - Ak
         int_Ht=int_Ht + sgn*square_lattice_dispersion(kt)*dt
      enddo
    end function int_Ht

  end subroutine guess_g0_sigma





  !********************************************************************
  !********************************************************************
  !********************************************************************





  !+-------------------------------------------------------------------+
  !PURPOSE  : BUild the 2^nd IPT sigma functions:
  !+-------------------------------------------------------------------+
  subroutine get_sigma()
    integer                               :: i,j,itau

    !Get SIgma:
    call msg("Get Sigma(t,t')")
    forall(i=0:nstep,j=0:nstep)
       Sless(i,j) = (U**2)*(G0less(i,j)**2)*G0gtr(j,i)
       Sgtr (i,j) = (U**2)*(G0gtr(i,j)**2)*G0less(j,i)
    end forall

    !Get impurity GF and use SPT method if required
    call obtain_Gimp

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
    return
  end subroutine Get_Sigma





  !********************************************************************
  !********************************************************************
  !********************************************************************







  !+-------------------------------------------------------------------+
  !PURPOSE  : Update the Weiss Fields G0^{>,<,Ret} using Dyson equation
  !+-------------------------------------------------------------------+
  subroutine update_weiss_field
    integer :: M,i,j,k,itau,jtau,NN
    real(8) :: R,deg
    real(8) :: w,A,An
    complex(8),dimension(0:nstep,0:nstep) :: Uno,GammaRet
    complex(8),dimension(0:nstep,0:nstep) :: G0ret,locGret,Sret
    complex(8),dimension(0:nstep,0:nstep) :: G0adv,locGadv,Sadv
    complex(8),dimension(0:nstep,0:nstep) :: dG0less,dG0gtr
    complex(8),dimension(0:nstep,0:nstep) :: G0kel,locGkel,Skel
    complex(8),dimension(0:nstep,0:nstep) :: locGtt,locGat,Stt,Sat
    complex(8),dimension(:,:),allocatable :: locGmat,Smat,G0mat,GammaMat,UnoMat

    call msg("Update WF: Dyson")
    if(update_wfftw)then
       include "update_G0_equilibrium.f90"
    else
       include "update_G0_nonequilibrium.f90"
    endif

    !Save data:
    call splot("G0less.data",G0less(0:nstep,0:nstep))
    call splot("G0gtr.data",G0gtr(0:nstep,0:nstep))
  end subroutine update_weiss_field





  !********************************************************************
  !********************************************************************
  !********************************************************************






  !+-------------------------------------------------------------------+
  !PURPOSE  : evaluate the impurity neq Green's functions
  !+-------------------------------------------------------------------+
  subroutine obtain_Gimp()
    integer                               :: i,j
    real(8)                               :: A,w
    complex(8),dimension(0:nstep,0:nstep) :: Uno,GammaRet,Gamma0Ret
    complex(8),dimension(0:nstep,0:nstep) :: dG0ret,dGret,dSret

    if(update_wfftw)then
       include "obtain_Gimp_equilibrium.f90"
    else
       include "obtain_Gimp_nonequilibrium.f90"
    endif

    !Save data:
    if(mpiID==0)then
       call splot("impGless.data",impGless(0:nstep,0:nstep))
       call splot("impGgtr.data",impGgtr(0:nstep,0:nstep))
    endif

  end subroutine obtain_Gimp





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
       forall(i=0:nstep)nt(i)=-xi*locGless(i,i)
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
    if(mpiID==0)then
       if(Efield/=0.d0)then
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
          forall(i=0:nstep)test_func(i)=-xi*locGless(i,i)
       endif
       converged=check_convergence(test_func(0:nstep),eps_error,Nsuccess,nloop)
    endif
  end function convergence_check


end module FUNX_NEQ
