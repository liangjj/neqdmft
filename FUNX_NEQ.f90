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
    complex(8),dimension(0:nstep,0:nstep) :: locGret,Sret
    complex(8),dimension(0:nstep,0:nstep) :: locGadv,Sadv
    complex(8),dimension(0:nstep,0:nstep) :: Uno,GammaRet
    complex(8),dimension(0:nstep,0:nstep) :: G0ret
    complex(8),dimension(0:nstep,0:nstep) :: G0adv
    !
    complex(8),dimension(:,:),allocatable :: mat_Delta
    complex(8),dimension(:,:),allocatable :: mat_Gamma
    complex(8),dimension(:,:),allocatable :: mat_G0,mat_Sigma,mat_locG
    !
    type(keldysh_contour_gf),save             :: G0_old

    if(G0_old%status.EQV..false.)call allocate_keldysh_contour_gf(G0_old,Nstep)
    G0_old=G0    

    call msg("Update WF: Dyson")
    if(update_wfftw)then
       call update_equilibrium_weiss_field
    else
       include "update_G0_nonequilibrium.f90"
    endif
    G0%less = weight*G0%less + (1.d0-weight)*G0_old%less
    G0%gtr  = weight*G0%gtr  + (1.d0-weight)*G0_old%gtr

    !Save data:
    if(mpiID==0)then
       call write_keldysh_contour_gf(G0,reg_filename(data_dir)//"/G0")
       if(plot3D)call plot_keldysh_contour_gf(G0,t(0:),"PLOT/G0")
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






  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : Build a guess for the initial Weiss Fields G0^{<,>} 
  ! !as non-interacting GFs. If required read construct it starting 
  ! !from a read seed (cf. routine for seed reading)
  ! !+-------------------------------------------------------------------+
  ! subroutine neq_guess_weiss_field
  !   integer    :: i,j,ik,redLk
  !   real(8)    :: en,intE,A
  !   complex(8) :: peso
  !   real(8)    :: nless,ngtr

  !   call msg("Get G0guess(t,t')",id=0)
  !   gf0=zero ; G0gtr=zero ; G0less=zero
  !   if(mpiID==0)then

  !      if(irdeq .OR. solve_eq)then            !Read from equilibrium solution
  !         call read_init_seed()
  !         do ik=1,irdL             !2*L
  !            en   = irdwr(ik)
  !            nless= fermi0(en,beta)
  !            ngtr = fermi0(en,beta)-1.d0
  !            A    = -aimag(irdG0w(ik))/pi*irdfmesh
  !            do i=-nstep,nstep
  !               peso=exp(-xi*en*t(i))
  !               gf0%less%t(i)=gf0%less%t(i) + xi*nless*A*peso
  !               gf0%gtr%t(i) =gf0%gtr%t(i)  + xi*ngtr*A*peso
  !            enddo
  !         enddo
  !         forall(i=0:nstep,j=0:nstep)
  !            G0less(i,j)=gf0%less%t(i-j)
  !            G0gtr(i,j) =gf0%gtr%t(i-j)
  !         end forall

  !      else

  !         if(equench)then
  !            do ik=1,Lk
  !               en   = epsik(ik)
  !               nless= fermi0(en,beta)
  !               ngtr = fermi0(en,beta)-1.d0
  !               do j=0,nstep
  !                  do i=0,nstep
  !                     intE=int_Ht(ik,i,j)
  !                     peso=exp(-xi*intE)
  !                     G0less(i,j)= G0less(i,j) + xi*nless*peso*wt(ik)
  !                     G0gtr(i,j) = G0gtr(i,j)  + xi*ngtr*peso*wt(ik)
  !                  enddo
  !               enddo
  !            enddo

  !         else

  !            do ik=1,Lk
  !               en   = epsik(ik)
  !               nless= fermi0(en,beta)
  !               ngtr = fermi0(en,beta)-1.d0
  !               A    = wt(ik)
  !               do i=-nstep,nstep
  !                  peso=exp(-xi*en*t(i))
  !                  gf0%less%t(i)=gf0%less%t(i) + xi*nless*A*peso
  !                  gf0%gtr%t(i) =gf0%gtr%t(i)  + xi*ngtr*A*peso
  !               enddo
  !            enddo
  !            forall(i=0:nstep,j=0:nstep)
  !               G0less(i,j)=gf0%less%t(i-j)
  !               G0gtr(i,j) =gf0%gtr%t(i-j)
  !            end forall

  !         endif
  !      endif

  !      call splot("guessG0less.data",G0less(0:nstep,0:nstep))
  !      call splot("guessG0gtr.data",G0gtr(0:nstep,0:nstep))
  !   endif
  !   call MPI_BCAST(G0less,(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  !   call MPI_BCAST(G0gtr,(nstep+1)*(nstep+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)

  ! contains

  !   function int_Ht(ik,it,jt)
  !     real(8)      :: int_Ht
  !     integer      :: i,j,ii,ik,it,jt,sgn
  !     type(vect2D) :: kt,Ak
  !     int_Ht=0.d0 ; if(it==jt)return
  !     sgn=1 ; if(jt > it)sgn=-1
  !     i=ik2ix(ik); j=ik2iy(ik)
  !     do ii=jt,it,sgn
  !        Ak=Afield(t(ii),Ek)
  !        kt=kgrid(i,j) - Ak
  !        int_Ht=int_Ht + sgn*square_lattice_dispersion(kt)*dt
  !     enddo
  !   end function int_Ht

  ! end subroutine neq_guess_weiss_field




  !********************************************************************
  !********************************************************************
  !********************************************************************



  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : Build a guess for the initial Weiss Fields G0^{<,>} 
  ! !as non-interacting GFs. If required construct it starting 
  ! !from a read seed (cf. routine for seed reading)
  ! !+-------------------------------------------------------------------+
  ! subroutine read_init_seed()
  !   logical :: control
  !   real(8) :: w1,w2
  !   integer :: ik,redLk
  !   real(8),allocatable :: rednk(:),redek(:)
  !   integer,allocatable :: orderk(:)
  !   real(8),allocatable :: uniq_rednk(:),uniq_redek(:)
  !   logical,allocatable :: maskk(:)

  !   !GO_realw:
  !   inquire(file=trim(irdG0file),exist=control)
  !   if(.not.control)call error("Can not find irdG0file")
  !   !Read the function WF
  !   irdL=file_length(trim(irdG0file))
  !   allocate(irdG0w(irdL),irdwr(irdL))
  !   call sread(trim(irdG0file),irdwr,irdG0w)
  !   !Get G0 mesh:
  !   irdfmesh=abs(irdwr(2)-irdwr(1))

  !   !n(k): A lot of work here to reshape the array
  !   inquire(file=trim(irdnkfile),exist=control)
  !   if(.not.control)call abort("Can not find irdnkfile")
  !   !Read the function nk.
  !   redLk=file_length(trim(irdnkfile))
  !   allocate(rednk(redLk),redek(redLk),orderk(redLk))
  !   call sread(trim(irdnkfile),redek,rednk)
  !   !work on the read arrays:
  !   !1 - sorting: sort the energies (X-axis), mirror on occupation (Y-axis) 
  !   !2 - delete duplicates energies (X-axis), mirror on occupation (Y-axis) 
  !   !3 - interpolate to the actual lattice structure (epsik,nk)
  !   call sort_array(redek,orderk)
  !   call reshuffle(rednk,orderk)
  !   call uniq(redek,uniq_redek,maskk)
  !   allocate(uniq_rednk(size(uniq_redek)))
  !   uniq_rednk = pack(rednk,maskk)
  !   allocate(irdnk(Lk))
  !   call linear_spline(uniq_rednk,uniq_redek,irdnk,epsik)

  !   !G0_iw:
  !   allocate(irdG0iw(L),irdG0tau(0:Ltau))
  !   call get_matsubara_gf_from_DOS(irdwr,irdG0w,irdG0iw,beta)
  !   call fftgf_iw2tau(irdG0iw,irdG0tau,beta)


  !   !Print out the initial conditions as effectively read from the files:
  !   call system("if [ ! -d InitialConditions ]; then mkdir InitialConditions; fi")
  !   call splot("InitialConditions/read_G0_realw.ipt",irdwr,irdG0w)
  !   call splot("InitialConditions/read_G0_iw.ipt",irdwm,irdG0iw)
  !   call splot("InitialConditions/read_G0_tau.ipt",tau,irdG0tau)
  !   call splot("InitialConditions/read_nkVSek.ipt",epsik,irdnk)

  !   call MPI_BCAST(irdG0w,irdL,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  !   call MPI_BCAST(irdG0iw,irdL,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  !   call MPI_BCAST(irdG0tau,Ltau+1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  !   call MPI_BCAST(irdNk,Lk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  ! end subroutine read_init_seed




  !********************************************************************
  !********************************************************************
  !********************************************************************





end module FUNX_NEQ
