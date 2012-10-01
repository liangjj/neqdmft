!###############################################################
!     PURPOSE  : A non-equilibrium IPT solver module. 
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_NEQ
  USE VARS_GLOBAL
  USE EQUILIBRIUM
  USE ELECTRIC_FIELD
  USE MATRIX
  implicit none
  private

  public  :: neq_init_run
  public  :: neq_solve_ipt

  integer :: Liw,Lw
  real(8),allocatable,dimension(:)    :: wr_,wm_

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Initialize the run guessing/reading/setting initial conditions
  !+-------------------------------------------------------------------+
  subroutine neq_init_run()
    integer                          :: i,j,ik
    integer                          :: iselect,irdL
    real(8)                          :: en,intE,A,fmesh_
    real(8)                          :: nless,ngtr,xmu_,beta_
    complex(8)                       :: peso
    logical                          :: checkS,checkG0,checkNk
    real(8) :: xgrid(2*Nstep)
    complex(8) :: G0ret_t(-nstep:nstep),G0ret_w(2*nstep),G0mat(0:Nstep,0:Nstep),G0tmp(0:nstep,0:nstep)

    call create_data_dir("InitialConditions")

    !Check if n(k) file exists.
    inquire(file=trim(irdnkfile),exist=checkNk)
    if(.not.checkNk)inquire(file=trim(irdNkfile)//".gz",exist=checkNk)
    if(checkNk)then
       allocate(eq_nk(Lk))
       call read_nkfile(eq_nk,trim(irdnkfile))
    else
       !Get non-interacting n(k):
       xmu_=xmu   ; if(iquench)xmu_=xmu0
       beta_=beta ; if(iquench)beta_=beta0
       allocate(eq_nk(Lk))
       do ik=1,Lk
          eq_nk(ik)=fermi0((epsik(ik)-xmu_),beta_)
       enddo
    endif
    if(mpiID==0)call splot("InitialConditions/ic_nkVSepsk.ipt",epsik,eq_nk)

    !Restart from a previous solution: check if Sigma^<,> file exists.
    checkS=inquire_keldysh_contour_gf(trim(irdSfile))

    if(checkS)then
       call msg("Reading self-energy guess from input file.",lines=1,id=0)
       call read_keldysh_contour_gf(Sigma,trim(irdSfile))

    else  !DEFAULT: no files read, start from non-interacting HF solution or G0_loc if required       

       if(.not.g0loc_guess)then
          call msg("Using Hartree-Fock for self-energy guess",id=0)
          call msg("G0less=G0gtr=zero",lines=1,id=0)
          G0=zero

       elseif(g0loc_guess)then
          if(equench)then
             call msg("Using G0_loc + electric field for self-energy guess",lines=1,id=0)
             do ik=1,Lk
                en   = epsik(ik)
                nless= fermi0(en,beta)
                ngtr = fermi0(en,beta)-1.d0
                do j=0,nstep
                   do i=0,nstep
                      intE=int_Ht(ik,i,j)
                      peso=exp(-xi*intE)
                      G0%less(i,j)= G0%less(i,j) + xi*nless*peso*wt(ik)
                      G0%gtr(i,j) = G0%gtr(i,j)  + xi*ngtr*peso*wt(ik)
                   enddo
                enddo
             enddo
          else
             call msg("Using G0_loc for self-energy guess",lines=1,id=0)
             do ik=1,Lk
                en   = epsik(ik)
                nless= fermi0(en,beta)
                ngtr = fermi0(en,beta)-1.d0
                A    = wt(ik)
                do i=0,nstep
                   do j=0,nstep
                      peso=exp(-xi*en*(t(i)-t(j)))
                      G0%less(i,j)=G0%less(i,j) + xi*nless*A*peso
                      G0%gtr(i,j) =G0%gtr(i,j)  + xi*ngtr*A*peso
                   enddo
                enddo
             enddo
          endif
       endif

       if(mpiID==0)then
          call write_keldysh_contour_gf(G0,"InitialConditions/guessG0")
          if(plot3D)call plot_keldysh_contour_gf(G0,t(0:),"PLOT/guessG0")
       endif

       call neq_solve_ipt()

    endif


    ! !Start from a non HF guess given the BATH: check if G0(w) file exists.
    ! inquire(file=trim(irdG0wfile),exist=checkG0)
    ! if(.not.checkG0)inquire(file=trim(irdG0wfile)//".gz",exist=checkG0)
    ! if(checkS)then              !S^<,>(t,t') file exists
    !    iselect=1
    ! elseif(checkG0)then         !G0(w) file exists
    !    iselect=2
    ! endif
    ! select case(iselect)
    !    !
    ! case default 
    ! case(2)
    !    call msg("Reading G0(w) from input file.",id=0)
    !    call msg("Using G0(w) to guess the self-energy.",lines=1,id=0)
    !    !
    !    Lw=file_length(trim(irdG0wfile))
    !    allocate(eq_G0w(Lw),wr_(Lw))
    !    call sread(trim(irdG0wfile),wr_,eq_G0w)
    !    fmesh_=abs(wr_(2)-wr_(1))
    !    !
    !    G0=zero
    !    do ik=1,Lw
    !       en   = wr_(ik)
    !       nless= fermi0(en,beta)
    !       ngtr = fermi0(en,beta)-1.d0
    !       A    = -aimag(eq_G0w(ik))/pi*fmesh_
    !       do i=0,nstep
    !          do j=0,nstep
    !             peso=exp(-xi*en*(t(i)-t(j)))
    !             G0%less(i,j)=G0%less(i,j) + xi*nless*A*peso
    !             G0%gtr(i,j) =G0%gtr(i,j)  + xi*ngtr*A*peso
    !          enddo
    !       enddo
    !    enddo
    !    deallocate(wr_)
    !    if(mpiID==0)then
    !       call write_keldysh_contour_gf(G0,"InitialConditions/guessG0")
    !       if(plot3D)call plot_keldysh_contour_gf(G0,t(0:),"PLOT/guessG0")
    !    endif
    !    call neq_solve_ipt()
    ! end select


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

    subroutine read_nkfile(irdnk,file)
      character(len=*)     :: file
      real(8),dimension(Lk):: irdnk
      integer              :: redLk
      real(8),allocatable  :: rednk(:),redek(:)
      integer,allocatable  :: orderk(:)
      real(8),allocatable  :: uniq_rednk(:),uniq_redek(:)
      logical,allocatable  :: maskk(:)
      !n(k): A lot of work here to reshape the array
      redLk=file_length(file)
      allocate(rednk(redLk),redek(redLk),orderk(redLk))
      call sread(file,redek,rednk)
      !
      !work on the read arrays:
      !1 - sorting: sort the energies (X-axis), mirror on occupation (Y-axis) 
      !2 - delete duplicates energies (X-axis), mirror on occupation (Y-axis) 
      !3 - interpolate to the actual lattice structure (epsik,nk)
      call sort_array(redek,orderk)
      call reshuffle(rednk,orderk)
      call uniq(redek,uniq_redek,maskk)
      allocate(uniq_rednk(size(uniq_redek)))
      uniq_rednk = pack(rednk,maskk)
      call linear_spline(uniq_rednk,uniq_redek,irdnk,epsik)
    end subroutine read_nkfile

  end subroutine neq_init_run




  !+-------------------------------------------------------------------+
  !PURPOSE  : BUild the 2^nd IPT sigma functions:
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt()
    integer      :: i,j,itau

    !Get SIgma:
    call msg("Get Sigma(t,t')")

    Sigma=zero
    forall(i=0:nstep,j=0:nstep)
       Sigma%gtr (i,j) = U**2*(G0%gtr(i,j)**2)*G0%less(j,i)
       Sigma%less(i,j) = U**2*(G0%less(i,j)**2)*G0%gtr(j,i)
    end forall

    ! !Get impurity GF and use SPT method if required
    ! call get_impuritygf
    ! if(method=="spt")then
    !    call msg("Recalculate Sigma using SPT")
    !    forall(i=0:nstep,j=0:nstep)
    !       Sig%less(i,j) = (U**2)*(impG%less(i,j)**2)*impG%gtr(j,i)
    !       Sig%gtr (i,j) = (U**2)*(impG%gtr(i,j)**2)*impG%less(j,i)                
    !    end forall
    ! end if

    !Save data:
    if(mpiID==0)then
       call write_keldysh_contour_gf(Sigma,reg_filename(data_dir)//"/Sigma")
       if(plot3D)call plot_keldysh_contour_gf(Sigma,t(0:),"PLOT/Sigma")
    endif

  end subroutine Neq_solve_ipt



  !********************************************************************
  !********************************************************************
  !********************************************************************




  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : evaluate the impurity neq Green's functions
  ! !+-------------------------------------------------------------------+
  ! subroutine get_impuritygf()
  !   integer                               :: i,j
  !   real(8)                               :: A,w
  !   complex(8),dimension(0:nstep,0:nstep) :: Uno,GammaRet,Gamma0Ret
  !   complex(8),dimension(0:nstep,0:nstep) :: dG0ret,dGret,dSret
  !   if(update_wfftw)then
  !      call get_equilibrium_impuritygf !not tested!
  !   else
  !      dSret=zero ; dG0ret=zero ; dGret=zero
  !      GammaRet=zero ; Gamma0Ret=zero
  !      !1 - get the Ret components of G_0 && \Sigma:
  !      forall(i=0:nstep,j=0:nstep)
  !         dG0ret(i,j)=heaviside(t(i)-t(j))*(G0gtr(i,j) - G0less(i,j))
  !         dSret(i,j) =heaviside(t(i)-t(j))*(Sgtr(i,j) - Sless(i,j))
  !      end forall
  !      !2 - get the  operator: \Gamma_0^R = \Id - \Sigma^R\circ G_0^R && invert it
  !      Uno=zero  ; forall(i=0:nstep)Uno(i,i)=One/dt
  !      Gamma0Ret(0:nstep,0:nstep) = Uno-matmul(dSret(0:nstep,0:nstep),dG0ret(0:nstep,0:nstep))*dt
  !      Gamma0Ret(0:nstep,0:nstep)=Gamma0Ret(0:nstep,0:nstep)*dt**2
  !      call mat_inversion_GJ(Gamma0Ret(0:nstep,0:nstep))
  !      !3 - get G_imp^R, G_imp^{>,<} using Dyson equations:
  !      dGret(0:nstep,0:nstep)    = matmul(dG0ret(0:nstep,0:nstep),Gamma0Ret(0:nstep,0:nstep))*dt 
  !      GammaRet(0:nstep,0:nstep) = Uno + matmul(dGret(0:nstep,0:nstep),dSret(0:nstep,0:nstep))*dt

  !      impGless(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),&
  !           matmul(G0less(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep)))))*dt**2 +&
  !           matmul(dGret(0:nstep,0:nstep),matmul(Sless(0:nstep,0:nstep),conjg(transpose(dGret(0:nstep,0:nstep)))))*dt**2

  !      impGgtr(0:nstep,0:nstep)  = matmul(GammaRet(0:nstep,0:nstep),&
  !           matmul(G0gtr(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep)))))*dt**2  +&
  !           matmul(dGret(0:nstep,0:nstep),matmul(Sgtr(0:nstep,0:nstep),conjg(transpose(dGret(0:nstep,0:nstep)))))*dt**2
  !   endif
  !   !Save data:
  !   if(mpiID==0)then
  !      call splot("impGless.data",impG%less(0:nstep,0:nstep))
  !      call splot("impGgtr.data",impG%gtr(0:nstep,0:nstep))
  !   endif
  ! end subroutine get_impuritygf





  !********************************************************************
  !********************************************************************
  !********************************************************************


end module IPT_NEQ
