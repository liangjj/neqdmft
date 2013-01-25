!###############################################################
!     PURPOSE  : A non-equilibrium IPT solver module. 
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_NEQ
  USE VARS_GLOBAL
  implicit none
  private

  public  :: neq_init_run
  public  :: neq_solve_ipt



contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Initialize the run guessing/reading self-energy
  !+-------------------------------------------------------------------+
  subroutine neq_init_run()
    logical :: init,checknk
    integer :: i,j,ik

    init = inquire_keldysh_contour_gf(trim(irdSFILE))
    if(init)then
       call msg(bold("Reading components of the input Self-energy and nk"))
       call read_keldysh_contour_gf(Sigma,trim(irdSFILE))
    else
       call msg(bold("Start from the Hartree-Fock self-energy"))
       Sigma=zero               !U*(n-0.5d0), n=0.5d0
    endif

    inquire(file=trim(irdnkfile),exist=checkNk)
    if(.not.checkNk)inquire(file=trim(irdNkfile)//".gz",exist=checkNk)
    if(checkNk)then
       call read_nkfile(eq_nk,trim(irdnkfile))
    else
       !Get non-interacting n(k):
       do ik=1,Lk
          eq_nk(ik)=fermi((epsik(ik))-xmu,beta)
       enddo
    endif
    if(mpiID==0)call splot("ic_nkVSepsk.ipt",epsik,eq_nk)

  contains

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
    real(8),dimension(0:nstep)            :: nt             !occupation(time)

    call msg("Get Sigma(t,t')")

    forall(i=0:nstep,j=0:nstep)
       Sigma%less(i,j) = (U**2)*(G0%less(i,j)**2)*G0%gtr(j,i)    !get Sigma^<(t,t`)
       Sigma%gtr (i,j) = (U**2)*(G0%gtr(i,j)**2)*G0%less(j,i)    !get Sigma^>(t,t`)
    end forall

    !Save data:
    if(mpiID==0)then
       call write_keldysh_contour_gf(Sigma,trim(data_dir)//"/Sigma")
       if(plot3D)call plot_keldysh_contour_gf(Sigma,t(0:),trim(plot_dir)//"/Sigma")
       forall(i=0:nstep)nt(i)=-xi*Sigma%less(i,i)
       call splot("nsVStime.ipt",t(0:),nt(0:),append=TT)
    endif
  end subroutine neq_solve_ipt



  !********************************************************************
  !********************************************************************
  !********************************************************************






end module IPT_NEQ
