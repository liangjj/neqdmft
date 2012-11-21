program getDATA
  USE VARS_GLOBAL
  USE ELECTRIC_FIELD
  USE BATH
  implicit none
  !Add here the extra variables 
  integer                                :: i,j,ik,loop,narg,iarg,k,ia,ir,irel,iave
  character(len=16)                      :: DIR
  complex(8),dimension(:,:),allocatable  :: locGret,Sret,impGret,G0ret
  !type(keldysh_contour_gf)               :: guessG0
  logical                                :: file
  !Wigner variables:
  real(8),dimension(:),allocatable       :: trel,tave
  complex(8),dimension(:,:),allocatable  :: wgnGgtr,wgnGless,wgnGret
  complex(8),dimension(:,:),allocatable  :: wgnSgtr,wgnSless,wgnSret
  complex(8),dimension(:,:),allocatable  :: gfret_wgn,sfret_wgn,gfless_wgn
  real(8),dimension(:,:),allocatable     :: nf_wgn

  call read_input_init("used.inputFILE.in")
  include "grid_setup.f90"  
  Lk   = square_lattice_dimension(Nx,Ny)
  allocate(epsik(Lk),wt(Lk))
  wt   = square_lattice_structure(Lk,Nx,Ny)
  epsik= square_lattice_dispersion_array(Lk,ts=ts)

  call set_efield_vector()

  DIR="RESULTS"
  if(Lkreduced<200)Lkreduced=200

  call get_data()

contains

  !+-------------------------------------------------------------------+
  subroutine get_data()
    logical :: control

    call create_data_dir(DIR)
    call msg("Results are in "//trim(adjustl(trim(DIR))),id=0)   

    call massive_allocation()

    call get_thermostat_bath()

    !Read the functions:
    call read_keldysh_contour_gf(G0,trim(data_dir)//"/G0")
    call read_keldysh_contour_gf(locG,trim(data_dir)//"/locG")
    call read_keldysh_contour_gf(Sigma,trim(data_dir)//"/Sigma")
    !call read_keldysh_contour_gf(guessG0,trim(data_dir)//"/guessG0")

    call sread(trim(data_dir)//"/nk.data",nk(0:nstep,1:Lk))

    if(fchi)then
       call sread("locChi_11.data",chi(1,1,0:nstep,0:nstep))
       call sread("locChi_12.data",chi(1,2,0:nstep,0:nstep))
       call sread("locChi_21.data",chi(2,1,0:nstep,0:nstep))
       call sread("locChi_22.data",chi(2,2,0:nstep,0:nstep))
    endif

    forall(i=0:nstep,j=0:nstep)
       locGret(i,j)   = heaviside(t(i)-t(j))*(locG%gtr(i,j) - locG%less(i,j))
       G0ret(i,j)     = heaviside(t(i)-t(j))*(G0%gtr(i,j)   - G0%less(i,j))
    end forall

    call get_plot_quasiWigner_functions( trim(adjustl(trim(DIR))) )

    call evaluate_print_observables( trim(adjustl(trim(DIR))) )

    !call plot_wigner_functions( trim(adjustl(trim(DIR))) )

    call massive_deallocation
    return
  end subroutine get_data
  !+-------------------------------------------------------------------+


  !+-------------------------------------------------------------------+
  subroutine massive_allocation
    call global_memory_allocation() !allocate memory as in the main code
    !+ something extra:
    !call allocate_keldysh_contour_gf(guessG0,nstep)
    allocate(locGret(0:nstep,0:nstep))
    allocate(Sret(0:nstep,0:nstep))
    allocate(G0ret(0:nstep,0:nstep))

    allocate(trel(-nstep:nstep),tave(0:nstep))
    allocate(wgnGless(-nstep:nstep,0:nstep),&
         wgnGgtr(-nstep:nstep,0:nstep),     &
         wgnGret(-nstep:nstep,0:nstep))
    allocate(wgnSless(-nstep:nstep,0:nstep),&
         wgnSgtr(-nstep:nstep,0:nstep),     &
         wgnSret(-nstep:nstep,0:nstep))
    allocate(gfret_wgn(0:nstep,-nstep:nstep),gfless_wgn(0:nstep,-nstep:nstep),sfret_wgn(0:nstep,-nstep:nstep))
    allocate(nf_wgn(0:nstep,-nstep:nstep))
    if(fchi)then
       allocate(chi(2,2,0:nstep,0:nstep))
    endif
  end subroutine massive_allocation
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  subroutine massive_deallocation   
    deallocate(locGret,Sret,G0ret)
    deallocate(nk)
    call deallocate_gf(gf0)
    call deallocate_gf(gf)
    call deallocate_gf(sf)
    deallocate(trel,tave)
    deallocate(wgnGless,wgnGgtr,wgnSless,wgnSgtr,wgnGret,wgnSret)
    deallocate(gfret_wgn,gfless_wgn,sfret_wgn)
    deallocate(nf_wgn)
    if(fchi)then
       deallocate(chi)
    endif
  end subroutine massive_deallocation
  !+-------------------------------------------------------------------+


  !+-------------------------------------------------------------------+
  subroutine evaluate_print_observables(dir)
    character(len=*)                          :: dir
    integer                                   :: i,ik,ix,iy,it,is,step
    complex(8)                                :: I1,Ib
    real(8)                                   :: Wtot
    real(8),dimension(Nstep) :: intJ
    type(vect2D)                              :: Ak,kt,Jk
    type(vect2D),dimension(0:nstep)           :: Jloc,Jheat,Xpos !local Current 
    type(vect2D),dimension(0:nstep,0:Nx,0:Ny) :: Jkvec,Tloc                  !current vector field
    real(8),dimension(0:nstep,Lk)             :: npi                    !covariant occupation n(\pi=\ka+\Ekt) 
    real(8),dimension(0:nstep)                :: nt,Jint,Stot   !occupation(time)
    real(8),dimension(0:Nx,0:Ny,0:nstep)      :: nDens                  !occupation distribution on the k-grid
    real(8),dimension(0:nstep)                :: Ekin,Epot,Eb,Etot,doble!energies and double occ.
    real(8),allocatable,dimension(:)          :: sorted_epsik
    integer,allocatable,dimension(:)          :: sorted_ik
    real(8),dimension(0:nstep,Lk)             :: sorted_nk,sorted_npi   !sorted arrays
    real(8),dimension(0:nstep,Lk)             :: epi,sorted_epi   !sorted arrays
    integer,dimension(:),allocatable          :: reduced_ik             !reduced arrays
    real(8),dimension(:),allocatable          :: reduced_epsik
    real(8),dimension(:,:),allocatable        :: reduced_nk,reduced_npi,reduced_epi
    real(8),dimension(2,2,0:nstep,0:nstep)    :: scond,sscond

    call msg("Print Out Results (may take a while)")

    !SORTING:
    call msg("Sorting")
    allocate(sorted_epsik(Lk),sorted_ik(Lk))

    sorted_epsik=epsik ; call sort_array(sorted_epsik,sorted_ik)

    forall(i=0:nstep,ik=1:Lk)sorted_nk(i,ik) = nk(i,sorted_ik(ik))


    !COVARIANT transformation: \ka --> \pi = \ka + \Ek*t
    call msg("Covariant transformation")
    call shift_kpoint(nk(0:nstep,1:Lk), npi(0:nstep,1:Lk))
    call shift_kpoint(sorted_nk(0:nstep,1:Lk), sorted_npi(0:nstep,1:Lk))


    !REDUCTION of the k-grid:
    call msg("Reducing BZ")
    step=Lk/Lkreduced; if(step==0)step=1
    call square_lattice_reduxGrid_dimension(Lk,step,Lkreduced)
    allocate(reduced_ik(Lkreduced),&
         reduced_epsik(Lkreduced), &
         reduced_epi(0:nstep,Lkreduced),&
         reduced_nk(0:nstep,Lkreduced),reduced_npi(0:nstep,Lkreduced))
    call square_lattice_reduxGrid_index(Lk,step,reduced_ik)
    call square_lattice_reduxGrid_dispersion_array(sorted_epsik,reduced_ik,reduced_epsik)
    forall(ik=1:Lkreduced)reduced_nk(0:nstep,ik)  = sorted_nk(0:nstep,reduced_ik(ik))
    forall(ik=1:Lkreduced)reduced_npi(0:nstep,ik) = sorted_npi(0:nstep,reduced_ik(ik))

    !Get the CURRENT Jloc(t)=\sum_\ka J_\ka(t) = -e n_\ka(t)*v_\ka(t)
    call msg("Current Field")
    Jloc=Vzero ; Jheat=Vzero   ;Stot=0.d0
    do ik=1,Lk
       ix=ik2ix(ik);iy=ik2iy(ik)
       do i=0,nstep
          Ak          = Afield(t(i),Ek)
          kt          = kgrid(ix,iy)-Ak
          epi(i,ik)   = square_lattice_dispersion(kt)
          Jk          = nk(i,ik)*square_lattice_velocity(kt)
          Jkvec(i,ix,iy)  = Jk
          Jloc(i)         = Jloc(i) +  wt(ik)*Jk !todo *2.d0 spin degeneracy
          Jheat(i)        = Jheat(i)+  wt(ik)*epi(i,ik)*Jk !todo *2.d0 spin degeneracy
          if(i>0)Stot(i)  = Stot(i) -  wt(ik)*(nk(i,ik)*log(nk(i,ik))+(1.d0-nk(i,ik))*log(1.d0-nk(i,ik)))
       enddo
    enddo
    Stot(0)=0.d0

    forall(i=0:nstep,ik=1:Lk)sorted_epi(i,ik) = epi(i,sorted_ik(ik))
    forall(ik=1:Lkreduced)reduced_epi(0:nstep,ik) = sorted_epi(0:nstep,reduced_ik(ik))

    !OCCUPATION :
    call msg("Get occupation")
    forall(i=0:nstep)nt(i)=-xi*locG%less(i,i)

    !OCCUPATION density:
    forall(i=0:nstep,ik=1:Lk)nDens(ik2ix(ik),ik2iy(ik),i)=npi(i,ik)

    !ENERGY kinetic
    Ekin=zero
    do ik=1,Lk
       ix=ik2ix(ik);iy=ik2iy(ik)
       do i=0,nstep
          Ak=Afield(t(i),Ek)
          Ekin(i) = Ekin(i) +  wt(ik)*square_lattice_dispersion(kgrid(ix,iy) - Ak)*nk(i,ik)
       enddo
    enddo


    !ENERGY potential && total 
    !Get Epot = <V>(it)= xi/2\lim_{t'-->t}\sum_k[xi\partial_t - h_{0,k}(t)] G^<_k(t,t') 
    ! by eq. of motion = xi/2\lim_{t'-->t}\sum_k{\delta^<(t,t')+\int_0^t S^R(t,z)*G^<_k(z,t')+\int_0^t' S^<(t,z)*G^A_k(z,t')}
    !                  = xi/2\lim_{t'-->t}{\delta^<(t,t')+\int_0^t S^R(t,z)*G^<_loc(z,t')+\int_0^t' S^<(t,z)*G^A_loc(z,t')}
    !                  = xi/2{1+\int_0^t [S^R(t,z)*G^<_loc(z,t) + S^<(t,z)*G^A_loc(z,t)]}
    ! cnst disregarded = xi/2\int_0^t {S^R(t,z)*G^<_loc(z,t) + S^<(t,z)*[G^R_loc(t,z)]^+}
    do it=0,nstep
       I1=zero; Ib=zero
       do i=1,it
          I1 = I1 + SretF(it,i)*locG%less(i,it) + Sigma%less(it,i)*conjg(GretF(it,i))
          Ib = Ib + S0retF(it,i)*locG%less(i,it) + S0%less(it,i)*conjg(GretF(it,i))
       enddo
       Epot(it)= -xi*I1*dt/2.d0
       Eb(it)= -xi*Ib*dt/2.d0
    enddo
    !Get Etot = Ekin + Epot
    Etot=Ekin+Epot

    Jint=modulo(Jloc)!Jloc%x + Jloc%y
    Wtot=sum(Jint(0:))*Efield*dt


    do i=0,nstep
       Xpos(i)=0.d0
       do ik=1,Lk
          ix=ik2ix(ik);iy=ik2iy(ik)          
          do k=0,i
             Jk=Jkvec(k,ix,iy)
             Xpos(i)=Xpos(i)+Jk*wt(ik)*dt
          enddo
       enddo
    enddo

    do i=1,Nstep
       intJ(i)=sum(Jint(0:i))*dt/t(i)
    enddo

    !Double OCCUPATION:
    doble= 0.5d0*(2.d0*nt) - 0.25d0 ; if(U/=0)doble = Epot/U + 0.5d0*(2.d0*nt)- 0.25d0

    if(fchi)then
       scond=0.d0
       do i=0,nstep
          do j=0,nstep
             do ik=j,nstep
                scond(1,1,i,j)=scond(1,1,i,j)+chi(1,1,i,ik)*dt
                scond(1,2,i,j)=scond(1,2,i,j)+chi(1,2,i,ik)*dt
                scond(2,1,i,j)=scond(2,1,i,j)+chi(2,1,i,ik)*dt
                scond(2,2,i,j)=scond(2,2,i,j)+chi(2,2,i,ik)*dt
             enddo
          enddo
       enddo

       sscond=zero
       do it=0,nstep
          do is=0,it
             sscond(:,:,it,it-is)=scond(:,:,it,is)
          enddo
       enddo

    endif

    !====================================================================================
    !PRINT:
    call msg("Print n(t)")
    call splot(dir//"/nVStime.ipt",t(0:nstep),2.d0*nt(0:nstep))


    call msg("Print J(t)")
    call splot(dir//"/JlocVStime.ipt",t(0:nstep),Jloc(0:nstep)%x,Jloc(0:nstep)%y,Jint(0:nstep))
    call splot(dir//"/intJVStime.ipt",t(1:nstep),intJ(1:nstep))
    call splot(dir//"/JheatVStime.ipt",t(0:nstep),Jheat(0:nstep)%x,Jheat(0:nstep)%y)
    call splot(dir//"/RposVStime.ipt",t(0:nstep),Xpos(0:nstep)%x,Xpos(0:nstep)%y,modulo(Xpos))

    call msg("Print Ex(t)")
    call splot(dir//"/EkinVStime.ipt",t(0:nstep),Ekin(0:nstep))
    call splot(dir//"/EpotVStime.ipt",t(0:nstep),Epot(0:nstep))
    call splot(dir//"/EhybVStime.ipt",t(0:nstep),Eb(0:nstep))
    call splot(dir//"/EtotVStime.ipt",t(0:nstep),Etot(0:nstep),Etot(0:nstep)+Eb(0:nstep))
    call splot(dir//"/WtotVSefield.ipt",Efield,Wtot)

    call msg("Print S(t)")
    call splot(dir//"/StotVStime.ipt",t(0:nstep),Stot(0:nstep))

    call msg("Print d(t)")
    call splot(dir//"/doccVStime.ipt",t(0:nstep),doble(0:nstep))


    !DISTRIBUTION:
    call msg("Print n(k,t)")
    call splot("nVStimeVSepsk3D.ipt",t(0:nstep),reduced_epsik,reduced_nk(0:nstep,:))
    do i=0,nstep
       call splot("nVSepi.ipt",reduced_epsik(:),reduced_npi(i,:),append=TT)
    enddo

    !Fermi Surface plot:
    call msg("Print FS(k,t)")
    call splot("3dFSVSpiVSt",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,nDens(0:Nx,0:Ny,0:Nstep))

    !Current Vector Field:
    !if(Efield/=0.d0 .AND. plotVF)call dplot_vector_field("vf_JfieldVSkVSt",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,Jkvec(0:nstep,:,:)%x,Jkvec(0:nstep,:,:)%y)

    if(fchi)then
       call splot("sigma_cond.ipt",t(0:nstep),t(0:nstep),scond(1,1,0:nstep,0:nstep))
       call splot("red_sigma_cond.ipt",t(0:nstep),t(0:nstep),sscond(1,1,0:nstep,0:nstep))
    endif

    !Local functions:
    !===========================================================================
    if(.not.plot3D)then         !functions have not be plotted during run, plot them now
       !if(Efield/=0.d0 .or. Vpd/=0.0)call dplot_3d_surface_animated("3dFSVSpiVSt","$k_x$","$k_y$","$FS(k_x,k_y)$",&
       !kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,nDens(0:Nx,0:Ny,0:nstep))
       !call plot_keldysh_contour_gf(guessG0,t(0:),trim(plot_dir)//"/guessG0")
       call plot_keldysh_contour_gf(G0,t(0:),trim(plot_dir)//"/G0")
       call plot_keldysh_contour_gf(locG,t(0:),trim(plot_dir)//"/locG")
       call plot_keldysh_contour_gf(Sigma,t(0:),trim(plot_dir)//"/Sigma")
    endif
    call system("rm -rf "//dir//"/3d* 2>/dev/null")
    call system("rm -rf "//dir//"/*3D 2>/dev/null")
    call system("rm -rf "//dir//"/gif_* 2>/dev/null")
    call system("mv -vf  gif_* vf_* 3d* *3D* nVSepi.ipt *cond* "//dir//"/ 2>/dev/null")
    return
  end subroutine evaluate_print_observables
  !+-------------------------------------------------------------------+



  !+-------------------------------------------------------------------+
  subroutine get_plot_quasiWigner_functions(dir)
    character(len=*)  :: dir
    integer           :: i,j
    real(8),dimension(2*nstep)         :: phi
    complex(8),dimension(-nstep:nstep) :: gtkel
    complex(8),dimension(2*nstep)      :: gfkel

    forall(i=0:nstep,j=0:nstep)
       gf%less%t(i-j) = locG%less(i,j)
       gf%gtr%t(i-j)  = locG%gtr(i,j)
       gf%ret%t(i-j)  = locGret(i,j)
    end forall
    if(heaviside(0.d0)==1.d0)gf%ret%t(0)=gf%ret%t(0)/2.d0
    call fftgf_rt2rw(gf%ret%t,gf%ret%w,nstep) ;    gf%ret%w=gf%ret%w*dt ; call swap_fftrt2rw(gf%ret%w)
    call splot(dir//"/locGless_t.ipt",t(-nstep:nstep),gf%less%t,append=TT)
    call splot(dir//"/locGgtr_t.ipt",t(-nstep:nstep),gf%gtr%t,append=TT)
    call splot(dir//"/locGret_t.ipt",t(-nstep:nstep),gf%ret%t,append=TT)
    call splot(dir//"/locGret_realw.ipt",wr,gf%ret%w,append=TT)
    call splot(dir//"/locDOS.ipt",wr,-aimag(gf%ret%w)/pi)

    forall(i=0:nstep,j=0:nstep)
       gf0%less%t(i-j)= G0%less(i,j)
       gf0%gtr%t(i-j) = G0%gtr(i,j)
       gf0%ret%t(i-j) = G0ret(i,j)
       sf%less%t(i-j) = Sigma%less(i,j)
       sf%gtr%t(i-j)  = Sigma%gtr(i,j)
       sf%ret%t(i-j)  = Sret(i,j)
    end forall
    if(heaviside(0.d0)==1.d0)gf0%ret%t(0)=gf0%ret%t(0)/2.d0 !; gf0%ret%t(0)=-xi
    if(heaviside(0.d0)==1.d0)sf%ret%t(0)=sf%ret%t(0)/2.d0
    if(heaviside(0.d0)==1.d0)gf%ret%t(0)=gf%ret%t(0)/2.d0   !; gf%ret%t(0)=-xi
    ! if(loop==1)then
    !    call splot(dir//"/guessG0less_t.ipt",t(-nstep:nstep),gf0%less%t)
    !    call splot(dir//"/guessG0gtr_t.ipt",t(-nstep:nstep),gf0%gtr%t)
    !    call splot(dir//"/guessG0ret_t.ipt",t(-nstep:nstep),gf0%ret%t)
    ! else
    call splot(dir//"/G0less_t.ipt",t(-nstep:nstep),gf0%less%t)
    call splot(dir//"/G0gtr_t.ipt",t(-nstep:nstep),gf0%gtr%t)
    call splot(dir//"/G0ret_t.ipt",t(-nstep:nstep),gf0%ret%t)
    ! endif
    call splot(dir//"/Sless_t.ipt",t(-nstep:nstep),sf%less%t)
    call splot(dir//"/Sgtr_t.ipt",t(-nstep:nstep),sf%gtr%t)
    call splot(dir//"/Sret_t.ipt",t(-nstep:nstep),sf%ret%t)


    !Obtain && plot Real frequency Functions:
    !===========================================================================
    call fftgf_rt2rw(gf0%ret%t,gf0%ret%t,nstep) ;    gf0%ret%t=gf0%ret%t*dt ; call swap_fftrt2rw(gf0%ret%t)
    call fftgf_rt2rw(gf%ret%t,gf%ret%t,nstep)   ;    gf%ret%t=gf%ret%t*dt   ; call swap_fftrt2rw(gf%ret%t)
    call fftgf_rt2rw(sf%ret%t,sf%ret%t,nstep)   ;    sf%ret%t=dt*sf%ret%t   ; call swap_fftrt2rw(sf%ret%t)
    ! if(loop==1)then
    !    call splot(dir//"/guessG0ret_realw.ipt",wr,gf0%ret%t)
    ! else
    call splot(dir//"/G0ret_realw.ipt",wr,gf0%ret%t)
    ! endif
    call splot(dir//"/Sret_realw.ipt",wr,sf%ret%t)
    call splot(dir//"/DOS.ipt",wr,-aimag(gf%ret%t)/pi)

    forall(i=0:nstep,j=0:nstep)gtkel(i-j) = locG%less(i,j)+locG%gtr(i,j)
    call fftgf_rt2rw(gtkel,gfkel,nstep) ; gfkel=gfkel*dt ; call swap_fftrt2rw(gfkel)
    call splot(dir//"/locGkel_realw.ipt",wr,gfkel)
    phi = xi*gfkel/aimag(gf%ret%w)/2.d0
    do i=1,2*nstep
       if(wr(i)>-1.d0)exit
    enddo
    call splot(dir//"/phi_realw.ipt",wr(i:(2*nstep-i)),phi(i:(2*nstep-i)))
    return
  end subroutine get_plot_quasiWigner_functions
  !+-------------------------------------------------------------------+



  !+-------------------------------------------------------------------+
  subroutine plot_wigner_functions(dir)
    character(len=*)  :: dir
    complex(8)        :: delta

    call init_trel(t,trel,nstep)
    call init_tave(t,tave,nstep)

    !Perform the Wigner Rotation:
    call msg("Perform Wigner Rotation")
    wgnGless= wigner_transform(locG%less,nstep)
    wgnGret = wigner_transform(locGret,nstep)
    wgnSless= wigner_transform(Sigma%less,nstep)
    wgnSret = wigner_transform(Sret,nstep)
    call system("if [ ! -d WIGNER ]; then mkdir WIGNER; fi")
    ! call plot_3D("wgnGless3D","X","Y","Z",trel(-nstep:nstep)/dt,tave(1:nstep)/dt,wgnGless(-nstep:nstep,1:nstep))
    ! call plot_3D("wgnSless3D","X","Y","Z",trel(-nstep:nstep)/dt,tave(1:nstep)/dt,wgnSless(-nstep:nstep,1:nstep))
    call system("mv wgn*3D WIGNER/")

    delta=(one+xi)/dble(nstep)
    do ia=0,nstep
       call fftgf_rt2rw(wgnGret(:,ia),gfret_wgn(ia,:),nstep)  ;gfret_wgn(ia,:)=gfret_wgn(ia,:)*dt;call swap_fftrt2rw(gfret_wgn(ia,:))
       call fftgf_rt2rw(wgnGless(:,ia),gfless_wgn(ia,:),nstep);gfless_wgn(ia,:)=gfless_wgn(ia,:)*dt;call swap_fftrt2rw(gfless_wgn(ia,:))
       call fftgf_rt2rw(wgnSret(:,ia),sfret_wgn(ia,:),nstep)  ;sfret_wgn(ia,:)=sfret_wgn(ia,:)*dt;call swap_fftrt2rw(sfret_wgn(ia,:))
       call splot("WIGNER/wgnDOS.ipt",wr,(-aimag(gfret_wgn(ia,:) - delta*dble(ia)*pi)/pi),append=TT)
       call splot("WIGNER/wgnSigma_realw.ipt",wr,(sfret_wgn(ia,:) + delta*dble(ia)),append=TT)
       call splot("WIGNER/wgnGless_realw.ipt",wr,(gfless_wgn(ia,:) + delta*dble(ia)),append=TT)
       ! nf_wgn(ia,:) = -aimag(gfless_wgn(ia,:))!/aimag(gfret_wgn(ia,:))/pi2
       nf_wgn(ia,:) = -xi*gfless_wgn(ia,:)/aimag(gfret_wgn(ia,:))
       call splot("n_wgnVSepi.ipt",wr(:),nf_wgn(ia,:),append=TT)
    enddo
    call splot("wgndosVSrealwVStime.ipt",tave(0:nstep),wr(-nstep:nstep),-aimag(gfret_wgn(0:nstep,-nstep:nstep))/pi)
    call splot("wgnnfVSrealwVStime.ipt",tave(0:nstep),wr(-nstep:nstep),nf_wgn(0:nstep,-nstep:nstep))
    call system("mv *wgn*  WIGNER/")
    !call plot_3D("wgnDOS3D","X","Y","Z",tave(0:nstep)/dt,wr(1:2*nstep),-aimag(gfret_wgn(0:nstep,1:2*nstep))/pi)

  end subroutine plot_wigner_functions
  !+-------------------------------------------------------------------+



  !+-------------------------------------------------------------------+
  subroutine shift_kpoint(arrayIN,arrayOUT)
    integer                 :: i,j,ik,ix,iy,jk,jx,jy
    real(8),dimension(0:,:) :: arrayIN,arrayOUT
    real(8),dimension(2)    :: pi_in
    integer,dimension(2)    :: pi_kcoord
    type(vect2D)            :: Ak
    arrayOUT=0.d0
    do i=0,nstep
       do ik=1,Lk
          ix=ik2ix(ik);iy=ik2iy(ik)
          !find the Xcoord of shifted point:
          Ak=Afield(t(i),Ek)
          pi_in(1)=kgrid(ix,iy)%x + Ak%x !+ (-t(i))*Ek%x
          do j=1,1000000
             if(pi_in(1) > pi) then
                pi_in(1)=pi_in(1) - pi2
             elseif(pi_in(1) < -pi) then
                pi_in(1)=pi_in(1) + pi2
             else
                exit
             endif
          enddo
          !find the Ycoord of shifted point:
          pi_in(2)=kgrid(ix,iy)%y + Ak%y !+ (-t(i))*Ek%y
          do j=1,1000000
             if(pi_in(2) > pi) then
                pi_in(2)=pi_in(2) - pi2
             elseif(pi_in(2) < -pi) then
                pi_in(2)=pi_in(2) + pi2
             else
                exit
             endif
          enddo
          !FIND the kgrid point corresponding to Pi-point
          call find2Dmesh(kgrid(0:Nx,1)%x,kgrid(1,0:Ny)%y,pi_in,pi_kcoord)
          jx=pi_kcoord(1)-1 ; jy=pi_kcoord(2)-1
          if(jx < 0  .or. jx > Nx)print*,"error jx=",jx
          if(jy < 0  .or. jy > Ny)print*,"error jy=",jy
          jk=kindex(jx,jy)
          arrayOUT(i,ik)=arrayIN(i,jk)
       enddo
    enddo
  end subroutine shift_kpoint
  !+-------------------------------------------------------------------+


  !+-------------------------------------------------------------------+
  function wigner_transform(Gin,M)
    integer                       :: M,ir,ia,i,j
    complex(8),dimension(0:M,0:M) :: Gin
    complex(8),dimension(-M:M,0:M):: Gout,wigner_transform
    do ir=-nstep,nstep
       do ia=0,nstep
          forall(j=0:nstep,i=0:nstep, i-j==ir .AND. (i+j)/2==ia)
             Gout(ir,ia)=Gin(i,j)
          end forall
       enddo
    enddo
    wigner_transform=Gout
  end function wigner_transform
  !+-------------------------------------------------------------------+


  !+-------------------------------------------------------------------+
  subroutine init_trel(time,tr,M)
    integer                            :: i,j,ir
    integer,intent(in)                 :: M
    real(8),dimension(-M:M),intent(in) :: time
    real(8),dimension(-M:M),intent(out):: tr
    do ir=-nstep,nstep
       forall(j=0:nstep,i=0:nstep, i-j==ir)tr(ir)=time(i)-time(j)
    enddo
    return
  end subroutine init_trel
  !+-------------------------------------------------------------------+


  !+-------------------------------------------------------------------+
  subroutine init_tave(time,ta,M)
    integer                            :: i,j,ia
    integer,intent(in)                 :: M
    real(8),dimension(-M:M),intent(in) :: time
    real(8),dimension(0:M),intent(out) :: ta
    do ia=0,nstep
       forall(j=0:nstep,i=0:nstep, (i+j)/2==ia)ta(ia)=(time(i)+time(j))/2.d0
    enddo
    return
  end subroutine init_tave
  !+-------------------------------------------------------------------+


  !+-------------------------------------------------------------------+
  function GretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: GretF
    GretF = heaviside(t(i)-t(j))*(locG%gtr(i,j)-locG%less(i,j))
  end function GretF
  !-------------------------------------------------------!


  !-------------------------------------------------------!
  function SretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: SretF
    SretF = heaviside(t(i)-t(j))*(Sigma%gtr(i,j)-Sigma%less(i,j))
  end function SretF
  !-------------------------------------------------------!


  !-------------------------------------------------------!
  function S0retF(i,j)
    integer,intent(in) :: i,j
    complex(8)         :: S0retF
    S0retF = heaviside(t(i)-t(j))*(S0%gtr(i,j)-S0%less(i,j))
  end function S0retF
  !-------------------------------------------------------!

end program getDATA
