  Lk   = square_lattice_dimension(Nx,Ny)
  allocate(epsik(Lk),wt(Lk))
  wt   = square_lattice_structure(Lk,Nx,Ny)
  epsik= square_lattice_dispersion_array(Lk,ts=ts)
  
  ! Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk))
  ! call bethe_lattice(wt,epsik,Lk,D_=D,eps_=eps)

  call get_free_dos(epsik,wt)
