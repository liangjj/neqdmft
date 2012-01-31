  !     Variables:
  dt      = 0.157080
  beta    = 100.0
  U       = 6.0
  Efield  = 0.0
  Vpd     = 0.0
  ts      = 1.0
  nstep   = 50
  nloop   = 30
  eps_error= 1.d-4
  Nsuccess = 2

  !     Efield:
  Ex  = 1.d0
  Ey  = 1.d0
  t0  = 0.d0
  t1  = dt*dble(nstep+10)
  tau0= 1.d0
  w0  = 20.d0
  field_profile='constant'

  !     Flags:
  irdeq         = .false.
  method        = 'ipt'
  update_wfftw  = .false.
  solve_wfftw   = .false.
  plotVF        = .false.
  plot3D        = .false.
  fchi          = .false.
  equench       = .false.

  !     Parameters:
  L          = 1024
  Ltau       = 32
  Lmu        = 2048
  Lkreduced  = 200
  wbath      = 10.0
  bath_type  = "constant"
  eps        = 0.05d0
  irdG0file  = "eqG0w.restart"
  irdnkfile  = "eqnk.restart"
  omp_num_threads =1      

  !     LatticeN
  Nx = 50
  Ny = 50    

  !     Quench
  iquench = .false.
  beta0   = 100.0
  U0      = 6.0
  xmu0    = 0.0    
