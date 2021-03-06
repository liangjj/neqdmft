  allocate(wr(2*nstep),t(-nstep:nstep))
  allocate(wm(L),tau(0:Ltau))
  wmax  = pi/dt
  t     = linspace(-dt*real(nstep,8),dt*real(nstep,8),2*nstep+1)
  wr    = linspace(-wmax,wmax,2*nstep,mesh=fmesh)-fmesh/2.d0
  !t     = wr(-nstep:nstep)/fmesh*dt
  wm    = pi/beta*real(2*arange(1,L)-1,8)
  tau   = linspace(0.d0,beta,Ltau+1,mesh=dtau)
  if(mpiID==0)then
     write(*,'(A,F12.6)')"dt   =",dt
     write(*,'(A,F12.6)')"dw   =",fmesh
     write(*,'(A,F12.6)')"wmax =",wmax
     write(*,'(A,F12.6)')"tmax =",dt*dble(nstep)
  endif
