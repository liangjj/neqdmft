  allocate(wr(2*nstep),t(-nstep:nstep))
  wmax = pi/dt
  tmax = dt*real(nstep,8)
  t    = linspace(-tmax,tmax,2*nstep+1)
  wr   = linspace(-wmax,wmax,2*nstep,mesh=fmesh)-fmesh/2.d0
  !t   = wr(-nstep:nstep)/fmesh*dt

  allocate(wm(L))
  wm   = pi/beta*real(2*arange(1,L)-1,8)

  allocate(tau(-Ltau:Ltau),ftau(-Ltau:Ltau))
  tau(0:)  = linspace(0.d0,beta,Ltau+1,mesh=dtau)
  ftau(0:) = upminterval(0.d0,beta,beta/2.d0,P,Q,type=0)
  forall(i=1:Ltau)tau(-i)  =-tau(i)
  forall(i=1:Ltau)ftau(-i) =-ftau(i)

  if(mpiID==0)then
     call msg("dt   ="//txtfy(dt))
     call msg("dw   ="//txtfy(fmesh))
     call msg("wmax ="//txtfy(wmax))
     call msg("tmax ="//txtfy(dt*dble(nstep)))
  endif
