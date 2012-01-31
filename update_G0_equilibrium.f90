  forall(i=0:nstep,j=0:nstep)
     gf%ret%t(i-j) = heaviside(t(i-j))*(locGgtr(i,j)-locGless(i,j))
     sf%ret%t(i-j) = heaviside(t(i-j))*(Sgtr(i,j)-Sless(i,j))
  end forall
  if(heaviside(0.d0)==1.d0)gf%ret%t(0)=gf%ret%t(0)/2.d0
  if(heaviside(0.d0)==1.d0)sf%ret%t(0)=sf%ret%t(0)/2.d0
  call fftgf_rt2rw(gf%ret%t,gf%ret%w,nstep) ; gf%ret%w=gf%ret%w*dt ; call swap_fftrt2rw(gf%ret%w)
  call fftgf_rt2rw(gf%ret%t,gf%ret%w,nstep) ; sf%ret%w=sf%ret%w*dt ; call swap_fftrt2rw(sf%ret%w)
  gf0%ret%w=one/(one/gf%ret%w + sf%ret%w)
  do i=1,2*nstep
     w = wr(i)
     A = -aimag(gf0%ret%w(i))/pi
     An= A*fermi(w,beta)
     gf0%less%w(i)= pi2*xi*An
     gf0%gtr%w(i) = pi2*xi*(An-A)
  enddo
  call fftgf_rw2rt(gf0%less%w,gf0%less%t,nstep) ; gf0%less%t=fmesh/pi2*gf0%less%t  ; gf0%less%t=gf0%less%t*exa
  call fftgf_rw2rt(gf0%gtr%w, gf0%gtr%t,nstep)  ; gf0%gtr%t =fmesh/pi2*gf0%gtr%t   ; gf0%gtr%t=gf0%gtr%t*exa
  forall(i=0:nstep,j=0:nstep)
     G0less(i,j)= gf0%less%t(i-j)
     G0gtr(i,j) = gf0%gtr%t(i-j)
  end forall

  !PLus this:
  forall(i=0:nstep,j=0:nstep)
     G0ret(i,j)=heaviside(t(i-j))*(G0gtr(i,j) - G0less(i,j))
     gf0%ret%t(i-j)=G0ret(i,j)
  end forall
  call fftgf_rt2rw(gf0%ret%t,gf0%less%w,nstep) ; gf0%less%w=gf0%less%w*dt ; call swap_fftrt2rw(gf0%less%w)

