  forall(i=0:nstep,j=0:nstep)
     gf0%ret%t(i-j)=heaviside(t(i-j))*(G0gtr(i,j) - G0less(i,j))
     sf%ret%t(i-j)=heaviside(t(i-j))*(Sgtr(i,j) - Sless(i,j))
  end forall
  if(heaviside(0.d0)==1.d0)gf0%ret%t(0)=gf0%ret%t(0)/2.d0
  if(heaviside(0.d0)==1.d0)sf%ret%t(0)=sf%ret%t(0)/2.d0
  call fftgf_rt2rw(gf0%ret%t,gf0%ret%w,nstep) ;  gf0%ret%w=gf0%ret%w*dt ; call swap_fftrt2rw(gf0%ret%w) !swap because F(t) are not oscillating in this formalism:
  call fftgf_rt2rw(sf%ret%t,sf%ret%w,nstep)  ;  sf%ret%w=dt*sf%ret%w   ; call swap_fftrt2rw(sf%ret%w)   !swap because F(t) are not oscillating in this formalism:
  gf%ret%w = one/(one/gf0%ret%w - sf%ret%w)
  do i=1,2*nstep
     w = wr(i)
     A=-aimag(gf%ret%w(i))/pi
     gf%less%w(i)= pi2*xi*fermi(w,beta)*A
     gf%gtr%w(i) = pi2*xi*(fermi(w,beta)-1.d0)*A
  enddo
  call fftgf_rw2rt(gf%less%w,gf%less%t,nstep)  ; gf%less%t=fmesh/pi2*gf%less%t ;  gf%less%t=gf%less%t*exa 
  call fftgf_rw2rt(gf%gtr%w,gf%gtr%t,nstep)   ; gf%gtr%t =fmesh/pi2*gf%gtr%t  ;  gf%gtr%t=gf%gtr%t*exa
  forall(i=0:nstep,j=0:nstep)
     impGless(i,j)= gf%less%t(i-j)
     impGgtr(i,j) = gf%gtr%t(i-j)
  end forall
