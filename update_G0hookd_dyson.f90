  !\Gamma(tau-tau') = \delta(tau-tau') + \int_0^beta dz S(tau-z)*G(z-tau')
  !u=z-tau'; du = dz; z=u+tau'
  !\Gamam(tau-tau') = \delta(tau-tau') + \int_0^beta du S(tau-tau'-u)*G(u)
  !tau-tau'=t
  !\Gamma(t) = \delta(t) + \int_0^beta du S(t-u)*G(u)
  !FFT:
  !\Gamma(w) = one + S(w)*G(w) ; \Gamma(w)^{-1} = one/Gamma(w)
  
  
  forall(i=0:Ltau,j=0:Ltau) 
     locGtau(i,j) = icGtau(i-j)
     Stau(i,j)    = icStau(i-j)
  end forall

  ! do i=0,Ltau
  !    do j=0,Ltau
  !       R=i-j ; deg=dble(L+1-abs(R))
  !       do k=0,Ltau
  !          GammatM(i-j)=GammatM(i-j) + Stau(i,k)*locGtau(k,j)*dtau/deg
  !       enddo
  !    enddo
  ! enddo
  ! GammatM(0)=one/dtau+GammatM(0)

  GammafM = One + icSiw*icGiw
  call splot("dysonGammaiw.ipt",wm,GammafM)
  GammafM=one/GammafM 
  call fftgf_iw2tau(GammafM,GammatM_,beta) 
  call extract(GammatM_,GammatM(0:Ltau))
  do i=0,L
     write(90,*)dble(i)*beta/dble(L),GammatM_(i)
  enddo
  forall(i=0:Ltau) GammatM(-i) = -GammatM(Ltau-i)
  call splot("dysonGammatM.ipt",tau,GammatM(0:Ltau))
  stop

  forall(i=0:Ltau,j=0:Ltau) GammaM(i,j)  = GammatM(i-j)

  !G0lceil  = GammaR^-1 * Glceil * GammaM^-1   -   gR * Slceil * gM
  Op1lceil=zero
  Op2lceil=zero
  do i=0,nstep
     do itau=0,Ltau          
        do jtau=0,Ltau
           Op1lceil(i,itau)=Op1lceil(i,itau) + locGlceil(i,jtau)*GammatM(jtau-itau)*dtau
           Op2lceil(i,itau)=Op2lceil(i,itau) + Slceil(i,jtau)*icGtau(jtau-itau)*dtau
        enddo
     enddo
  enddo

  G0lceil(0:nstep,0:Ltau) = matmul(GammaRet(0:nstep,0:nstep),Op1lceil(0:nstep,0:Ltau))*dt -&
       matmul(G0ret(0:nstep,0:nstep),Op2lceil(0:nstep,0:Ltau))*dt
  forall(i=0:Ltau)G0rceil(i,:)=conjg(G0lceil(:,Ltau-i))

  !G0less = G0less - gR * Slceil * grceil - GammaR^-1 * Glceil * (Srceil * gA + Sm * grceil)
  !G0gtr  = G0gtr  - gR * Slceil * grceil - GammaR^-1 * Glceil * (Srceil * gA + Sm * grceil)
  Op1less(0:Ltau,0:nstep) = matmul(Srceil(0:Ltau,0:nstep),&
       conjg(transpose(G0ret(0:nstep,0:nstep))))*dt
  Op2less=zero
  do itau=0,Ltau
     do i=0,nstep
        do jtau=0,Ltau
           Op2less(itau,i) = Op2less(itau,i) + icStau(itau-jtau)*G0rceil(jtau,i)*dtau
        enddo
     enddo
  enddo
  Op3less(0:nstep,0:Ltau) = matmul(GammaRet(0:nstep,0:nstep),locGlceil(0:nstep,0:Ltau))*dt
  Op1less=Op1less+Op2less

  G0less(0:nstep,0:nstep) = G0less(0:nstep,0:nstep) + matmul(G0ret(0:nstep,0:nstep),matmul(Slceil(0:nstep,0:Ltau),G0rceil(0:Ltau,0:nstep)))*dt*dtau-&
       matmul(Op3less(0:nstep,0:Ltau),Op1less(0:Ltau,0:nstep))*dtau
  G0gtr(0:nstep,0:nstep) = G0gtr(0:nstep,0:nstep) + matmul(G0ret(0:nstep,0:nstep),matmul(Slceil(0:nstep,0:Ltau),G0rceil(0:Ltau,0:nstep)))*dt*dtau-&
       matmul(Op3less(0:nstep,0:Ltau),Op1less(0:Ltau,0:nstep))*dtau

  ! call plot_3D("dysonGtau3D","$\Delta t","$\tau","Z",tau,tau,locGtau(0:Ltau,0:Ltau))
  ! call plot_3D("dysonStau3D","$\Delta t","$\tau","Z",tau,tau,Stau(0:Ltau,0:Ltau))
  call plot_3D("dysonG0lceil3D","$\Delta t","$\tau","Z",t(0:nstep)/dt,tau,G0lceil(0:nstep,0:Ltau))
  call plot_3D("dysonG0less3D","$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,G0less(0:nstep,0:nstep))
