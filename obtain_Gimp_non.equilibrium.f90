  dSret=zero ; dG0ret=zero ; dGret=zero
  GammaRet=zero ; Gamma0Ret=zero
  !1 - get the Ret components of G_0 && \Sigma:
  forall(i=0:nstep,j=0:nstep)
     dG0ret(i,j)=heaviside(t(i)-t(j))*(G0gtr(i,j) - G0less(i,j))
     dSret(i,j) =heaviside(t(i)-t(j))*(Sgtr(i,j) - Sless(i,j))
  end forall
  !2 - get the  operator: \Gamma_0^R = \Id - \Sigma^R\circ G_0^R && invert it
  Uno=zero  ; forall(i=0:nstep)Uno(i,i)=One/dt
  Gamma0Ret(0:nstep,0:nstep) = Uno-matmul(dSret(0:nstep,0:nstep),dG0ret(0:nstep,0:nstep))*dt
  Gamma0Ret(0:nstep,0:nstep)=Gamma0Ret(0:nstep,0:nstep)*dt**2
  call mat_inversion_GJ(Gamma0Ret(0:nstep,0:nstep))
  !3 - get G_imp^R, G_imp^{>,<} using Dyson equations:
  dGret(0:nstep,0:nstep)    = matmul(dG0ret(0:nstep,0:nstep),Gamma0Ret(0:nstep,0:nstep))*dt 
  GammaRet(0:nstep,0:nstep) = Uno + matmul(dGret(0:nstep,0:nstep),dSret(0:nstep,0:nstep))*dt

  impGless(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),&
       matmul(G0less(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep)))))*dt**2 +&
       matmul(dGret(0:nstep,0:nstep),matmul(Sless(0:nstep,0:nstep),conjg(transpose(dGret(0:nstep,0:nstep)))))*dt**2

  impGgtr(0:nstep,0:nstep)  = matmul(GammaRet(0:nstep,0:nstep),&
       matmul(G0gtr(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep)))))*dt**2  +&
       matmul(dGret(0:nstep,0:nstep),matmul(Sgtr(0:nstep,0:nstep),conjg(transpose(dGret(0:nstep,0:nstep)))))*dt**2
