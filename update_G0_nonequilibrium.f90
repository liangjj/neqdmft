  !=======Component by component inversion==========================
  if(TT)then
     forall(i=0:nstep,j=0:nstep)
        locGret(i,j)= heaviside(t(i)-t(j))*(locG%gtr(i,j) - locG%less(i,j))
        Sret(i,j)   = heaviside(t(i)-t(j))*(Sig%gtr(i,j) - Sig%less(i,j))
     end forall
     locGadv=conjg(transpose(locGret))
     Sadv=conjg(transpose(Sret))

     ! !G0ret = [\11 + Gret * Sret]^-1 * Gret
     Uno=zero  ; forall(i=0:nstep)Uno(i,i)=One/dt
     GammaRet(0:nstep,0:nstep) = Uno+matmul(locGret(0:nstep,0:nstep),Sret(0:nstep,0:nstep))*dt
     GammaRet(0:nstep,0:nstep) = GammaRet(0:nstep,0:nstep)*dt**2
     call mat_inversion(GammaRet(0:nstep,0:nstep))
     G0ret(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),locGret(0:nstep,0:nstep))*dt
     forall(i=0:nstep)G0ret(i,i)=-xi !???
     G0adv=conjg(transpose(G0ret))

     !G0less = GammaR^-1 * Gless * GammaA^-1  -  gR * Sless * gA
     G0%less(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),matmul(locG%less(0:nstep,0:nstep),&
          conjg(transpose(GammaRet(0:nstep,0:nstep))))*dt)*dt -&
          matmul(G0ret(0:nstep,0:nstep),matmul(Sig%less(0:nstep,0:nstep),G0adv(0:nstep,0:nstep))*dt)*dt

     !G0gtr  = GammaR^-1 * Ggtr * GammaA^-1   -  gR * Sgtr * gA
     G0%gtr(0:nstep,0:nstep)  = matmul(GammaRet(0:nstep,0:nstep),matmul(locG%gtr(0:nstep,0:nstep),&
          conjg(transpose(GammaRet(0:nstep,0:nstep))))*dt)*dt  -&
          matmul(G0ret(0:nstep,0:nstep),matmul(Sig%gtr(0:nstep,0:nstep),G0adv(0:nstep,0:nstep))*dt)*dt
  endif





