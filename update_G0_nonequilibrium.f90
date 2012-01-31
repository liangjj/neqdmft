  !=======Component by component inversion==========================
  if(TT)then
     forall(i=0:nstep,j=0:nstep)
        locGret(i,j)= heaviside(t(i)-t(j))*(locGgtr(i,j) - locGless(i,j))
        Sret(i,j)   = heaviside(t(i)-t(j))*(Sgtr(i,j) - Sless(i,j))
     end forall
     !forall(i=0:nstep)locGret(i,i)=-xi!locGless(i,i)

     locGadv=conjg(transpose(locGret))
     Sadv=conjg(transpose(Sret))

     ! !G0ret = [\11 + Gret * Sret]^-1 * Gret
     Uno=zero  ; forall(i=0:nstep)Uno(i,i)=One/dt
     GammaRet(0:nstep,0:nstep) = Uno+matmul(locGret(0:nstep,0:nstep),Sret(0:nstep,0:nstep))*dt
     GammaRet(0:nstep,0:nstep) = GammaRet(0:nstep,0:nstep)*dt**2
     call mat_inversion(GammaRet(0:nstep,0:nstep))
     G0ret(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),locGret(0:nstep,0:nstep))*dt
     forall(i=0:nstep)G0ret(i,i)=-xi
     G0adv=conjg(transpose(G0ret))

     !G0less = GammaR^-1 * Gless * GammaA^-1  -  gR * Sless * gA
     G0less(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),matmul(locGless(0:nstep,0:nstep),&
          conjg(transpose(GammaRet(0:nstep,0:nstep))))*dt)*dt -&
          matmul(G0ret(0:nstep,0:nstep),matmul(Sless(0:nstep,0:nstep),G0adv(0:nstep,0:nstep))*dt)*dt

     !G0gtr  = GammaR^-1 * Ggtr * GammaA^-1   -  gR * Sgtr * gA
     G0gtr(0:nstep,0:nstep)  = matmul(GammaRet(0:nstep,0:nstep),matmul(locGgtr(0:nstep,0:nstep),&
          conjg(transpose(GammaRet(0:nstep,0:nstep))))*dt)*dt  -&
          matmul(G0ret(0:nstep,0:nstep),matmul(Sgtr(0:nstep,0:nstep),G0adv(0:nstep,0:nstep))*dt)*dt
  endif




  !=======Inversion of the Keldysh-Schwinger Matrix==========================
  if(FF)then
     !1) build time/antitime-ordered GF: G^t && G^at; S^t && S^at
     forall(i=0:nstep,j=0:nstep,i>=j)
        locGtt(i,j) = locGgtr(i,j)
        locGat(i,j) = locGless(i,j)
        Stt(i,j) = Sgtr(i,j)
        Sat(i,j) = Sless(i,j)
     end forall
     forall(i=0:nstep,j=0:nstep,i<j)
        locGtt(i,j) = locGless(i,j)
        locGat(i,j) = locGgtr(i,j)
        Stt(i,j) = Sless(i,j)
        Sat(i,j) = Sgtr(i,j)
     end forall
     !call plot_3D("dGtt3D","X","Y","Z",t(0:nstep),t(0:nstep),locGtt(0:nstep,0:nstep))
     !call plot_3D("dStt3D","X","Y","Z",t(0:nstep),t(0:nstep),Stt(0:nstep,0:nstep))

     !2) Build the KS matrices: \FF = {{FF^t , FF^>}, {-FF^<, -FF^at}}
     NN=nstep+1                  !Size of the KS matrices
     allocate(locGmat(1:2*NN,1:2*NN),Smat(1:2*NN,1:2*NN))
     allocate(G0mat(1:2*NN,1:2*NN),GammaMat(1:2*NN,1:2*NN),UnoMat(1:2*NN,1:2*NN))
     forall(i=1:NN,j=1:NN)
        locGmat(i,j)       = locGtt(i-1,j-1)   !++
        locGmat(i,NN+j)    = locGgtr(i-1,j-1)  !+-
        locGmat(NN+i,j)    =-locGless(i-1,j-1) !-+
        locGmat(NN+i,NN+j) =-locGat(i-1,j-1)   !--
        !
        Smat(i,j)       = Stt(i-1,j-1)   !++
        Smat(i,NN+j)    = Sgtr(i-1,j-1)  !+-
        Smat(NN+i,j)    =-Sless(i-1,j-1) !-+
        Smat(NN+i,NN+j) =-Sat(i-1,j-1)   !--
     end forall

     !3)Begin inversion of Dyson equation: 
     !G = g + g*S*G ; G = g*[I + S*G] \== g*\G ==> g = G * {\G}^-1
     UnoMat=zero   ; forall(i=1:2*NN)UnoMat(i,i)=One/dt !Form the delta function
     GammaMat = UnoMat+matmul(Smat,locGmat)*dt          !Form \G operator
     GammaMat = GammaMat*dt**2                          !Prepare for the inversion     
     call mat_inversion_GJ(GammaMat)                       !Inversion
     G0mat    = matmul(locGmat,GammaMat)*dt             !Update G0 operators:

     !4) Extract the bigger&&lesser components
     forall(i=1:NN,j=1:NN)
        G0gtr(i-1,j-1)  = G0mat(i,j+NN)
        G0less(i-1,j-1) =-G0mat(NN+i,j)
     end forall

     forall(i=0:nstep,j=0:nstep)
        G0ret(i,j)=heaviside(t(i-j))*(G0gtr(i,j) - G0less(i,j))
        gf0%ret%t(i-j)=G0ret(i,j)
     end forall
     call fftgf_rt2rw(gf0%ret%t,gf0%ret%w,nstep) ; gf0%ret%w=gf0%ret%w*dt ; call swap_fftrt2rw(gf0%ret%w)

     deallocate(locGmat,Smat,G0mat,GammaMat,UnoMat)
  endif

