  ! !=======Component by component inversion==========================
  ! forall(i=0:nstep,j=0:nstep)
  !    locGret(i,j)= heaviside(t(i)-t(j))*(locG%gtr(i,j) - locG%less(i,j))
  !    Sret(i,j)   = heaviside(t(i)-t(j))*(Sigma%gtr(i,j) - Sigma%less(i,j))
  ! end forall
  ! !forall(i=0:nstep)locGret(i,i)=-xi!locG%less(i,i)
  ! !
  ! !G0ret = [\11 + Gret * Sret]^-1 * Gret
  ! Uno=zero  ; forall(i=0:nstep)Uno(i,i)=One/dt
  ! GammaRet(0:nstep,0:nstep) = Uno+matmul(locGret(0:nstep,0:nstep),Sret(0:nstep,0:nstep))*dt
  ! GammaRet(0:nstep,0:nstep) = GammaRet(0:nstep,0:nstep)*dt**2
  ! call mat_inversion(GammaRet(0:nstep,0:nstep))
  ! G0ret(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),locGret(0:nstep,0:nstep))*dt
  ! !### COMMENTING THIS LINE THE RESULTS ARE IDENTICAL WITH THE TWO METHODS OF UPDATE ###
  ! !forall(i=0:nstep)G0ret(i,i)=-xi !???
  ! !#####################################################################################
  ! G0adv=conjg(transpose(G0ret))
  ! !
  ! !G0%less = GammaR^-1 * Gless * GammaA^-1  -  gR * Sigma%less * gA
  ! G0%less(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),matmul(locG%less(0:nstep,0:nstep),&
  !      conjg(transpose(GammaRet(0:nstep,0:nstep))))*dt)*dt -&
  !      matmul(G0ret(0:nstep,0:nstep),matmul(Sigma%less(0:nstep,0:nstep),G0adv(0:nstep,0:nstep))*dt)*dt
  ! !
  ! !G0%gtr  = GammaR^-1 * Ggtr * GammaA^-1   -  gR * Sigma%gtr * gA
  ! G0%gtr(0:nstep,0:nstep)  = matmul(GammaRet(0:nstep,0:nstep),matmul(locG%gtr(0:nstep,0:nstep),&
  !      conjg(transpose(GammaRet(0:nstep,0:nstep))))*dt)*dt  -&
  !      matmul(G0ret(0:nstep,0:nstep),matmul(Sigma%gtr(0:nstep,0:nstep),G0adv(0:nstep,0:nstep))*dt)*dt
  
  
  
  ! !Matrix update, from testKELDYSHMATGF3
  ! !Build Gloc matrix
  ! allocate(mat_locG(0:2*nstep+1,0:2*nstep+1))
  ! mat_locG = build_keldysh_matrix_gf(locG,nstep)
  ! !
  ! !Build Sigma matrix
  ! allocate(mat_Sigma(0:2*nstep+1,0:2*nstep+1))
  ! mat_Sigma = build_keldysh_matrix_gf(Sig,nstep)
  ! !
  ! !Allocate space for other matrices:
  ! allocate(mat_Delta(0:2*nstep+1,0:2*nstep+1))
  ! allocate(mat_Gamma(0:2*nstep+1,0:2*nstep+1))
  ! allocate(mat_G0(0:2*nstep+1,0:2*nstep+1))
  ! !
  ! mat_Delta=zero ; forall(i=0:2*nstep+1)mat_Delta(i,i)=One/dt
  ! mat_Gamma = mat_Delta + matmul(mat_Sigma,mat_locG)*dt
  ! mat_Gamma = mat_Gamma*dt**2
  ! call mat_inversion(mat_Gamma)
  ! mat_G0  = matmul(mat_locG,mat_Gamma)*dt
  ! !
  ! G0%less = -mat_G0(0:Nstep,Nstep+1:2*Nstep+1)
  ! G0%gtr  =  mat_G0(Nstep+1:2*Nstep+1,0:Nstep)
  ! !
  ! deallocate(mat_locG,mat_Sigma,mat_G0,mat_Delta,mat_Gamma)
  
  
  
  allocate(mat_locG(t1min:t3max,t1min:t3max))
  mat_locG(0:,0:) = build_kbm_matrix_gf(locG,Nstep,Ltau)

  allocate(mat_Sigma(t1min:t3max,t1min:t3max))
  mat_Sigma(0:,0:) = build_kbm_matrix_gf(Sigma,Nstep,Ltau)

  forall(i=t1min:t3max,j=t1min:t3max)mat_locG(i,j)=dtloc(i)*mat_locG(i,j)*dtloc(j)
  call matrix_inverse(mat_locG(0:,0:))

  allocate(mat_calG(t1min:t3max,t1min:t3max))
  mat_calG(0:,0:) = mat_locG(0:,0:) + mat_Sigma(0:,0:)

  forall(i=t1min:t3max,j=t1min:t3max)mat_calG(i,j)=dtloc(i)*mat_calG(i,j)*dtloc(j)
  call matrix_inverse(mat_calG(0:,0:))

  call scatter_kbm_matrix_gf(mat_calG(0:,0:),Nstep,Ltau,G0)
  
