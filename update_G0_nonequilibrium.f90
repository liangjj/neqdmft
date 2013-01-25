  !=======Component by component inversion==========================
  if(FF)then
     call msg("update with method 1: equations for <,>")
     forall(i=0:nstep,j=0:nstep)
        locGret(i,j)= heaviside(t(i)-t(j))*(locG%gtr(i,j) - locG%less(i,j))
        Sret(i,j)   = heaviside(t(i)-t(j))*(Sigma%gtr(i,j) - Sigma%less(i,j))
     end forall

     ! !G0ret = [\11 + Gret * Sret]^-1 * Gret
     Uno=zero  ; forall(i=0:nstep)Uno(i,i)=One/dt
     GammaRet(0:nstep,0:nstep) = Uno+matmul(locGret(0:nstep,0:nstep),Sret(0:nstep,0:nstep))*dt
     GammaRet(0:nstep,0:nstep) = GammaRet(0:nstep,0:nstep)*dt**2
     call matrix_inverse(GammaRet(0:nstep,0:nstep))
     G0ret(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),locGret(0:nstep,0:nstep))*dt
     !### COMMENTING THIS LINE THE RESULTS ARE IDENTICAL WITH THE THREE METHODS OF UPDATE ###
     !forall(i=0:nstep)G0ret(i,i)=-xi !???
     !#####################################################################################
     G0adv=conjg(transpose(G0ret))

     !G0less = GammaR^-1 * Gless * GammaA^-1  -  gR * Sless * gA
     G0%less(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),matmul(locG%less(0:nstep,0:nstep),&
          conjg(transpose(GammaRet(0:nstep,0:nstep))))*dt)*dt -&
          matmul(G0ret(0:nstep,0:nstep),matmul(Sigma%less(0:nstep,0:nstep),G0adv(0:nstep,0:nstep))*dt)*dt

     !G0gtr  = GammaR^-1 * Ggtr * GammaA^-1   -  gR * Sgtr * gA
     G0%gtr(0:nstep,0:nstep)  = matmul(GammaRet(0:nstep,0:nstep),matmul(locG%gtr(0:nstep,0:nstep),&
          conjg(transpose(GammaRet(0:nstep,0:nstep))))*dt)*dt  -&
          matmul(G0ret(0:nstep,0:nstep),matmul(Sigma%gtr(0:nstep,0:nstep),G0adv(0:nstep,0:nstep))*dt)*dt
  endif



  !Matrix update, from testKELDYSHMATGF3
  if(FF)then
     call msg("update with method 2: inversion and matrix-multiplication")
     !Build Gloc matrix
     allocate(mat_locG(0:2*nstep+1,0:2*nstep+1))
     mat_locG(0:,0:) = build_keldysh_matrix_gf(locG,nstep)

     !Build Sigma matrix
     allocate(mat_Sigma(0:2*nstep+1,0:2*nstep+1))
     mat_Sigma(0:,0:) = build_keldysh_matrix_gf(Sigma,nstep)

     !Allocate space for other matrices:
     allocate(mat_Delta(0:2*nstep+1,0:2*nstep+1))
     allocate(mat_Gamma(0:2*nstep+1,0:2*nstep+1))
     allocate(mat_G0(0:2*nstep+1,0:2*nstep+1))

     mat_Delta(0:,0:)=zero ; forall(i=0:2*nstep+1)mat_Delta(i,i)=One/dt
     mat_Gamma(0:,0:) = mat_Delta(0:,0:) + matmul(mat_Sigma(0:,0:),mat_locG(0:,0:))*dt
     mat_Gamma(0:,0:) = mat_Gamma(0:,0:)*dt**2
     call matrix_inverse(mat_Gamma(0:,0:))
     mat_G0(0:,0:)  = matmul(mat_locG(0:,0:),mat_Gamma(0:,0:))*dt

     G0%less(0:,0:) = -mat_G0(0:Nstep,Nstep+1:2*Nstep+1)
     G0%gtr(0:,0:)  =  mat_G0(Nstep+1:2*Nstep+1,0:Nstep)

     deallocate(mat_locG,mat_Sigma,mat_G0,mat_Delta,mat_Gamma)
  endif


  !Matrix update, from testKELDYSHMATGF4
  if(TT)then
     call msg("update with method 3: direct inversion of keldysh matrix GF")
     !Build Gloc matrix
     allocate(mat_locG(0:2*nstep+1,0:2*nstep+1))
     mat_locG(0:,0:) = build_keldysh_matrix_gf(locG,nstep)

     !Build Sigma matrix
     allocate(mat_Sigma(0:2*nstep+1,0:2*nstep+1))
     allocate(mat_SigmaHF(0:2*nstep+1,0:2*nstep+1))
     mat_Sigma(0:,0:) = build_keldysh_matrix_gf(Sigma,nstep)

     forall(i=0:nstep)
        mat_SigmaHF(i,i)                 = SigmaHF(i)/dt
        mat_SigmaHF(nstep+1+i,nstep+1+i )= SigmaHF(i)/dt
     end forall

     !Get G0 matrix:
     allocate(mat_G0(0:2*nstep+1,0:2*nstep+1))
     mat_locG(0:,0:) = mat_locG(0:,0:)*dt**2
     call matrix_inverse(mat_locG(0:,0:))
     mat_G0(0:,0:) = mat_locG(0:,0:) + mat_Sigma(0:,0:) - mat_SigmaHF(0:,0:)
     mat_G0(0:,0:) = mat_G0(0:,0:)*dt**2
     call matrix_inverse(mat_G0(0:,0:))

     G0%less(0:,0:) = -mat_G0(0:Nstep,Nstep+1:2*Nstep+1)
     G0%gtr(0:,0:)  =  mat_G0(Nstep+1:2*Nstep+1,0:Nstep)

     deallocate(mat_locG,mat_Sigma,mat_G0)
  endif




