  do i=1,command_argument_count()
     nml_var=get_cmd_variable(i)
     select case(nml_var%name) !all variables here are UPPER case
        !variables:
     case("DT")		;read(nml_var%value,*)dt
     case("BETA")	;read(nml_var%value,*)beta
     case("U")		;read(nml_var%value,*)U
     case("EFIELD")	;read(nml_var%value,*)Efield
     case("VPD")	;read(nml_var%value,*)Vpd
     case("TS")		;read(nml_var%value,*)ts
     case("NSTEP")	;read(nml_var%value,*)nstep
     case("NLOOP")	;read(nml_var%value,*)nloop
     case("EPS_ERROR")	;read(nml_var%value,*)eps_error
     case("NSUCCESS")	;read(nml_var%value,*)Nsuccess
     case("WEIGHT")	;read(nml_var%value,*)weight
        !Efield:
     case("EX")		;read(nml_var%value,*)Ex
     case("EY")		;read(nml_var%value,*)Ey
     case("T0")		;read(nml_var%value,*)t0
     case("T1")		;read(nml_var%value,*)t1
     case("TAU0")	;read(nml_var%value,*)tau0
     case("W0")		;read(nml_var%value,*)w0
     case("OMEGA0")	;read(nml_var%value,*)omega0
     case("E1")	        ;read(nml_var%value,*)E1
     case("FIELD_PROFILE");read(nml_var%value,*)field_profile
        !flags:
     case("IRDEQ")	        ;read(nml_var%value,*)irdeq
     case("METHOD")		;read(nml_var%value,*)method
     case("UPDATE_WFFTW")	;read(nml_var%value,*)update_wfftw
     case("SOLVE_WFFTW")	;read(nml_var%value,*)solve_wfftw
     case("PLOTVF")		;read(nml_var%value,*)plotVF
     case("PLOT3D")		;read(nml_var%value,*)plot3D
     case("FCHI")		;read(nml_var%value,*)fchi
     case("EQUENCH")		;read(nml_var%value,*)equench
        !parameters:
     case("L")		;read(nml_var%value,*)L
     case("LTAU")	;read(nml_var%value,*)Ltau
     case("LMU")	;read(nml_var%value,*)Lmu
     case("LKREDUCED")	;read(nml_var%value,*)Lkreduced
     case("WBATH")	;read(nml_var%value,*)wbath
     case("BATH_TYPE")	;read(nml_var%value,*)bath_type
     case("EPS")	;read(nml_var%value,*)eps
     case("IRDG0wFILE")	;read(nml_var%value,*)irdG0wfile
     case("IRDNKFILE")	;read(nml_var%value,*)irdNkfile
     case("IRDSFILE")	;read(nml_var%value,*)irdSfile
        !LatticeN:
     case("NX")		;read(nml_var%value,*)Nx
     case("NY")		;read(nml_var%value,*)Ny
        !Quench:
     case("IQUENCH")	;read(nml_var%value,*)iquench
     case("BETA0")	;read(nml_var%value,*)beta0
     case("U0")         ;read(nml_var%value,*)U0
     case("XMU0")       ;read(nml_var%value,*)xmu0
     case("SOLVE_EQ")   ;read(nml_var%value,*)solve_eq
     case("G0LOC_GUESS");read(nml_var%value,*)g0loc_guess
     case default
        print*,"No corresponding variable in NML"
     end select
  enddo
