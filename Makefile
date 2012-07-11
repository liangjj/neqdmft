FC=mpif90
EXE   = neqDMFT
DIREXE= $(HOME)/.bin

#=========================================================================
include $(SFDIR)/etc/lib.mk
include $(SFDIR)/etc/libdmft.mk
#=========================================================================

OBJS     =  VARS_GLOBAL.o ELECTRIC_FIELD.o BATH.o EQUILIBRIUM.o IPT_NEQ.o FUNX_NEQ.o KADANOFBAYM.o
OBJS_OPT = VARS_GLOBAL_OPT.o ELECTRIC_FIELD_OPT.o BATH_OPT.o EQUILIBRIUM_OPT.o IPT_NEQ_OPT.o FUNX_NEQ_OPT.o KADANOFBAYM_OPT.o
OBJS_DEB = VARS_GLOBAL_DEB.o ELECTRIC_FIELD_DEB.o BATH_DEB.o EQUILIBRIUM_DEB.o IPT_NEQ_DEB.o FUNX_NEQ_DEB.o KADANOFBAYM_DEB.o

#=================STANDARD COMPILATION====================================
all: 	version $(OBJS)
	@echo " ........... compile: normal ........... "
	$(FC)  $(STD) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT) $(SFMODS) $(SFLIBS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)



#================OPTIMIZED COMPILATION====================================
opt: 	version $(OBJS_OPT)
	@echo " ........... compile: optimized   ........... "
	$(FC) $(OPT) $(OBJS_OPT) $(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT) $(SFMODS) $(SFLIBS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)



#================DEBUGGIN COMPILATION=====================================
debug: 	version $(OBJS_DEB)
	@echo " ........... compile: debug   ........... "
	$(FC) $(DEB) $(OBJS_DEB) $(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT_DEB) $(SFMODS_DEB) $(SFLIBS_DEB)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)



#==============DATA EXTRACTION======================================
data: 	version $(OBJS_DEB)
	@echo " ........... compile: getdata ........... "
	${FC} ${DEB} $(OBJS_DEB) get_data_$(EXE).f90 -o ${DIREXE}/get_data_$(EXE) $(LIBDMFT_DEB) ${DSL_MODS} ${SFMODS_DEB} ${SFLIBS_DEB} ${DSL_LIBS}
	@echo ""
	@echo ""
	@echo " ...................... done .............................. "



VARS_GLOBAL.o: VARS_GLOBAL.f90
	$(FC) $(STD) -c VARS_GLOBAL.f90 $(SFMODS)
BATH.o: BATH.f90 
	$(FC) $(STD) -c BATH.f90 $(SFMODS)
ELECTRIC_FIELD.o: ELECTRIC_FIELD.f90 
	$(FC) $(STD) -c ELECTRIC_FIELD.f90 $(SFMODS)
FUNX_NEQ.o: FUNX_NEQ.f90
	$(FC) $(STD) -c FUNX_NEQ.f90 $(SFMODS)
IPT_NEQ.o: IPT_NEQ.f90
	$(FC) $(STD) -c IPT_NEQ.f90 $(SFMODS)
KADANOFBAYM.o: KADANOFBAYM.f90 
	$(FC) $(STD) -c KADANOFBAYM.f90 $(SFMODS)
EQUILIBRIUM.o: EQUILIBRIUM.f90 
	$(FC) $(STD) -c EQUILIBRIUM.f90  $(SFMODS)



VARS_GLOBAL_OPT.o: VARS_GLOBAL.f90
	$(FC) $(OPT) -c VARS_GLOBAL.f90 $(SFMODS) -o VARS_GLOBAL_OPT.o
BATH_OPT.o: BATH.f90 
	$(FC) $(OPT) -c BATH.f90 $(SFMODS) -o BATH_OPT.o
ELECTRIC_FIELD_OPT.o: ELECTRIC_FIELD.f90 
	$(FC) $(OPT) -c ELECTRIC_FIELD.f90 $(SFMODS) -o ELECTRIC_FIELD_OPT.o
FUNX_NEQ_OPT.o: FUNX_NEQ.f90 
	$(FC) $(OPT) -c FUNX_NEQ.f90 $(SFMODS) -o FUNX_NEQ_OPT.o
IPT_NEQ_OPT.o: IPT_NEQ.f90 
	$(FC) $(OPT) -c IPT_NEQ.f90 $(SFMODS) -o IPT_NEQ_OPT.o
KADANOFBAYM_OPT.o: KADANOFBAYM.f90 
	$(FC) $(OPT) -c KADANOFBAYM.f90 $(SFMODS) -o KADANOFBAYM_OPT.o
EQUILIBRIUM_OPT.o: EQUILIBRIUM.f90 
	$(FC) $(OPT) -c EQUILIBRIUM.f90 $(SFMODS) -o EQUILIBRIUM_OPT.o



VARS_GLOBAL_DEB.o: VARS_GLOBAL.f90
	$(FC) $(DEB) -c VARS_GLOBAL.f90 $(SFMODS_DEB) -o VARS_GLOBAL_DEB.o
BATH_DEB.o: BATH.f90 
	$(FC) $(DEB) -c BATH.f90 $(SFMODS_DEB) -o BATH_DEB.o
ELECTRIC_FIELD_DEB.o: ELECTRIC_FIELD.f90 
	$(FC) $(DEB) -c ELECTRIC_FIELD.f90 $(SFMODS_DEB) -o ELECTRIC_FIELD_DEB.o
FUNX_NEQ_DEB.o: FUNX_NEQ.f90 
	$(FC) $(DEB) -c FUNX_NEQ.f90 $(SFMODS_DEB) -o FUNX_NEQ_DEB.o
IPT_NEQ_DEB.o: IPT_NEQ.f90 
	$(FC) $(DEB) -c IPT_NEQ.f90 $(SFMODS_DEB) -o IPT_NEQ_DEB.o
KADANOFBAYM_DEB.o: KADANOFBAYM.f90 
	$(FC) $(DEB) -c KADANOFBAYM.f90 $(SFMODS_DEB) -o KADANOFBAYM_DEB.o
EQUILIBRIUM_DEB.o: EQUILIBRIUM.f90 
	$(FC) $(DEB) -c EQUILIBRIUM.f90 $(SFMODS) -o EQUILIBRIUM_DEB.o


#=============CLEAN ALL===================================================
clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc


#=========================================================================
include $(SFDIR)/etc/version.mk
#=========================================================================
