EXE   = neqDMFT
DIREXE= $(HOME)/.bin

#=========================================================================

include $(HOME)/lib/lib.mk
#=========================================================================

OBJS     =  VARS_GLOBAL.o ELECTRIC_FIELD.o BATH.o FUNX_NEQ.o KADANOFBAYM.o
OBJS_OPT = VARS_GLOBAL_OPT.o ELECTRIC_FIELD_OPT.o BATH_OPT.o FUNX_NEQ_OPT.o KADANOFBAYM_OPT.o
OBJS_DEB = VARS_GLOBAL_DEB.o ELECTRIC_FIELD_DEB.o BATH_DEB.o FUNX_NEQ_DEB.o KADANOFBAYM_DEB.o

#=================STANDARD COMPILATION====================================
all: 	version $(OBJS)
	@echo " ........... compile: normal ........... "
	$(FC)  $(STD) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE) $(MODS) $(LIBS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)



#================OPTIMIZED COMPILATION====================================
opt: 	version $(OBJS_OPT)
	@echo " ........... compile: optimized   ........... "
	$(FC) $(OPT) $(OBJS_OPT) $(EXE).f90 -o $(DIREXE)/$(EXE) $(MODS) $(LIBS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)



#================DEBUGGIN COMPILATION=====================================
debug: 	version $(OBJS_DEB)
	@echo " ........... compile: debug   ........... "
	$(FC) $(DEB) $(OBJS_DEB) $(EXE).f90 -o $(DIREXE)/$(EXE) $(MODS_DEB) $(LIBS_DEB)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)



#==============DATA EXTRACTION======================================
data: 	version $(OBJS)
	@echo " ........... compile: getdata ........... "
	${FC} ${STD} $(OBJS) get_data_$(EXE).f90 -o ${DIREXE}/get_data_$(EXE) ${DSL_MODS} ${MODS} ${DSL_LIBS} ${X11} ${LIBS}
	@echo ""
	@echo ""
	@echo " ...................... done .............................. "



#==============FIT CURRENT======================================
fit: 	$(OBJS_DEB)
	@echo " ........... compile: getdata ........... "
	${FC} ${DEB} $(OBJS_DEB) testFITJ.f90 -o ${DIREXE}/testFITJ ${DSL_MODS} ${MODS} ${DSL_LIBS} ${X11} ${LIBS}
	@echo ""
	@echo ""
	@echo " ...................... done .............................. "


VARS_GLOBAL.o: VARS_GLOBAL.f90
	$(FC) $(STD) -c VARS_GLOBAL.f90 $(MODS)
BATH.o: BATH.f90 
	$(FC) $(STD) -c BATH.f90 $(MODS)
ELECTRIC_FIELD.o: ELECTRIC_FIELD.f90 
	$(FC) $(STD) -c ELECTRIC_FIELD.f90 $(MODS)
FUNX_NEQ.o: FUNX_NEQ.f90
	$(FC) $(STD) -c FUNX_NEQ.f90 $(MODS)
KADANOFBAYM.o: KADANOFBAYM.f90 
	$(FC) $(STD) -c KADANOFBAYM.f90 $(MODS)


VARS_GLOBAL_OPT.o: VARS_GLOBAL.f90
	$(FC) $(OPT) -c VARS_GLOBAL.f90 $(MODS) -o VARS_GLOBAL_OPT.o
BATH_OPT.o: BATH.f90 
	$(FC) $(OPT) -c BATH.f90 $(MODS) -o BATH_OPT.o
ELECTRIC_FIELD_OPT.o: ELECTRIC_FIELD.f90 
	$(FC) $(OPT) -c ELECTRIC_FIELD.f90 $(MODS) -o ELECTRIC_FIELD_OPT.o
FUNX_NEQ_OPT.o: FUNX_NEQ.f90 
	$(FC) $(OPT) -c FUNX_NEQ.f90 $(MODS) -o FUNX_NEQ_OPT.o
KADANOFBAYM_OPT.o: KADANOFBAYM.f90 
	$(FC) $(OPT) -c KADANOFBAYM.f90 $(MODS) -o KADANOFBAYM_OPT.o



VARS_GLOBAL_DEB.o: VARS_GLOBAL.f90
	$(FC) $(DEB) -c VARS_GLOBAL.f90 $(MODS_DEB) -o VARS_GLOBAL_DEB.o
BATH_DEB.o: BATH.f90 
	$(FC) $(DEB) -c BATH.f90 $(MODS_DEB) -o BATH_DEB.o
ELECTRIC_FIELD_DEB.o: ELECTRIC_FIELD.f90 
	$(FC) $(DEB) -c ELECTRIC_FIELD.f90 $(MODS_DEB) -o ELECTRIC_FIELD_DEB.o
FUNX_NEQ_DEB.o: FUNX_NEQ.f90 
	$(FC) $(DEB) -c FUNX_NEQ.f90 $(MODS_DEB) -o FUNX_NEQ_DEB.o
KADANOFBAYM_DEB.o: KADANOFBAYM.f90 
	$(FC) $(DEB) -c KADANOFBAYM.f90 $(MODS_DEB) -o KADANOFBAYM_DEB.o



#=============CLEAN ALL===================================================
clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc


#=============REVISION ===================================================
include $(HOME)/lib/version.mk
