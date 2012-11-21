#=========================================================================
include lib.mk
#=========================================================================
FC=$(SFMPI)/mpif90
EXE   = neqDMFT
DIREXE= $(HOME)/.bin

.SUFFIXES: .f90 
OBJS =  CONTOUR_GF.o VARS_GLOBAL.o EQUILIBRIUM.o ELECTRIC_FIELD.o BATH.o  IPT_NEQ.o UPDATE_WF.o KADANOFBAYM.o

ARGS=$(LIBDMFT) $(SFMODS) $(SFLIBS)
ARGS_DEB=$(LIBDMFT_DEB) $(SFMODS_DEB) $(SFLIBS_DEB)
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)

#=================STANDARD COMPILATION====================================
all: FLAG=$(STD)
all: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC)  $(FLAG) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)

#================OPTIMIZED COMPILATION====================================
opt: FLAG=$(OPT)
opt: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FLAG) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)

#================DEBUGGIN COMPILATION=====================================
debug: FLAG=$(DEB)
debug: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC)  $(FLAG) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS_DEB)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)


#==============DATA EXTRACTION======================================
data: 	version $(OBJS)
	@echo " ........... compile: getdata ........... "
	${FC} ${STD} $(OBJS) get_data_$(EXE).f90 -o ${DIREXE}/get_data_$(EXE)_$(BRANCH) $(LIBDMFT) ${DSL_MODS} ${SFMODS} ${SFLIBS} ${DSL_LIBS}
	@echo ""
	@echo " ...................... done .............................. "


.f90.o:	
	$(FC) $(FLAG) -c $< $(SFMODS) 



#=============CLEAN ALL===================================================
clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc


#=========================================================================
include version.mk
#=========================================================================
