FC=mpif90
EXE   = neqDMFT
DIREXE= $(HOME)/.bin

#=========================================================================
include lib.mk
#=========================================================================

.SUFFIXES: .f90 
OBJS =  CONTOUR_GF.o VARS_GLOBAL.o ELECTRIC_FIELD.o BATH.o EQUILIBRIUM.o IPT_NEQ.o FUNX_NEQ.o KADANOFBAYM.o

#=================STANDARD COMPILATION====================================
all:	FLAG=$(STD)
	ARGS=$(LIBDMFT) $(SFMODS) $(SFLIBS)
	BRANCH=  $(shell git rev-parse --abbrev-ref HEAD)
all: 	compile

#================OPTIMIZED COMPILATION====================================
opt: 	FLAG=$(OPT)
	ARGS=$(LIBDMFT) $(SFMODS) $(SFLIBS)
	BRANCH=  $(shell git rev-parse --abbrev-ref HEAD)
opt: 	compile

#================DEBUGGIN COMPILATION=====================================
debug:	FLAG=$(DEB)
	ARGS=$(LIBDMFT_DEB) $(SFMODS_DEB) $(SFLIBS_DEB)
	BRANCH=  $(shell git rev-parse --abbrev-ref HEAD)
debug:	compile


compile: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC)  $(FLAG) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
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
