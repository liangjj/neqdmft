#=========================================================================
include sfmake.inc
#=========================================================================
FC=$(SFMPI)/mpif90
EXE   = neqDMFT
DIREXE= $(HOME)/.bin
BRANCH=  $(shell git rev-parse --abbrev-ref HEAD)
.SUFFIXES: .f90 
OBJS =  CONTOUR_GF.o VARS_GLOBAL.o ELECTRIC_FIELD.o BATH.o EQUILIBRIUM.o IPT_NEQ.o UPDATE_WF.o KADANOFBAYM.o RESULTS.o

#=================STANDARD COMPILATION====================================
all:FLAG=$(STD)
    ARGS=$(LIBDMFT) $(SFMODS) $(SFLIBS)
all:compile

#================OPTIMIZED COMPILATION====================================
opt:FLAG=$(OPT)
    ARGS=$(LIBDMFT) $(SFMODS) $(SFLIBS)
opt:compile

#================DEBUGGIN COMPILATION=====================================
debug:FLAG=$(DEB)
      ARGS=$(LIBDMFT_DEB) $(SFMODS_DEB) $(SFLIBS_DEB)
debug:compile


compile: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FLAG) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)



#==============DATA EXTRACTION======================================
data:	FLAG=$(STD)
	ARGS=$(LIBDMFT) $(SFMODS) $(SFLIBS) #$(DSL_MODS) $(DSL_LIBS)
	BRANCH=  $(shell git rev-parse --abbrev-ref HEAD)
data: 	version $(OBJS)
	@echo " ........... compile: getdata ........... "
	${FC} ${FLAG} $(OBJS) get_data_$(EXE).f90 -o ${DIREXE}/get_data_$(EXE)_$(BRANCH) $(ARGS) 
	@echo ""
	@echo " ...................... done .............................. "


.f90.o:	
	$(FC) $(FLAG) -c $< $(SFMODS) 

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)
