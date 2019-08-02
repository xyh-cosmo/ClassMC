#   ============================================================================
#   Makefile for ClassMC
#   Youhua Xu
#
#	This Makefile is modified from the one included in CLASS.
#   ============================================================================

MDIR    := $(shell pwd)
WRKDIR   = $(MDIR)/build
BINDIR   = $(MDIR)/bin

CLASS_DIR = $(MDIR)/source/class_DDE
WMAP7_DIR = $(MDIR)/source/wmap_likelihood_v4p1
WMAP9_DIR = $(MDIR)/source/wmap_likelihood_v5


.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR); fi;
	touch $(WRKDIR)/.base;
	if ! [ -e $(BINDIR) ]; then mkdir $(BINDIR); fi;
	touch $(BINDIR)/.base;

vpath %.c main:source:source/datasets:source/likelihoods:source/Wrapper:source/priorss:source/test
vpath %.cc main:source:source/datasets:source/likelihoods:source/Wrapper:source/priors:source/test
vpath %.cpp main:source:source/datasets:source/likelihoods:source/Wrapper:source/priors:source/test
vpath %.o build
vpath .base build

################################################################################
#	flags:
#	CCFLAG 		-- compilation flags
#	LDFLAG 		-- library link flags
#	OPTFLAG		-- optimization flags
#	DEBUGFLAG	-- debug flags
CCFLAG      = #-g
LDFLAG      =
OPTFLAG     =

DEBUGFLAG	=
#DEBUGFLAG	+= -D_DEBUG_CHI2_

CLASSMC_CCFLAG	=
# CLASSMC_CCFLAG   += -D_DEBUG_PAIR_MCMC_PARAMS_


# which w(z) and weff(z) approximation method is to be used ?
#CLASS_WEFF_CCFLAG = -D_DDE_APPROX_OPT_0_	# use linear interpolation approximation for w(z)
CLASS_WEFF_CCFLAG = -D_DDE_APPROX_OPT_1_	# use piece-wise constant approximation for w(z)

# for debug weff(z)
#CLASS_WEFF_CCFLAG += -D_DEBUG_WEFF_

################################################################################
# Compiler options: if you want to use intel compiler, set this flag to T
ERR = $(shell which icpc >/dev/null; echo $$?)
ifeq "$(ERR)" "0"
	HAS_INTEL	= T
else
	HAS_INTEL	= F
endif

# MK = Makefile
ifeq ($(shell uname -s),Linux)
CC			= gcc
CXX			= g++
FC			= gfortran
CLASS_MK	= Makefile.gnu
WMAP7_MK    = Makefile.gfortran
WMAP9_MK    = Makefile.gfortran
LIBS        = -lgfortran
endif

ifeq ($(shell uname -s),Darwin)
CC			= clang
CXX			= clang++
FC			= gfortran
CLASS_MK	= Makefile.Mac
WMAP7_MK    = Makefile.gfortran
WMAP9_MK    = Makefile.gfortran
LIBS        = -lgfortran
endif

################################
# force usage of intel compilers
################################
ifeq ($(HAS_INTEL),T)
$(info **********************************)
$(info **** DETECTED INTEL COMPILERS ****)
$(info **********************************)
CC			= icc
CXX			= icpc
FC			= ifort
CLASS_MK	= Makefile.intel
WMAP7_MK    = Makefile.intel
WMAP9_MK    = Makefile.intel
endif

# by default we use mpic++ to compile ClassMC
MPI_CC      = mpic++

CCFLAG  	+= $(CLASSMC_CCFLAG) $(CLASS_WEFF_CCFLAG)
# CCFLAG  	+= -Wall
# CCFLAG		+= -openmp

#########################################################
# uncomment the following to suppress noise warnings !!!!
#########################################################
ifneq ($(HAS_INTEL),T)
    ifeq ($(shell uname -s),Darwin)
        CCFLAG	+= -Wno-literal-range
    else
        CCFLAG	+= -Wno-literal-suffix
    endif
    CCFLAG		+= -DUSE_GFORTRAN   # this flag is used in WMAP9
#	CCFLAG 		+= -Wno-strict-overflow -Wno-format-overflow
else
    CCFLAG		+= -DUSE_INTEL      # this flag is used in WMAP9
endif

CCFLAG		+= -Wno-write-strings # -nofor-main #-g
# CCFLAG 		+= -Wformat-overflow=5 -Wstrict-overflow=5
# CCFLAG		+= -fno-stack-protector
OPTFLAG		+= -O3
#OPTFLAG     += -ffast-math #( not recognized by intel compiler )

#########################################################################################################
# LDFLAG is set to avoid the following warning:
# ld: warning: could not create compact unwind for ___wmap_likelihood_9yr_MOD_wmap_likelihood_compute:
# stack subq instruction is too different from dwarf stack size
# solution is found from: https://github.com/ericmandel/pyds9/issues/17

ifeq ($(shell uname -s),Darwin)
LDFLAG		= -Wl,-no_compact_unwind
endif

####################################################################################################
# which WMAP likelihood to be used?

#--------------------------------------------------------------------------------------------------
# if you want to use WMAP7, then uncomment the following 4 lines and comment out those for WMAP9
#--------------------------------------------------------------------------------------------------
CCFLAG += -D_USE_WMAP7_
WMAP_LIB = libwmap7.a
OBJ_Data_WMAP= Data_WMAP7.o
OBJ_Likelihood_WMAP = Likelihood_WMAP7.o

#--------------------------------------------------------------------------------------------------
# if you want to use WMAP9, then uncomment the following 4 lines and comment out those for WMAP7
#--------------------------------------------------------------------------------------------------
# CCFLAG += -D_USE_WMAP9_
# WMAP_LIB = libwmap9.a
# OBJ_Data_WMAP= Data_WMAP9.o
# OBJ_Likelihood_WMAP = Likelihood_WMAP9.o


# if you choose MCMC engine other than IMCMC, uncomment the one of ther other lines
LIBS		+= -limcmc -L${CLIK_PATH}/lib/ -lclik -larmadillo -lgsl -lgslcblas -llapack -lblas -lcfitsio
#LIBS		+= -limcmc -L${CLIK_PATH}/lib/ -lclik -larmadillo -lgsl -lgslcblas -llapack -lblas -lcfitsio -lmultinest
#LIBS		+= -limcmc -L${CLIK_PATH}/lib/ -lclik -larmadillo -lgsl -lgslcblas -llapack -lblas -lcfitsio -lmultinest_mpi -lchord

INCLUDES = -I ${CLASS_DIR}/include
INCLUDES += -I ${CLIK_PATH}/include
INCLUDES += -I ../include

#   ClassMC is written in C++ and it needs mpic++ for compilation
%.o: %.cpp .base
	cd $(WRKDIR); $(MPI_CC) $(CCFLAG) $(OPTFLAG) $(DEBUGFLAG) $(INCLUDES) -c ../$< -o $*.o

#   Wrappers are placed in *.cc files, need C++ compiler, non-parallel
%.o: %.cc .base
	cd $(WRKDIR); $(MPI_CC) $(CCFLAG) $(OPTFLAG) $(DEBUGFLAG) $(INCLUDES) -c ../$< -o $*.o

%.o: %.c .base
	cd $(WRKDIR); $(MPI_CC) $(CCFLAG) $(OPTFLAG) $(DEBUGFLAG) $(INCLUDES) -c ../$< -o $*.o

# wrapper to class
CLASSWRAPPER = 	ClassEngine.o 	\
				Engine.o 		\
				Cosmology.o		\
				CosmologyTest.o

# date structures
DATA = 	Data_Planck2015.o 	\
		Data_SNE_JLA.o 		\
		Data_SNE_JLA_Mock.o \
		Data_SNE_UNION.o 	\
		Data_SNE_SNLS.o 	\
		Data_SNE_WFIRST.o 	\
		Data_Hubble.o 		\
		Data_HST.o 			\
		Data_RSD.o		\
		Data_Age.o 			\
		Data_BAO.o			\
		$(OBJ_Data_WMAP)		\
		DataList.o

# likelihood functions or wrappers
LIKE = 	Likelihood_Planck2015.o 	\
		Likelihood_Age.o 			\
		Likelihood_BAO.o			\
		Likelihood_Hubble.o 		\
		Likelihood_HST.o 			\
		Likelihood_RSD.o		\
		Likelihood_SNE_Union2.1.o	\
		Likelihood_SNE_SNLS3.o 		\
		Likelihood_SNE_JLA.o 		\
		Likelihood_SNE_WFIRST.o 	\
		$(OBJ_Likelihood_WMAP)			\
		Likelihoods.o

# priors
PRIOR = Prior_DDE_CPZ.o CMB_Dist_Prior.o

# some simple mathematical tools
MATH = Math.o

#
MISC = 	DataList.o 	\
		ParamList.o \
		Misc.o 		\
		JLA_ini.o

CLASS_LIB = ${CLASS_DIR}/libclass.a
WMAP7_LIB = ${WMAP7_DIR}/libwmap7.a
WMAP9_LIB = ${WMAP9_DIR}/libwmap9.a


# choose proper version you want (do not forget to change LIBS shown above):
#1) MCMC engine is provided by IMCMC
CLASSMC_OBJS = ClassMC.o ${CLASSWRAPPER} ${DATA} ${LIKE} ${PRIOR} ${MATH} ${MISC}

#2) MCMC engine is provided by MultiNest
#CLASSMC_OBJS = ClassMC-nest.o ${CLASSWRAPPER} ${DATA} ${LIKE} ${PRIOR} ${MATH} ${MISC}

#3) MCMC engine is provided by PolyChord
#CLASSMC_OBJS = ClassMC-chord.o ${CLASSWRAPPER} ${DATA} ${LIKE} ${PRIOR} ${MATH} ${MISC}

all: ClassMC  libclass.a $(WMAP_LIB)

libclass.a:
	cd ${CLASS_DIR}; make DDEFLAG="$(CLASS_WEFF_CCFLAG)" --makefile=${CLASS_MK} ; cd -; cp ${CLASS_LIB} ${WRKDIR};

libwmap7.a:
	cd ${WMAP7_DIR}; make --makefile=${WMAP7_MK} libwmap7.a; cd -; cp ${WMAP7_LIB} ${WRKDIR};

libwmap9.a:
	cd ${WMAP9_DIR}; make --makefile=${WMAP9_MK} libwmap9.a; cd -; cp ${WMAP9_LIB} ${WRKDIR};

######## The following are the required executables #########

ClassMC:$(CLASSMC_OBJS) libclass.a $(WMAP_LIB)
	$(MPI_CC) $(CCFLAG) $(OPTFLAG) $(LDFLAG) -o ${BINDIR}/$@ $(addprefix build/,$(notdir $^)) ${LIBS}

#   ================================================================================================
.PHONY:clean tidy run
clean:
	rm -rf $(WRKDIR);
clean_wmap7:
	cd $(WMAP7_DIR); make distclean; cd -;
clean_wmap9:
	cd $(WMAP9_DIR); make distclean; cd -;
clean_class:
	cd $(CLASS_DIR); make clean; cd -;
tidy:
	make clean;
	rm -rf $(BINDIR);
	make clean_wmap7; make clean_wmap9;
	make clean_class
