#Makefile --- WMAP-9 Likelihood code...

#CFITSIO=/usr/local
CFITSIO=/home/xyh/Tools/

#
#			Non-test object files.
#
WMAPLIB = libwmap9.a
OBJS = read_archive_map.o \
	read_fits.o \
	healpix_types.o \
	br_mod_dist.o \
	WMAP_9yr_options.o \
	WMAP_9yr_util.o \
	WMAP_9yr_gibbs.o \
	WMAP_9yr_tt_pixlike.o \
	WMAP_9yr_tt_beam_ptsrc_chisq.o \
	WMAP_9yr_teeebb_pixlike.o \
	WMAP_9yr_tetbeebbeb_pixlike.o \
	WMAP_9yr_likelihood.o
#
#			General Commands.
#
DIFF = diff -w
RM = rm -f

# See the CHANGES files for a description of these options
WMAPFLAGS  = -DOPTIMIZE
#WMAPFLAGS += -DUSE_LOWELL_TBEB    # turns on low-l maxlike TB/EB
#WMAPFLAGS += -DUSE_HIGHELL_TB     # turns on high-l master TB
#WMAPFLAGS += -DFASTERTT           # speed up matrix low-l TT
#WMAPFLAGS += -DTIMING             # print out timing stats for profiling

## Linux/Intel compiler and MKL libraries

#F90    = ifort
#FFLAGS = -O2 -fpic $(WMAPFLAGS) 
#INCS   = -I. -I$(CFITSIO)/includei
#INTEL_LIB=/export/apps/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64/
#LIBS = -L. -L$(INTEL_LIB) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lmkl_mc3 -lmkl_def -lmkl_lapack -L$(CFITSIO) -lcfitsio
#LIBS = -L. -L$(INTEL_LIB) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lmkl_mc3 -lmkl_def -lmkl_lapack95_lp64 -L$(CFITSIO) -lcfitsio

F90		= gfortran-7
FFLAGS	= -O2 -fpic $(WMAPFLAGS)
LIBS	= -lcfitsio -llapack -lblas -lgfortran


PROGRAMS = test 

all: $(PROGRAMS) $(WMAPLIB) 

check: test
	./test

$(WMAPLIB): $(OBJS)
	ar r $@ $^

%: $(OBJS) %.o
	$(F90) $(FFLAGS) -o $@ $^ $(LIBS)

%.o: %.f90
	$(F90) $(FFLAGS) $(INCS) -c -o $@ $<

%.o: %.F90
	$(F90) $(FFLAGS) $(INCS) -c -o $@ $<

clean:
	$(RM) *.o *.mod *.log *~ *.a 

distclean: clean
	$(RM) *.a $(PROGRAMS)




