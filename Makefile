SHELL = /bin/bash

## OPTIONS  ##
OPT 	+= -DNFWC_DUFFY08	# alternate fit to concentr. param

OPT += -DBETA=0.54

#OPT     += -DPARABOLA       # merge in a parabola
OPT	+= -DCOMET			# merge like a comet, ball+tail (recommended)

#OPT	+= -DDOUBLE_BETA_COOL_CORES

#OPT 	+= -DGIVEPARAMS		# more merger parameters in .par file

#OPT		+= -DNO_RCUT_IN_T		# set Rcut very large

#OPT	+= -DSUBSTRUCTURE	# add substructure
#OPT += -DSUBHOST=1		# host halos
#OPT	+= -DSLOW_SUBSTRUCTURE	# put subhalos on Hernquist orbits
#OPT += -DREPORTSUBHALOS		# print info about all subhaloes

#OPT += -DADD_THIRD_SUBHALO  # manually set the first subhalo mass, pos, vel
#OPT  += -DTHIRD_HALO_ONLY

#OPT += -DSPH_CUBIC_SPLINE 	# for use with Gadget2

## Target Computer ##
ifndef SYSTYPE
SYSTYPE := $(shell hostname)
endif

# standard systypes
CC       = gcc
OPTIMIZE = -Wall -g -O2
GSL_INCL = $(CPPFLAGS)
GSL_LIBS = $(LDFLAGS)

ifeq ($(SYSTYPE),DARWIN)
CC      	=  icc
OPTIMIZE	=-fast -g -m64  -xhost
GSL_INCL 	= $(CPPFLAGS)
GSL_LIBS	= -L/Users/jdonnert/Dev/lib
endif

ifeq ($(SYSTYPE),MSI)
CC      	= icc
OPTIMIZE	= -Wall -g  -O3 -xhost
GSL_INCL 	= 
GSL_LIBS	= 
FFTW_LIBS 	= 
FFTW_INCL 	=
endif

ifeq ($(SYSTYPE),mach64.ira.inaf.it)
CC      	= gcc
OPTIMIZE	= -O2 -Wall -g  -m64 -march=native -mtune=native -mprefer-avx128 -fopenmp  -minline-all-stringops -fprefetch-loop-arrays --param prefetch-latency=300 -funroll-all-loops
GSL_INCL 	= -I/homes/donnert/Libs/include
GSL_LIBS	= -L/homes/donnert/Libs/lib
FFTW_LIBS 	= 
FFTW_INCL 	=
endif

## TARGET ##

EXEC = Toycluster

## FILES ##

SRCDIR	= src/

OBJFILES = main.o aux.o positions.o velocities.o temperature.o \
		   magnetic_field.o io.o unit.o cosmo.o setup.o  tree.o \
		   sph.o wvt_relax.o substructure.o ids.o sort.o peano.o

OBJS	= $(addprefix $(SRCDIR),$(OBJFILES))

INCLFILES = globals.h proto.h io.h tree.h sph.h macro.h sort.h peano.h \
			../Makefile

INCL	= $(addprefix $(SRCDIR),$(INCLFILES))

CFLAGS 	= -std=c99 -fopenmp $(OPTIMIZE) $(OPT) $(GSL_INCL) $(FFTW_INCL)

LINK	= -lm -lgsl -lgslcblas $(GSL_LIBS) $(FFTW_LIBS)

## RULES ## 

$(EXEC)	: $(OBJS)
	@echo SYSTYPE=$(SYSTYPE)
	$(CC) $(CFLAGS) $(OBJS) $(LINK) -o $(EXEC)
	@cd src && ctags *.[ch]

$(OBJS)	: $(INCL)

clean	: 
	rm -f  $(OBJS) $(EXEC)
