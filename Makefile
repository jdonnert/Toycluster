SHELL = /bin/bash

## OPTIONS  ##
OPT 	+= -DNFWC_DUFFY08	# alternate fit to concentr. param

OPT     += -DPARABOLA       # merge in a parabula
#OPT		+= -DCOMET			# merge like a comet, ball+tail (recommended)

#OPT 	+= -DGIVEPARAMS		# more merger parameters in .par file

#OPT	+= -DSUBSTRUCTURE	# add substructure
#OPT	+= -DSLOW_SUBSTRUCTURE	# put subhalos on Hernquist orbits
#OPT += -DREPORTSUBHALOS		# print info about all subhaloes

#OPT += -DADD_THIRD_SUBHALO  # manually set the first subhalo mass, pos, vel
#OPT  += -DTHIRD_HALO_ONLY 

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
CC      	=  gcc
OPTIMIZE	= -O3 -m64 -fopenmp -mtune=native 
GSL_INCL 	= $(CPPFLAGS)
GSL_LIBS	= -L/Users/julius/Devel/lib
endif

ifeq ($(SYSTYPE),MPA64)
CC      	=  gcc
OPTIMIZE	= -O3 -Wall -g 
GSL_INCL 	= -I/afs/mpa/home/jdonnert/Libs/@sys/include
GSL_LIBS	= -L/afs/mpa/home/jdonnert/Libs/@sys/lib
FFTW_LIBS	= 
FFTW_INCL 	=
endif

ifeq ($(SYSTYPE),coma.msi.umn.edu)
CC      	= icc
OPTIMIZE	= -O2 -Wall -g -xhost 
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
