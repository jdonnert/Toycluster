# This is Toycluster, a code that generates an artificial cluster merger.
# Clusters consist of a DM halo defined by its density profile and a gaseous
# ICM defined by a beta model. Under the assumption of hydrostatic equillibrium
# all other quantities follow. (Donnert 2014, Donnert et al in prep.)

SHELL = /bin/bash

OPT += -DRCUT_R200_RATIO=1.2 # values >1 make cluster unstable, use with care

#OPT += -DPARABOLA       # merge in a parabola
OPT	+= -DCOMET			 # merge like a comet, ball+tail (recommended)
						 # if nothing is selected, merge as ball with R_Sample

#OPT 	+= -DGIVEPARAMS		 # set beta models in parameter file

#OPT	+= -DNO_RCUT_IN_T	 # set Rcut very large in U calculation

#
OPT += -DSUBSTRUCTURE		 # add a population of galaxy-like subhalos
OPT += -DSUBHOST=1			 # host subhalos in this cluster
OPT += -DSLOW_SUBSTRUCTURE	 # put subhalos on Hernquist orbits
OPT += -DREPORTSUBHALOS	 # print info about all subhaloes

#OPT += -DADD_THIRD_SUBHALO  # manually set the first subhalo mass, pos, vel
#OPT  += -DTHIRD_HALO_ONLY

#OPT += -DSPH_CUBIC_SPLINE 	 # for use with Gadget2

#OPT	+= -DDOUBLE_BETA_COOL_CORES # cool cores as double beta model

OPT 	+= -DNFWC_DUFFY08	 # alternate fit to concentr. param

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
OPTIMIZE	= -fast -g
GSL_INCL 	= $(CPPFLAGS)
GSL_LIBS	= -L/Users/jdonnert/Dev/lib
endif

ifeq ($(SYSTYPE),MSI)
CC      	= icc
OPTIMIZE	= -Wall -g  -O3 -xhost -ipo4
GSL_INCL 	= -I/home/jonestw/donne219/Libs/$(shell hostname)/include
GSL_LIBS	= -L/home/jonestw/donne219/Libs/$(shell hostname)/lib
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

EXEC = Toycluster_sub

## FILES ##

SRCDIR	= src/
 
SRCFILES := ${shell find $(SRCDIR) -name \*.c -print} # all .c files in SRCDIR
OBJFILES = $(SRCFILES:.c=.o)

INCLFILES := ${shell find src -name \*.h -print} # all .h files in SRCDIR
INCLFILES += Makefile

CFLAGS 	= -std=c99 -fopenmp $(OPTIMIZE) $(OPT) $(GSL_INCL) $(FFTW_INCL)

LINK	= $(GSL_LIBS) -lm -lgsl -lgslcblas 

## RULES ## 

%.o : %.c
	@echo [CC] $@
	@$(CC) $(CFLAGS)  -o $@ -c $<

$(EXEC)	: $(OBJFILES)
	@echo SYSTYPE=$(SYSTYPE)
	$(CC) $(CFLAGS) $(OBJFILES) $(LINK) -o $(EXEC)
	@ctags -w $(SRCFILES) $(INCLFILES)

$(OBJFILES)	: $(INCLFILES) $(SRCFILES)

clean	: 
	rm -f  $(OBJFILES) $(EXEC)
