SHELL = /bin/bash

## OPTIONS  ##
OPT 	+= -DNFWC_DUFFY08	# alternate fit to concentr. param

OPT += -DBETA=0.54

#OPT     += -DPARABOLA       # merge in a parabola
OPT	+= -DCOMET			# merge like a comet, ball+tail (recommended)

#OPT	+= -DDOUBLE_BETA_COOL_CORES

#OPT 	+= -DGIVEPARAMS		# more merger parameters in .par file

OPT	+= -DNO_RCUT_IN_T		# set Rcut very large

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
OPTIMIZE	= -m64  -O0 #-xhost -fast
GSL_INCL 	= $(CPPFLAGS)
GSL_LIBS	= -L/Users/jdonnert/Dev/lib
endif

ifeq ($(SYSTYPE),MSI)
CC      	= icc
OPTIMIZE	= -Wall -g  -O3 -xhost
GSL_INCL 	=
GSL_LIBS	=
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
