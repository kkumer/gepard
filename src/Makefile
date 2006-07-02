# $Id$
# Makefile for gepard Fortran and C sources
# 2006-04-13 

# Interesting targets:
# 	radcorr  -  program for producing Fig. 1 in letter
# 	scaledep -  program for producing Fig. 2 in letter
# 	gepard.exe - Mathematica-installable package
#
# For debugging/profiling call like this:
#  make DEBUG='-g -pg' and change -ladacf -> -ladacf_prof


# Where and how gepard.exe should be installed
OWNER = $(USER) 
GROUP = $(USER)
INSTALL = install
BINDIR = $(HOME)/.mma/Applications/gepard.exe/Linux/


# Compilers
MATHLINKDIR=/usr/local/Wolfram/Mathematica/5.2/AddOns/MathLink/DeveloperKit/Linux/CompilerAdditions
MPREP = $(MATHLINKDIR)/mprep

# Compiler options
# DEBUG = -g -pg
FFLAGS = -O2 $(DEBUG)
CFLAGS = -O2 $(DEBUG)

# Libraries and include files that are used
#LIBS =  -L/home/kkumer/local/lib -ladacf -lslatec -L$(MATHLINKDIR) -lML -lpthread
LIBS =  -L$(HOME)/local/lib -ladacf -lslatec 
MLLIBS = -L$(MATHLINKDIR) -lML -lpthread
INCS =  -I$(MATHLINKDIR)

# Directories for storing object files
FODIR = Fobj
ODIR = obj

# common Fortran object files
FOBJS = Fcommon.o Fbeta.o Fspecial.o Fdcomplex.o Fastrong.o Flambda.o\
Fprojectors.o Frnnlo.o Ferfunc.o Fevol.o Fcdvcs.o Fmsbar.o Fparwav.o Fcff.o\
Fdvcs.o Ff2.o
_FOBJS = $(patsubst %,$(FODIR)/%,$(FOBJS))

# rule for compiling Fortran files
$(FODIR)/F%.o: %.f 
	$(FC) -c -o $@ $< $(FFLAGS)

# common C object files
OBJS =  ALLtm.o main.o mlink.o common.o astrong.o lambda.o projectors.o\
rnnlo.o erfunc.o cdvcs.o evol.o msbar.o 
_OBJS = $(patsubst %,$(ODIR)/%,$(OBJS))

# rule for compiling C files
$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(INCS)

# Template files
TMS = astrong.tm cdvcs.tm common.tm erfunc.tm evol.tm lambda.tm\
msbar.tm projectors.tm rnnlo.tm

# targets

all: radcorr scaledep gepard.exe fit test

ALLtm.c: $(TMS)
	cat $(TMS) > ALL.tm
	$(MPREP) ALL.tm -o $@

radcorr: $(_FOBJS) $(FODIR)/Fradcorr.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) 

scaledep: $(_FOBJS) $(FODIR)/Fscaledep.o $(FODIR)/Fderiv.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) 

gepard.exe: $(_FOBJS) $(_OBJS) 
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) $(MLLIBS)

fit: $(_FOBJS) $(FODIR)/Ffit.o $(ODIR)/flushout.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) -lpacklib -lkernlib

test: $(_FOBJS) $(FODIR)/Ftest.o 
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS)

install: gepard.exe
	$(INSTALL) -o $(OWNER) -g $(GROUP) -m 755 $< $(BINDIR)

.PHONY: clean
clean:
	-rm -f $(FODIR)/F*.o $(ODIR)/*.o FIG*DAT ALL*tm*  radcorr scaledep gepard.exe
