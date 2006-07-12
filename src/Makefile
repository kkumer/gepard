# Makefile for gepard Fortran and C sources

# Interesting targets:
# 	radcorr  -  program for producing Fig. 1 in letter
# 	scaledep -  program for producing Fig. 2 in letter
# 	test     -  tests both DVCS and DIS routines
# 	fit      -  fitting GPD ansatz to DVCS and DIS data
#
# For debugging/profiling call like this:
#  make DEBUG='-g -pg' and change -ladacf -> -ladacf_prof
# For optimization call like this:
#  make DEBUG='-O2'


# Compiler options
# DEBUG = -g -pg
FFLAGS =  $(DEBUG)
CFLAGS =  $(DEBUG)

# Libraries and include files that are used
#LIBS =  -L/home/kkumer/local/lib -ladacf -L$(MATHLINKDIR) -lML  -lpthread
LIBS =  -L$(HOME)/local/lib -ladacf 

# Directories for storing object files
FODIR = Fobj
ODIR = obj

# common Fortran object files
FOBJS = Finit.o Fbeta.o Fspecial.o Fdcomplex.o Fastrong.o Flambda.o\
Fprojectors.o Frnnlo.o Ferfunc.o Fevol.o Fcdvcs.o Fmsbar.o Fparwav.o Fcff.o\
Fdvcs.o Ff2.o Freadpar.o Ferror.o
_FOBJS = $(patsubst %,$(FODIR)/%,$(FOBJS))

# rule for compiling Fortran files
$(FODIR)/F%.o: %.f 
	$(FC) -c -o $@ $< $(FFLAGS)

# common C object files
OBJS = 
_OBJS = $(patsubst %,$(ODIR)/%,$(OBJS))

# rule for compiling C files
$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(INCS)

# targets

all: radcorr scaledep fit test auxtest

ALLtm.c: $(TMS)
	cat $(TMS) > ALL.tm
	$(MPREP) ALL.tm -o $@

radcorr: $(_FOBJS) $(FODIR)/Fradcorr.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) 

scaledep: $(_FOBJS) $(FODIR)/Fscaledep.o $(FODIR)/Fderiv.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) 

fit: $(_FOBJS) $(FODIR)/Ffit.o $(ODIR)/flushout.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) -lpacklib -lkernlib

auxtest: $(_FOBJS) $(FODIR)/Fauxtest.o 
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS)

test: $(_FOBJS) $(FODIR)/Ftest.o 
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS)

.PHONY: clean
clean:
	-rm -f $(FODIR)/F*.o $(ODIR)/*.o FIG*DAT radcorr scaledep fit test auxtest
