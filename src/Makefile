# Makefile for gepard Fortran and C sources

# Interesting targets:
# 	radcorr  -  program for producing Fig. 1 in letter
# 	scaledep -  program for producing Fig. 2 in letter
# 	test     -  tests both DVCS and DIS routines
# 	fit      -  fitting GPD ansatz to DVCS and DIS data
#
# For debugging/profiling call like this:
#  make DEBUG='-g -pg'
# For optimization call like this:
#  make DEBUG='-O2'


# Compiler options
FFLAGS =  $(DEBUG)

# Where are CERNLIB's kernlib and packlib
LIBS =  -L$(HOME)/local/lib

# Directories for storing object files
FODIR = Fobj
ODIR = obj
ADACFODIR = adacf/obj

# common Fortran object files
FOBJS = Finit.o Fbeta.o Fspecial.o Fdcomplex.o Fastrong.o Flambda.o\
Fprojectors.o Frnnlo.o Ferfunc.o Fevol.o Fcdvcs.o Fmsbar.o Fparwav.o Fcff.o\
Fdvcs.o Ff2.o Freadpar.o Ferror.o
_FOBJS = $(patsubst %,$(FODIR)/%,$(FOBJS))

# rule for compiling Fortran files
$(FODIR)/F%.o: %.f
	$(FC) -c -o $@ $< $(FFLAGS)

# common "adacf" Fortran object files
ADACFOBJS = Harmonic.o Mellin.o WgammaVNSP0F.o WgammaVNSP1F.o\
WgammaVNSP2F.o WgammaVQQ0F.o WgammaVQQ1F.o WgammaVQQ2F.o WgammaVQG0F.o\
WgammaVQG1F.o WgammaVQG2F.o WgammaVGQ0F.o WgammaVGQ1F.o WgammaVGQ2F.o\
WgammaVGG0F.o WgammaVGG1F.o WgammaVGG2F.o WcVF2NSP0F.o WcVF2NSP1F.o\
WcVF2NSP2F.o WcVFLNSP0F.o WcVFLNSP1F.o WcVFLNSP2F.o WcVF2Q0F.o WcVF2Q1F.o\
WcVF2Q2F.o WcVFLQ0F.o WcVFLQ1F.o WcVFLQ2F.o WcVF2G0F.o WcVF2G1F.o\
WcVF2G2F.o WcVFLG0F.o WcVFLG1F.o WcVFLG2F.o
_ADACFOBJS = $(patsubst %,$(ADACFODIR)/%,$(ADACFOBJS))

# rule for compiling "adacf" Fortran files
$(ADACFODIR)/%.o: adacf/%.f 
	$(FC) -c -o $@ $< $(FFLAGS)

# common C object files
OBJS = flushout.o
_OBJS = $(patsubst %,$(ODIR)/%,$(OBJS))

# rule for compiling C files
$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(INCS)

# targets

all: auxtest test radcorr scaledep fit

auxtest: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Fauxtest.o 
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS)

test: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Ftest.o 
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS)

radcorr: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Fradcorr.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) 

scaledep: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Fscaledep.o $(FODIR)/Fderiv.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) 

fit: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Ffit.o $(ODIR)/flushout.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) -lpacklib -lkernlib


.PHONY: clean
clean:
	-rm -f $(ADACFODIR)/*.o $(FODIR)/F*.o $(ODIR)/*.o\
		   	auxtest test radcorr scaledep fit
