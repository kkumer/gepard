# Makefile for gepard Fortran and C sources

# Interesting targets:
# 	radcorr    -  program for producing Fig. 1 in letter
# 	scaledep   -  program for producing Fig. 2 in letter
# 	test       -  tests both DVCS and DIS routines
# 	fit        -  fitting GPD ansatz to DVCS and DIS data
#
# To remove PGPLOT (plotting) compile with
# 	make NOPGPLOT=1 fit 
#      or
#   make fit_nopgplot
#

ifndef NOPGPLOT
ALLPGPLOTLIBS = $(PGPLOTLIBS) $(X11LIBS)
else
CPPFLAGS = '-D NOPGPLOT'
endif

# Compiler options
FFLAGS =  -I. -O 
CFLAGS =  -O2

# Subdirectories for storing example source files
EXDIR = ex

# Directories for storing object files
ADACFODIR = adacf/obj
FODIR = Fobj
ODIR = obj
MMAODIR = Mobj
EXFODIR = ex/obj

# "adacf" package Fortran object files
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


# Fortran object files
FOBJS = Finit.o Fbeta.o Fspecial.o Fdcomplex.o Fastrong.o Flambda.o\
Fprojectors.o Frnnlo.o Ferfunc.o Fevol.o Fcdvcs.o Fmsbar.o \
Fdvcs.o Ff2.o Freadpar.o Ferror.o Fgauss.o Fnd.o Fbigcns.o Fintegra.o Fxspace.o\
Finitc.o Fevolc.o Fgetmbgpd.o Fhj.o Fsplice.o
_FOBJS = $(patsubst %,$(FODIR)/%,$(FOBJS))

# rule for compiling fcn.F
$(FODIR)/Ffcn.o: fcn.F header.f
	$(FC) -c -o $@ $< $(FFLAGS) $(CPPFLAGS)
$(FODIR)/Ffcn_nopgplot.o: fcn.F header.f
	$(FC) -c -o $@ $< $(FFLAGS) '-D NOPGPLOT'

# rule for compiling *.F
$(FODIR)/Fcff.o: cff.F header.f
	$(FC) -c -o $@ $< $(FFLAGS) $(CPPFLAGS)

# rule for compiling Fortran .f  files
# (prefix "F" is used to avoid possible future clashes with
#  C-wrappers with the same name)
$(FODIR)/F%.o: %.f  header.f
	$(FC) -c -o $@ $< $(FFLAGS)

# C object files
OBJS = flushout.o atest.o
_OBJS = $(patsubst %,$(ODIR)/%,$(OBJS))

# rule for compiling C files
$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(INCS)

# rules for compiling MathLink C files
$(MMAODIR)/%.tm.c: %.tm
	$(MPREP) -o $@ $< 
$(MMAODIR)/%.tm.o: $(MMAODIR)/%.tm.c $(MMAODIR)/%.o
	$(CC) -c -o $@ $< -I$(MLINCDIR) $(CFLAGS) $(MLEXTRA)
$(MMAODIR)/%.o: %.c  header.h
	$(CC) -c -o $@ $< -I$(MLINCDIR) $(CFLAGS)


# rule for compiling example Fortran programs
$(EXFODIR)/%.o: $(EXDIR)/%.f header.f
	$(FC) -c -o $@ $< $(FFLAGS)

# mma targets

gepard.exe: FFLAGS := $(FFLAGS) '-D MMA'
gepard.exe: $(_FOBJS) $(_ADACFOBJS) $(MMAODIR)/fit.tm.o $(MMAODIR)/fit.o $(ODIR)/flushout.o $(FODIR)/Ffcn.o $(FODIR)/Fmninterface.o $(FODIR)/Fcff.o
	$(FC) -o ../fits/$@/$(SYS)/$@ $^ $(FFLAGS) $(CERNLIBS) -L$(MLLIBDIR) -l$(MLLIB) $(ALLPGPLOTLIBS) $(MLEXTRA)

# targets
atest: $(_FOBJS) $(_ADACFOBJS) $(ODIR)/atest.o $(FODIR)/Fcff.o
	$(FC) -o $@ $^ $(FFLAGS)

auxtest: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Fauxtest.o $(FODIR)/Fcff.o
	$(FC) -o $@ $^ $(FFLAGS)

test: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Ftest.o $(FODIR)/Fcff.o
	$(FC) -o $@ $^ $(FFLAGS)

radcorr: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Fradcorr.o $(FODIR)/Fcff.o
	$(FC) -o $@ $^ $(FFLAGS)

scaledep: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Fscaledep.o $(FODIR)/Fderiv.o $(FODIR)/Fcff.o
	$(FC) -o $@ $^ $(FFLAGS)

fit: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Ffit.o $(FODIR)/Ffcn.o $(ODIR)/flushout.o $(FODIR)/Fcff.o
	$(FC) -o $@ $^ $(FFLAGS) $(CERNLIBS) $(ALLPGPLOTLIBS)

fit_nopgplot: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Ffit.o $(FODIR)/Ffcn_nopgplot.o $(ODIR)/flushout.o $(FODIR)/Fcff.o
	$(FC) -o $@ $^ $(FFLAGS) $(CERNLIBS)

houches: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Fhouches.o $(FODIR)/Fcff.o
	$(FC) -o $@ $^ $(FFLAGS)

accuracy: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Faccuracy.o $(FODIR)/Fcff.o
	$(FC) -o $@ $^ $(FFLAGS)

# targets for examples
#
auxns: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/auxns.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

auxsi: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/auxsi.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

anatomyNS: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/anatomyNS.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

anatomy: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/anatomy.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

radNLONS: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/radNLONS.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

radNLO: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/radNLO.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

radQ: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/radQ.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

evolutNS: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/evolutNS.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

evolut: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/evolut.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

radNNLONS: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/radNNLONS.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

radNNLO: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/radNNLO.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

scalesNS: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/scalesNS.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

scales: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/scales.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

scalesNNLO: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/scalesNNLO.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

fitres: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/fitres.o
	$(FC) -o ../fits/$@ $^ $(FFLAGS)

fitpdfs: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/fitpdfs.o
	$(FC) -o ../fits/$@ $^ $(FFLAGS)

slope: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/slope.o $(FODIR)/Fderiv.o
	$(FC) -o ../fits/$@ $^ $(FFLAGS)

.PHONY: clean
clean:
	-rm -f $(ADACFODIR)/*.o $(FODIR)/F*.o $(ODIR)/*.o $(EXFODIR)/*.o $(MMAODIR)/*
	-rm -f $(SRCTARGETS)
	-rm -f $(patsubst %,%.exe,$(SRCTARGETS))
	-rm -f ../fits/gepard.exe/$(SYS)/gepard.exe
	-rm -f $(patsubst %,../ex/%,$(EXTARGETS))
	-rm -f $(patsubst %,../ex/%.exe,$(EXTARGETS))
