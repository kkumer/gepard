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
FFLAGS =  -O
CFLAGS =  -O2

# Subdirectories for storing source files
EXDIR = ex

# Directories for storing object files
FODIR = Fobj
ODIR = obj
ADACFODIR = adacf/obj
EXFODIR = ex/obj

# Fortran object files
FOBJS = Finit.o Fbeta.o Fspecial.o Fdcomplex.o Fastrong.o Flambda.o\
Fprojectors.o Frnnlo.o Ferfunc.o Fevol.o Fcdvcs.o Fmsbar.o Fparwav.o Fcff.o\
Fdvcs.o Ff2.o Freadpar.o Ferror.o Fgauss.o Fnd.o
_FOBJS = $(patsubst %,$(FODIR)/%,$(FOBJS))

# rule for compiling example Fortran programs
$(EXFODIR)/%.o: $(EXDIR)/%.f
	$(FC) -c -o $@ $< $(FFLAGS)

# rule for compiling fit.F
$(FODIR)/Ffit.o: fit.F
	$(FC) -c -o $@ $< $(FFLAGS) $(CPPFLAGS)
$(FODIR)/Ffit_nopgplot.o: fit.F
	$(FC) -c -o $@ $< $(FFLAGS) '-D NOPGPLOT'

# rule for compiling Fortran .f  files
# (prefix "F" is used to avoid possible future clashes with
#  C-wrappers with the same name)
$(FODIR)/F%.o: %.f
	$(FC) -c -o $@ $< $(FFLAGS)

# "adacf" Fortran object files
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

# C object files
OBJS = flushout.o
_OBJS = $(patsubst %,$(ODIR)/%,$(OBJS))

# rule for compiling C files
$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(INCS)

# targets

auxtest: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Fauxtest.o
	$(FC) -o $@ $^ $(FFLAGS)

test: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Ftest.o
	$(FC) -o $@ $^ $(FFLAGS)

radcorr: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Fradcorr.o
	$(FC) -o $@ $^ $(FFLAGS)

scaledep: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Fscaledep.o $(FODIR)/Fderiv.o
	$(FC) -o $@ $^ $(FFLAGS)

fit: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Ffit.o $(ODIR)/flushout.o
	$(FC) -o $@ $^ $(FFLAGS) $(CERNLIBS) $(ALLPGPLOTLIBS)

fit_nopgplot: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Ffit_nopgplot.o $(ODIR)/flushout.o
	$(FC) -o $@ $^ $(FFLAGS) $(CERNLIBS)

houches: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Fhouches.o 
	$(FC) -o $@ $^ $(FFLAGS)

accuracy: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Faccuracy.o 
	$(FC) -o $@ $^ $(FFLAGS)

evoldep: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Fevoldep.o
	$(FC) -o $@ $^ $(FFLAGS)

auxns: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Fauxns.o
	$(FC) -o $@ $^ $(FFLAGS)

aux: $(_FOBJS) $(_ADACFOBJS) $(FODIR)/Faux.o
	$(FC) -o $@ $^ $(FFLAGS)

# targets for examples
radNLONS: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/radNLONS.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

radNLO: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/radNLO.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

radQ: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/radQ.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

radQNS: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/radQNS.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

radNNLO: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/radNNLO.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

scalesNS: $(_FOBJS) $(_ADACFOBJS) $(EXFODIR)/scalesNS.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

.PHONY: clean
clean:
	-rm -f $(SRCTARGETS)
	-rm -f $(patsubst %,%.exe,$(SRCTARGETS))
	-rm -f $(ADACFODIR)/*.o $(FODIR)/F*.o $(ODIR)/*.o $(EXFODIR)/*.o 
	-rm -f $(patsubst %,../ex/%,$(EXTARGETS))
	-rm -f $(patsubst %,../ex/%.exe,$(EXTARGETS))
