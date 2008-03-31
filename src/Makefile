# Makefile for gepard Fortran and C sources

# Interesting targets:
# 	test       -  tests both DVCS and DIS routines
# 	radcorr    -  program for producing Fig. 1 in letter
# 	scaledep   -  program for producing Fig. 2 in letter
# 	fit        -  fitting GPD ansatz to DVCS and DIS data
#
# To remove PGPLOT (plotting) compile with
# 	make NOPGPLOT=1 fit 
#      or
#   make fit_nopgplot
#


# Compiler options
FFLAGS =  -I. -O 
CFLAGS =  -O2
CPPFLAGS = -D NOPGPLOT

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
Finitc.o Fevolc.o Fgetmbgpd.o Fhj.o Fsplice.o Fcff.o Fderiv.o Fbca.o
_FOBJS = $(patsubst %,$(FODIR)/%,$(FOBJS))


# rule for compiling Fortran .f  files
# (prefix "F" is used to avoid possible future clashes with
#  C-wrappers with the same name)
$(FODIR)/F%.o: %.f  header.f
	$(FC) -c -o $@ $< $(FFLAGS)
$(FODIR)/F%.o: %.F  header.f
	$(FC) -c -o $@ $< $(FFLAGS) $(CPPFLAGS)

# C object files
OBJS = flushout.o
_OBJS = $(patsubst %,$(ODIR)/%,$(OBJS))

# rule for compiling C files
$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(INCS)

# All object files together
export ALLOBJS = $(_OBJS) $(_ADACFOBJS) $(_FOBJS)

# Mathlink C object files
MMAOBJS = fit.tm.o fit.o 
_MMAOBJS = $(patsubst %,$(MMAODIR)/%,$(MMAOBJS))

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

# test and fit targets

$(TESTTARGETS): %: $(ALLOBJS) $(FODIR)/F%.o
	$(FC) -o ../Tests/$@ $^ $(FFLAGS)

#fit_nopgplot: $(ALLOBJS) $(FODIR)/Ffit.o
	#$(FC) -o $@ $^ $(FFLAGS) $(CERNLIBS)

# main fit target  (a variant withouth PGPLOT)
fit: CPPFLAGS :=
fit_nopgplot: ALLPGPLOTLIBS :=
fit fit_nopgplot:  $(ALLOBJS) $(FODIR)/Ffcn.o $(FODIR)/Ffit.o
	$(FC) -o ../fits/$@ $^ $(FFLAGS) $(CERNLIBS) $(ALLPGPLOTLIBS)

# targets for processing/presenting GPDs resulting from fit
$(FITTARGETS): %: $(ALLOBJS) $(EXFODIR)/%.o
	$(FC) -o ../fits/$@ $^ $(FFLAGS)

# mma targets

gepard.exe: CPPFLAGS := -D MMA
gepard.exe: $(ALLOBJS) $(FODIR)/Ffcn.o $(_MMAOBJS)  $(FODIR)/Fmninterface.o
	$(FC) -o ../fits/$@/$(SYS)/$@ $^ $(FFLAGS) $(CERNLIBS) -L$(MLLIBDIR) -l$(MLLIB) $(ALLPGPLOTLIBS) $(MLEXTRA)

# targets for examples
#
$(EXTARGETS): %: $(ALLOBJS) $(EXFODIR)/%.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)


.PHONY: clean
clean:
	-rm -f $(ADACFODIR)/*.o $(FODIR)/F*.o $(ODIR)/*.o $(EXFODIR)/*.o $(MMAODIR)/*
	-rm -f $(patsubst %,../Tests/%,$(TESTTARGETS))
	-rm -f $(patsubst %,../Tests/%.exe,$(TESTTARGETS))
	-rm -f $(patsubst %,../fits/%,$(FITTARGETS))
	-rm -f $(patsubst %,../fits/%.exe,$(FITTARGETS))
	-rm -f ../fits/fit ../fits/fit.exe ../fits/fit_nopgplot ../fits/fit_nopgplot.exe
	-rm -f $(patsubst %,../ex/%,$(EXTARGETS))
	-rm -f $(patsubst %,../ex/%.exe,$(EXTARGETS))
	-rm -f ../fits/gepard.exe/$(SYS)/gepard.exe
