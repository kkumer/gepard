# Makefile for gepard Fortran and C sources
PYSUFFIX = $(shell python-config --extension-suffix)

# Subdirectories for storing example source files
TESTDIR = test
EXDIR = ex

# Directories for storing object files
ADACFODIR = adacf/obj
FODIR = Fobj
ODIR = obj
MMAODIR = Mobj
TESTFODIR = test/obj
EXFODIR = ex/obj


# "adacf" package Fortran object files
ADACFOBJS = Harmonic.o WgammaVNSP0F.o WgammaVNSP1F.o WgammaVNSM1F.o\
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
Fprojectors.o Frnnlo.o Ferfunc.o Fevol.o Fwccsbar.o Fwcmsbar.o \
Fdvcs.o Ff2.o Ferror.o Fgauss.o Fnd.o Fintegra.o Freadpar.o\
Fevolc.o Fgetmbgpd.o Fhj.o Fej.o Fcff.o Ftff.o Fderiv.o\
Fwcdvmp.o Fc1dvmp.o
_FOBJS = $(patsubst %,$(FODIR)/%,$(FOBJS))


# rule for compiling Fortran .F  files
# (prefix "F" in "F%.o" is used to avoid possible future clashes with
#  C-wrappers with the same name)
#
$(FODIR)/F%.o: %.F  header.F
	$(FC) -c -o $@ $< $(FFLAGS) $(CPPFLAGS)

# rule for compiling C files
$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(INCS)

# All object files together
export ALLOBJS = $(_ADACFOBJS) $(_FOBJS)

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

# rule for compiling test Fortran programs
$(TESTFODIR)/%.o: $(TESTDIR)/%.F header.F
	$(FC) -c -o $@ $< $(FFLAGS)

# rule for compiling example Fortran programs
$(EXFODIR)/%.o: $(EXDIR)/%.F header.F
	$(FC) -c -o $@ $< $(FFLAGS)

# Fortran prep files (needed for f2py which ignores preprocessing)
FPREPDIR = Fprep
FPREPS = init.f beta.f special.f dcomplex.f astrong.f lambda.f\
projectors.f rnnlo.f erfunc.f evol.f wccsbar.f wcmsbar.f \
dvcs.f f2.f error.f gauss.f nd.f integra.f readpar.f\
evolc.f getmbgpd.f hj.f ej.f cff.f\
tff.f wcdvmp.f c1dvmp.f
_FPREPS = $(patsubst %,$(FPREPDIR)/%,$(FPREPS))

$(FPREPDIR)/%.f: %.F  header.F
	$(FC) -E $<  > $@


# main fit target
fit: $(ALLOBJS) $(FODIR)/Ffcn.o $(FODIR)/Ffit.o
	$(FC) -o ../fits/$@ $^ $(FFLAGS) $(CERNLIBS)

# targets for processing/presenting GPDs resulting from fit
$(FITTARGETS): %: $(ALLOBJS) $(EXFODIR)/%.o
	$(FC) -o ../fits/$@ $^ $(FFLAGS)

# mma targets

gepard.exe: CPPFLAGS += -DMMA
gepard.exe: $(ALLOBJS) $(FODIR)/Ffcn.o $(_MMAOBJS)  $(FODIR)/Fmninterface.o
	$(FC) -o ../fits/$@/$(SYS)/$@ $^ $(FFLAGS) $(CERNLIBS) -L$(MLLIBDIR) -l$(MLLIB) $(MLEXTRA)


# pygepard targets
#
pygepard: pygepard.so

pygepards: pygepard.so pygepard2.so pygepard3.so pygepard4.so pygepard5.so pygepard6.so pygepard7.so

pygepard.so: $(_FPREPS)
	f2py --fcompiler=gnu95  -c -m pygepard $^ adacf/*f -lgomp 
	mv pygepard$(PYSUFFIX) ../pype/pygepard.so

pygepard%.so: $(_FPREPS)
	f2py --quiet --fcompiler=gnu95 -c -m pygepard$* $^ adacf/*f -lgomp
	mv pygepard$*$(PYSUFFIX) ../pype/pygepard$*.so


# targets for examples
#
$(EXTARGETS): %: $(ALLOBJS) $(EXFODIR)/%.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)

$(TESTTARGETS): %: $(ALLOBJS) $(TESTFODIR)/%.o
	$(FC) -o ../test/$@ $^ $(FFLAGS)


.PHONY: clean
clean:
	-rm -f $(ADACFODIR)/*.o $(FODIR)/F*.o $(ODIR)/*.o $(EXFODIR)/*.o $(MMAODIR)/*
	-rm -f *.s
	-rm -f $(FPREPDIR)/*.f
	-rm -f $(patsubst %,../test/%,$(TESTTARGETS))
	-rm -f $(patsubst %,../test/%.exe,$(TESTTARGETS))
	-rm -f $(patsubst %,../fits/%,$(FITTARGETS))
	-rm -f $(patsubst %,../fits/%.exe,$(FITTARGETS))
	-rm -f ../fits/fit ../fits/fit.exe  ../fits/pygepard.so
	-rm -f $(patsubst %,../ex/%,$(EXTARGETS))
	-rm -f $(patsubst %,../ex/%.exe,$(EXTARGETS))
	-rm -f ../ifaces/gepard.exe/$(SYS)/gepard.exe
	-rm -f ../pype/pygepar*.so
