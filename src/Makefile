# Makefile for gepard Fortran and C sources

# Subdirectories for storing example source files
EXDIR = ex

# Directories for storing object files
ADACFODIR = adacf/obj
FODIR = Fobj
ODIR = obj
MMAODIR = Mobj
EXFODIR = ex/obj

ifndef NOPGPLOT
  USEPGPLOT := -DPGPLOT
endif

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

# Fortran object files needed for DVEM
FDVEMOBJS = Fhyp2F1.o Fgauss.o Fintegra.o Fint2F1.o Fbeta.o Fspecial.o Fcdvem.o 
_FDVEMOBJS = $(patsubst %,$(FODIR)/%,$(FDVEMOBJS))

# rule for compiling Fortran .F  files
# (prefix "F" in "F%.o" is used to avoid possible future clashes with
#  C-wrappers with the same name)
$(FODIR)/Fhyp2F1.o: hyp2F1.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(CPPFLAGS)

$(FODIR)/F%.o: %.F  header.F
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

MMADVEMOBJS = dvem.tm.o dvem.o 
_MMADVEMOBJS = $(patsubst %,$(MMAODIR)/%,$(MMADVEMOBJS))

# rules for compiling MathLink C files
$(MMAODIR)/%.tm.c: %.tm
	$(MPREP) -o $@ $< 
$(MMAODIR)/%.tm.o: $(MMAODIR)/%.tm.c $(MMAODIR)/%.o
	$(CC) -c -o $@ $< -I$(MLINCDIR) $(CFLAGS) $(MLEXTRA)
$(MMAODIR)/%.o: %.c  header.h
	$(CC) -c -o $@ $< -I$(MLINCDIR) $(CFLAGS)


# rule for compiling example Fortran programs
$(EXFODIR)/%.o: $(EXDIR)/%.F header.F
	$(FC) -c -o $@ $< $(FFLAGS)

# Fortran prep files (needed for f2py which ignores preprocessing)
FPREPDIR = Fprep
FPREPS = init.f beta.f special.f dcomplex.f astrong.f lambda.f\
projectors.f rnnlo.f erfunc.f evol.f cdvcs.f msbar.f \
dvcs.f f2.f error.f gauss.f nd.f bigcns.f integra.f\
initc.f evolc.f getmbgpd.f hj.f splice.f cff.f 
_FPREPS = $(patsubst %,$(FPREPDIR)/%,$(FPREPS))

$(FPREPDIR)/%.f: %.F  header.F
	$(FC) -E $<  > $@

# test and fit targets

$(TESTTARGETS): %: $(ALLOBJS) $(FODIR)/F%.o
	$(FC) -o ../Tests/$@ $^ $(FFLAGS)

# main fit target  (a variant withouth PGPLOT)
fit: CPPFLAGS += $(USEPGPLOT)
fit: $(ALLOBJS) $(FODIR)/Ffcn.o $(FODIR)/Ffit.o
	$(FC) -o ../fits/$@ $^ $(FFLAGS) $(CERNLIBS) $(ALLPGPLOTLIBS)

# targets for processing/presenting GPDs resulting from fit
$(FITTARGETS): %: $(ALLOBJS) $(EXFODIR)/%.o
	$(FC) -o ../fits/$@ $^ $(FFLAGS) $(ALLPGPLOTLIBS)

# mma targets

gepard.exe: CPPFLAGS += -DMMA  $(USEPGPLOT)
gepard.exe: $(ALLOBJS) $(FODIR)/Ffcn.o $(_MMAOBJS)  $(FODIR)/Fmninterface.o
	$(FC) -o ../fits/$@/$(SYS)/$@ $^ $(FFLAGS) $(CERNLIBS) -L$(MLLIBDIR) -l$(MLLIB) $(ALLPGPLOTLIBS) $(MLEXTRA)

int2f1.exe: $(FODIR)/Fhyp2F1.o $(FODIR)/Fgauss.o $(FODIR)/Fintegra.o  \
             $(FODIR)/Fint2F1F.o $(_MMADVEMOBJS)
	$(FC) -o ../fits/$@/$(SYS)/$@ $^ $(FFLAGS)  -L$(MLLIBDIR) -l$(MLLIB)  $(MLEXTRA)

dvem.exe: $(ADACFODIR)/Harmonic.o $(_FDVEMOBJS) $(_MMADVEMOBJS)
	$(FC) -o ../fits/$@/$(SYS)/$@ $^ $(FFLAGS)  -L$(MLLIBDIR) -l$(MLLIB)  $(MLEXTRA)
# pygepard target

pygepard.so: $(_FPREPS)
	f2py --fcompiler=$(FC) -c -m pygepard $^ adacf/*f
	mv pygepard.so ../fits/

# targets for examples
#
$(EXTARGETS): %: $(ALLOBJS) $(EXFODIR)/%.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS)


.PHONY: clean
clean:
	-rm -f $(ADACFODIR)/*.o $(FODIR)/F*.o $(ODIR)/*.o $(EXFODIR)/*.o $(MMAODIR)/*
	-rm -f $(FPREPDIR)/*.f
	-rm -f $(patsubst %,../Tests/%,$(TESTTARGETS))
	-rm -f $(patsubst %,../Tests/%.exe,$(TESTTARGETS))
	-rm -f $(patsubst %,../fits/%,$(FITTARGETS))
	-rm -f $(patsubst %,../fits/%.exe,$(FITTARGETS))
	-rm -f ../fits/fit ../fits/fit.exe  ../fits/pygepard.so
	-rm -f $(patsubst %,../ex/%,$(EXTARGETS))
	-rm -f $(patsubst %,../ex/%.exe,$(EXTARGETS))
	-rm -f ../fits/gepard.exe/$(SYS)/gepard.exe
	-rm -f ../fits/int2f1.exe/$(SYS)/int2f1.exe
