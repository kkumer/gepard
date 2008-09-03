# Makefile for gepard Fortran and C sources

# Subdirectories for storing example source files
EXDIR = ex

# Directories for storing object files
FODIR = Fobj
ODIR = obj
MMAODIR = Mobj
EXFODIR = ex/obj

ifndef NOPGPLOT
  USEPGPLOT := -D PGPLOT
endif

# Fortran object files
FOBJS = Finit.o Fbeta.o Fspecial.o Fdcomplex.o Fastrong.o Flambda.o\
Fprojectors.o Frnnlo.o Ferfunc.o Fevol.o Fcdvcs.o Fmsbar.o \
Fdvcs.o Ff2.o Freadpar.o Ferror.o Fgauss.o Fnd.o Fbigcns.o Fintegra.o Fxspace.o\
Finitc.o Fevolc.o Fgetmbgpd.o Fhj.o Fsplice.o Fcff.o Fderiv.o Fbca.o
_FOBJS = $(patsubst %,$(FODIR)/%,$(FOBJS))


# rule for compiling Fortran .F  files
# (prefix "F" in "F%.o" is used to avoid possible future clashes with
#  C-wrappers with the same name)
$(FODIR)/F%.o: %.F  header.F
	$(FC) -c -o $@ $< $(FFLAGS) $(CPPFLAGS)

# C object files
OBJS = flushout.o
_OBJS = $(patsubst %,$(ODIR)/%,$(OBJS))

# rule for compiling C files
$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(INCS)

# All object files together
export ALLOBJS = $(_OBJS) $(_FOBJS)

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
$(EXFODIR)/%.o: $(EXDIR)/%.F header.F
	$(FC) -c -o $@ $< $(FFLAGS)

# test and fit targets

$(TESTTARGETS): %: $(ALLOBJS) $(FODIR)/F%.o
	$(FC) -o ../Tests/$@ $^ $(FFLAGS) $(ADACFLIB)

# main fit target  (a variant withouth PGPLOT)
fit: CPPFLAGS = $(USEPGPLOT)
fit: $(ALLOBJS) $(FODIR)/Ffcn.o $(FODIR)/Ffit.o
	$(FC) -o ../fits/$@ $^ $(FFLAGS) $(CERNLIBS) $(ADACFLIB) $(ALLPGPLOTLIBS)

# targets for processing/presenting GPDs resulting from fit
$(FITTARGETS): %: $(ALLOBJS) $(EXFODIR)/%.o
	$(FC) -o ../fits/$@ $^ $(FFLAGS)

# mma targets

gepard.exe: CPPFLAGS := -D MMA  $(USEPGPLOT)
gepard.exe: $(ALLOBJS) $(FODIR)/Ffcn.o $(_MMAOBJS)  $(FODIR)/Fmninterface.o
	$(FC) -o ../fits/$@/$(SYS)/$@ $^ $(FFLAGS) $(CERNLIBS) $(ADACFLIB) -L$(MLLIBDIR) -l$(MLLIB) $(ALLPGPLOTLIBS) $(MLEXTRA)

# targets for examples
#
$(EXTARGETS): %: $(ALLOBJS) $(EXFODIR)/%.o
	$(FC) -o ../ex/$@ $^ $(FFLAGS) $(ADACFLIB)


.PHONY: clean
clean:
	-rm -f $(FODIR)/F*.o $(ODIR)/*.o $(EXFODIR)/*.o $(MMAODIR)/*
	-rm -f $(patsubst %,../Tests/%,$(TESTTARGETS))
	-rm -f $(patsubst %,../Tests/%.exe,$(TESTTARGETS))
	-rm -f $(patsubst %,../fits/%,$(FITTARGETS))
	-rm -f $(patsubst %,../fits/%.exe,$(FITTARGETS))
	-rm -f ../fits/fit ../fits/fit.exe 
	-rm -f $(patsubst %,../ex/%,$(EXTARGETS))
	-rm -f $(patsubst %,../ex/%.exe,$(EXTARGETS))
	-rm -f ../fits/gepard.exe/$(SYS)/gepard.exe
