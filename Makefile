# Master Makefile for gepard project 2008-07-10 

# Interesting targets:
# 	radcorr  -  program for producing Fig. 1 in letter
# 	scaledep -  program for producing Fig. 2 in letter
# 	auxns    -  comparing MSBAR and CSBAR schemes
# 	test     -  tests both DVCS and DIS routines
# 	auxtest  -  tests DVCS \mathcal{H} calculation
# 	fit      -  fitting GPD ansatz to DVCS and DIS data
# 	fit_nopgplot -  same as fit, but without plotting
# 	gepard.exe - Mathematica interface
# 	html     -  HTML documentation
# 	pdf      -  LaTeX -> PDF documentation
#
# 	houches  -  comparison to Les Houches benchmark
# 	accuracy -  Analysis of accuracy and SPEED
#
# For compiling 'fit' without PGPLOT (no plotting) do
# 	  make NOPGPLOT=1 fit
# 	(This is what fit_nopgplot target does, but with
# 	 different name for executable.)
#
# For debugging/profiling call like this:
#     make FFLAGS='-g -pg' 
#
# For compiling without cygwin.dll dependency
# 	  make FFLAGS='-O -mno-cygwin' CFLAGS='-O2 -mno-cygwin' test
#

# ------------------------------------------------------------------------  
# ---- BEGIN of system dependent stuff (fix it by hand where needed!) ----
# ------------------------------------------------------------------------  

# -- 1. MINUIT related things
#
# Location and links to CERNLIB's kernlib and packlib
export CERNLIBS =  -L$(HOME)/local/lib -lpacklib -lkernlib

# -- 2. MathLink related things
#
# Version of Mathematica
export MMAVERSION=5.2
ifdef WINDIR
  export SYS = Windows
  export MLDIR=/cygdrive/c/Program\ Files/Wolfram\ Research/Mathematica/$(MMAVERSION)/AddOns/MathLink/DeveloperKit/$(SYS)/CompilerAdditions/mldev32
  export MPREP = $(MLDIR)/bin/mprep
  export MLINCDIR = $(MLDIR)/include
  export MLLIBDIR = $(MLDIR)/lib 
  export MLLIB = ml32i2w
  export MLEXTRA = -mwindows -DWIN32_MATHLINK
else
  ifdef NOT64
    export SYS = Linux
  else 
    export SYS = Linux-x86-64
  endif
  ifeq '$(MMAVERSION)' '6.0'
    export MLDIR = /usr/local/Wolfram/Mathematica/$(MMAVERSION)/SystemFiles/Links/MathLink/DeveloperKit/$(SYS)/CompilerAdditions
  else
    export MLDIR=/usr/local/Wolfram/Mathematica/$(MMAVERSION)/AddOns/MathLink/DeveloperKit/$(SYS)/CompilerAdditions
  endif
  export MPREP = $(MLDIR)/mprep
  export MLINCDIR = $(MLDIR)
  export MLLIBDIR = $(MLDIR)
  export MLLIB = ML
  export MLEXTRA = -lpthread
endif


# -- 3. PGPLOT related things
#  
# Location and links to pgplot libs (if you have them. 
# If not, compile fit_nopgplot instead of fit'.)
export PGPLOTLIBS = -L$(HOME)/local/lib/pgplot -lpgplot
# If you have pgplot libs, but without /XSERVE driver, comment 
# out the next three lines (should be automatic on Windows)
ifndef WINDIR
  export X11LIBS = -L/usr/X11R6/lib -lX11 -lpng
endif
# All PGPLOT libs together for easier reference:
export ALLPGPLOTLIBS = $(PGPLOTLIBS) $(X11LIBS)

# ------------------------------------------------------------------------  
# ---- END of system dependent stuff                                  ----
# ------------------------------------------------------------------------  


# targets
export TESTTARGETS = radcorr scaledep test auxtest houches accuracy
export EXTARGETS = auxsi auxns anatomyNS anatomy radNLONS radNLO evolutNS evolut radQ \
                   radNNLONS radNNLO scalesNS scales scalesNNLO
export FITTARGETS = fitres fitpdfs slope contours
export MMATARGETS = gepard.exe

.PHONY: $(TESTTARGETS) $(EXTARGETS) $(FITTARGETS) $(MMATARGETS)
DOCTARGETS = pdf html htmlnocss

all: $(TESTTARGETS) $(EXTARGETS) fit fit_nopgplot  $(FITTARGETS) $(MMATARGETS) $(DOCTARGETS)

tests: $(TESTTARGETS)

$(TESTTARGETS) $(EXTARGETS) fit fit_nopgplot $(FITTARGETS) $(MMATARGETS):
	$(MAKE) -C src $@

doc: html pdf

pdf: tex
	$(MAKE) -C doc/tex allpdfs

tex:
	robodoc --rc doc/robodoc.rc --latex --singledoc --toc --index --doc ./doc/tex/gepard-api
	
html:
	robodoc --rc doc/robodoc.rc  --html --multidoc --index --doc ./doc/html --css ./doc/gepard.css

htmlnocss:
	robodoc --rc doc/robodoc.rc  --html --multidoc --index --doc ./doc/html

.PHONY: rmfig
rmfig:
	-rm -f Tests/*dat
	-rm -f ex/*dat

.PHONY: clean
clean:
	$(MAKE) -C src clean
	$(MAKE) -C doc/tex clean
	-rm -rf doc/html/*
	-rm -f fits/*.{mnt,out,ps,eps,out}
	-rm -f Tests/*{dat,out,eps}
	-rm -f ex/*{dat,eps}
