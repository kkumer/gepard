# Master Makefile for gepard project 2010-06-25 

# Interesting targets:
# 	test     -  tests both DVCS and DIS routines
# 	fit      -  fitting GPD ansatz to DVCS and DIS data
# 	gepard.exe - Mathematica interface
#   pygepard.so - python interface
# 	html     -  HTML documentation
# 	pdf      -  LaTeX -> PDF documentation
# 	radcorr  -  program for producing Fig. 1 in hep-ph/0605237
# 	scaledep -  program for producing Fig. 2 in hep-ph/0605237
# 	auxns    -  comparing MSBAR and CSBAR schemes
# 	auxtest  -  tests DVCS \mathcal{H} calculation
#
# 	houches  -  comparison to Les Houches DIS benchmark (broken)
# 	accuracy -  Analysis of accuracy and SPEED of DVCS calculation
#
# For compiling 'fit' and 'gepard.exe' without PGPLOT (no plotting) do
# 	  make NOPGPLOT=1 fit
#
# For debugging/profiling call like this:
#     make DEBUG=1 <target> 
#     make PROFILE=1 <target> 
#
# For compiling on Windows to create binary without cygwin.dll dependency
# 	  make NOPGPLOT=1 FFLAGS='-O -mno-cygwin' CFLAGS='-O2 -mno-cygwin' gepard.exe
#

# ------------------------------------------------------------------------  
# ---- BEGIN of system dependent stuff (fix it by hand where needed!) ----
# ------------------------------------------------------------------------  


# -- 0. Compilation options
#

# Fortran compiler dependent options (uncomment only one)
#
## [--1--] GNU g77  (also cygwin)
#FC = g77
#CMP_FFLAGS = -Wall
#OPT_FFLAGS = -O3
#OPT_CFLAGS = -O3 
#
## [--2--] GNU gfortran
#FC = gfortran
#CMP_FFLAGS = -Wall
#OPT_FFLAGS = -O3 -ffast-math -funroll-all-loops  -ftree-vectorize -fPIC
#OPT_CFLAGS = -O3 -fPIC
#
## [--2b--] GNU gfortran + OpenMP parallelization
FC = gfortran
CC = gcc
CMP_FFLAGS = -Wall
OPT_FFLAGS = -O3 -fopenmp  -fPIC
OPT_CFLAGS = -O3

## [--3--] INTEL + OpenMP parallelization
#FC = ifort
#CC = icc
#CMP_FFLAGS = -warn all -nofor_main 
#OPT_FFLAGS = -O3 -ipo -openmp
#OPT_CFLAGS = -O3 -ipo

## [--4--] NAG f95  
#FC = f95
#CMP_FFLAGS = -colour -w=x77 -w=obs
#OPT_FFLAGS = -O3
#OPT_CFLAGS = -O3 
#CPPFLAGS = -DNAGf90Fortran

# optimized, debug and profiling modes
#
ifdef DEBUG
  OPT_FFLAGS = -g
  OPT_CFLAGS = -g
endif
ifdef PROFILE
  OPT_FFLAGS = -g -pg
  OPT_CFLAGS = -g -pg
endif
# adding everything + compiler-specific options
FFLAGS = -I. $(OPT_FFLAGS) $(CMP_FFLAGS)
CFLAGS = $(OPT_CFLAGS)

export FC CC FFLAGS CFLAGS CPPFLAGS


# -- 1. MINUIT related things
#
# Location and links to CERNLIB's kernlib and packlib
# FIXME: there is compiler dependence in these libs!
export CERNLIBS =  -L$(HOME)/local/lib -lpacklib_$(FC) -lkernlib_$(FC) 
ifdef DEBUG
  export CERNLIBS =  -L$(HOME)/local/lib -lpacklib_$(FC)_dbg -lkernlib_$(FC)_dbg
endif
ifdef PROFILE
  export CERNLIBS =  -L$(HOME)/local/lib -lpacklib_$(FC)_prof -lkernlib_$(FC)_prof
endif

# -- 2. MathLink related things
#
# Put your version of Mathematica here and it's root dir (final slash needed!)
export MMAVERSION=7.0
export MMAROOT = /usr/local/Wolfram/Mathematica/
# export MMAROOT = /cygdrive/c/Program\ Files/Wolfram\ Research/Mathematica/
#export MMAROOT = /psi/math-
ifdef WINDIR
  export SYS = Windows
  ifeq '$(MMAVERSION)' '5.0'
	export MLDIR=$(MMAROOT)$(MMAVERSION)/AddOns/MathLink/DeveloperKit/$(SYS)/CompilerAdditions/mldev32
  else
	export MLDIR=$(MMAROOT)$(MMAVERSION)/SystemFiles/Links/MathLink/DeveloperKit/$(SYS)/CompilerAdditions/cygwin
  endif
  export MPREP = $(MLDIR)/bin/mprep
  export MLINCDIR = $(MLDIR)/include
  export MLLIBDIR = $(MLDIR)/lib 
  ifeq '$(MMAVERSION)' '5.0'
    export MLLIB = ml32i2w
  else
    export MLLIB = ML32i3
  endif
  export MLEXTRA = -mwindows -DWIN32_MATHLINK
else
  ifdef NOT64
    export SYS = Linux
  else 
    export SYS = Linux-x86-64
  endif
  ifeq '$(MMAVERSION)' '5.0'
    export MLDIR=$(MMAROOT)$(MMAVERSION)/AddOns/MathLink/DeveloperKit/$(SYS)/CompilerAdditions
  else
    export MLDIR=$(MMAROOT)$(MMAVERSION)/SystemFiles/Links/MathLink/DeveloperKit/$(SYS)/CompilerAdditions
  endif
  export MPREP = $(MLDIR)/mprep
  export MLINCDIR = $(MLDIR)
  export MLLIBDIR = $(MLDIR)
  ifeq '$(MMAVERSION)' '5.0'
	export MLLIB = ML
	export MLEXTRA = -lpthread
  else
	ifdef NOT64
	  export MLLIB = ML32i3
	else
	  export MLLIB = ML64i3
	endif
	export MLEXTRA = -lpthread -lrt
  endif
endif


# -- 3. PGPLOT related things
#  
# Location and links to pgplot libs (if you have them. 
# If not, compile fit and gepard.exe with NOPGPLOT=1'.)
export PGPLOTLIBS = -L$(HOME)/local/lib/pgplot -lpgplot
# If you have pgplot libs, but without /XSERVE driver, comment 
# out the next three lines (should be automatic on Windows)
ifndef WINDIR
  export X11LIBS = -L/usr/X11R6/lib -lX11 -lpng
endif
# All PGPLOT libs together for easier reference:
ifndef NOPGPLOT
  export ALLPGPLOTLIBS = $(PGPLOTLIBS) $(X11LIBS)
endif


# ------------------------------------------------------------------------  
# ---- END of system dependent stuff                                  ----
# ------------------------------------------------------------------------  


# targets
export TESTTARGETS = radcorr scaledep test auxtest houches accuracy
export EXTARGETS = auxsi auxns anatomyNS anatomy radNLONS radNLO evolutNS evolut radQ \
                   radNNLONS radNNLO scalesNS scales scalesNNLO
export FITTARGETS = fitres fitres2 fitpdfs slope contours
export MMATARGETS = gepard.exe int2f1.exe dvem.exe
export PYTARGETS = pygepard.so

.PHONY: $(TESTTARGETS) $(EXTARGETS) $(FITTARGETS) $(MMATARGETS) $(PYTARGETS)
DOCTARGETS = pdf html htmlnocss

all: $(TESTTARGETS) $(EXTARGETS) fit $(FITTARGETS) $(MMATARGETS) $(PYTARGETS) $(DOCTARGETS)

tests: $(TESTTARGETS)

$(TESTTARGETS) $(EXTARGETS) fit $(FITTARGETS) $(MMATARGETS) $(PYTARGETS):
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
