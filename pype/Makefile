# Master Makefile for pype


# ------------------------------------------------------------------------  
# ---- BEGIN of system dependent stuff (fix it by hand where needed!) ----
# ------------------------------------------------------------------------  


# -- 0. Compilation options

FC = gfortran
CC = gcc
OPT_CFLAGS = -O3 -fPIC
PYSUFFIX = $(shell python3-config --extension-suffix)

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
CFLAGS = $(OPT_CFLAGS)


ifdef WINDIR
	#PYINCLUDES = -I c:\Python25\include
	PYINCLUDES = -I /cygdrive/c/Python25/include

	#PYLIBS = -L c:\Python25\libs -lpython25 -lm 
	PYLIBS = -L /cygdrive/c/Python25/libs -lpython25 -lm
else
	PYINCLUDES = -I/usr/include/python2.6 -I/usr/include/python2.6
	PYLIBS = -lpthread -ldl -lutil -lm -lpython2.6
endif

export CC CFLAGS PYINCLUDES PYLIBS

# -- 2. MathLink related things
#
# Put your version of Mathematica here and it's root dir (final slash needed!)
export MMAVERSION=7.0
export MMAROOT = /usr/local/Wolfram/Mathematica/
#export MMAROOT = /cygdrive/c/Program\ Files/Wolfram\ Research/Mathematica/
# export MMAROOT = /psi/math-
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



# ------------------------------------------------------------------------  
# ---- END of system dependent stuff                                  ----
# ------------------------------------------------------------------------  


# targets
export TESTTARGETS = test_driver
export MMATARGETS = xs_mma.exe

.PHONY: $(TESTTARGETS)  $(MMATARGETS) 

all: $(TESTTARGETS) $(MMATARGETS)

tests: $(TESTTARGETS)

mma: $(MMATARGETS)


$(TESTTARGETS) $(MMATARGETS):
	$(MAKE) -C ifaces $@


optModel.so: optModel.f
	f2py --fcompiler=$(FC) -c -m optModel $^ 
	mv optModel$(PYSUFFIX) optModel.so


.PHONY: clean 
clean:
	$(MAKE) -C ifaces clean
	-rm -f optModel.so
