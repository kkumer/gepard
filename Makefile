# Master Makefile for gepard project

# Interesting targets:
# 	radcorr  -  program for producing Fig. 1 in letter
# 	scaledep -  program for producing Fig. 2 in letter
# 	test     -  tests both DVCS and DIS routines
# 	auxtest  -  tests DVCS \mathcal{H} calculation
# 	fit      -  fitting GPD ansatz to DVCS and DIS data
# 	fit_nopgplot -  same as fit, but without plotting
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


# Location and links to CERNLIB's kernlib and packlib
export CERNLIBS =  -L$(HOME)/local/lib -lpacklib -lkernlib

# Location and links to pgplot libs (if you have them. 
# If not, compile with 'make NOPGPLOT=1 fit'.)
export PGPLOTLIBS = -L$(HOME)/local/lib/pgplot -lpgplot
# If you have pgplot libs, but without /XSERVE driver, comment 
# out the next three lines (should be automatic on Windows,
# but it's an ugly hack)
ifndef WINDIR
export X11LIBS = -L/usr/X11R6/lib -lX11
endif

# targets
SRCTARGETS = radcorr scaledep fit test auxtest fit_nopgplot houches accuracy
.PHONY: $(SRCTARGETS)
DOCTARGETS = pdf html

all: $(SRCTARGETS) $(DOCTARGETS)

$(SRCTARGETS):
	$(MAKE) -C src $@

doc: html pdf

pdf: tex
	$(MAKE) -C doc/tex allpdfs

tex:
	robodoc --rc doc/robodoc.rc --latex --singledoc --toc --index --doc ./doc/tex/gepard-api
	
html:
	robodoc --rc doc/robodoc.rc  --html --multidoc --index --doc ./doc/html

.PHONY: rmfig
rmfig:
	-rm FIG*DAT

.PHONY: clean
clean:
	$(MAKE) -C src clean
	$(MAKE) -C doc/tex clean
	-rm doc/html/*
	-rm -f $(SRCTARGETS)
	-rm -f $(patsubst %,src/%,$(SRCTARGETS))
	-rm -f $(patsubst %,%.exe,$(SRCTARGETS))
	-rm -f $(patsubst %,src/%.exe,$(SRCTARGETS))
	-rm fits/*.{min,out,ps} fits/gmon.out
	-rm Tests/FIG*DAT Tests/gmon.out
