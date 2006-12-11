# Master Makefile for gepard project

# Interesting targets:
# 	radcorr  -  program for producing Fig. 1 in letter
# 	scaledep -  program for producing Fig. 2 in letter
# 	auxns    -  comparing MSBAR and CSBAR schemes
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
export SRCTARGETS = radcorr scaledep fit test auxtest fit_nopgplot houches accuracy
export EXTARGETS = aux auxns anatomyNS anatomy evoldep\
                   radNLONS radNLO radQ radQNS radNNLO scales scalesNS scalesNNLO

.PHONY: $(SRCTARGETS) $(EXTARGETS)
DOCTARGETS = pdf html

all: $(SRCTARGETS) $(EXTARGETS) $(DOCTARGETS)

$(SRCTARGETS):
	$(MAKE) -C src $@

$(EXTARGETS):
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
	-rm -f Tests/*dat
	-rm -f ex/*dat

.PHONY: clean
clean:
	$(MAKE) -C src clean
	$(MAKE) -C doc/tex clean
	-rm -rf doc/html/*
	-rm -f fits/*.{min,out,ps} fits/gmon.out
	-rm -f Tests/*dat Tests/gmon.out
	-rm -f ex/*dat ex/*eps
