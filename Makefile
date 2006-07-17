# Master Makefile for gepard project

# Interesting targets:
# 	radcorr  -  program for producing Fig. 1 in letter
# 	scaledep -  program for producing Fig. 2 in letter
# 	test     -  tests both DVCS and DIS routines
# 	fit      -  fitting GPD ansatz to DVCS and DIS data
# 	html     -  HTML documentation
# 	pdf      -  LaTeX -> PDF documentation
#
# For debugging/profiling call like this:
#  make DEBUG='-g -pg' and change -ladacf -> -ladacf_prof
# For optimization call like this:
#  make DEBUG='-O2'



# targets
SRCTARGETS = radcorr scaledep fit test auxtest
DOCTARGETS = pdf html
.PHONY: $(SRCTARGETS)

all: $(SRCTARGETS) $(DOCTARGETS)

$(SRCTARGETS):
	$(MAKE) -C src $@
	mv src/$@ .


doc: html pdf

pdf: tex
	$(MAKE) -C doc/tex pdf

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
	-rm fits/FIT.OUT fits/gmon.out fits/*.ps
	-rm FIG*DAT gmon.out

.PHONY: distclean
distclean:
	$(MAKE) -C src clean
	$(MAKE) -C doc/tex clean
	-rm fits/FIT.OUT fits/gmon.out fits/*.ps
	-rm FIG*DAT gmon.out
	-rm -rf .svn Tests/.svn doc/.svn doc/html/.svn doc/tex/.svn \
		   	fits/.svn mma/.svn mma/.svn mma/gepard.exe/.svn \
		mma/gepard.exe/Linux/.svn mma/gepard.exe/Windows/.svn \
		src/.svn src/Fobj/.svn src/adacf/.svn  src/adacf/obj/.svn \
		src/obj/.svn
