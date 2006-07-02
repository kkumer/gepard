# $Id$
# Master Makefile for gepard project
# 2006-04-15 

# Interesting targets:
# 	radcorr  -  program for producing Fig. 1 in letter
# 	scaledep -  program for producing Fig. 2 in letter
# 	gepard.exe - Mathematica-installable package
# 	install - putting gepard.exe in Mathematica path


# Where and how gepard.exe should be installed
OWNER = $(USER) 
GROUP = $(USER)
INSTALL = install
BINDIR = $(HOME)/.mma/Applications/gepard.exe/Linux/


# targets
SRCTARGETS = radcorr scaledep gepard.exe fit plotsigma
DOCTARGETS = pdf html
.PHONY: $(SRCTARGETS)

all: $(SRCTARGETS) $(DOCTARGETS)

$(SRCTARGETS):
	$(MAKE) -C src $@
	mv src/$@ .

install: gepard.exe
	$(INSTALL) -o $(OWNER) -g $(GROUP) -m 755 $< $(BINDIR)


doc: html pdf

pdf: tex
	$(MAKE) -C doc/tex pdf

tex:
	ln -s doc/robodoc.rc .
	robodoc  --latex --singledoc --doc ./doc/tex/gepard 
	rm robodoc.rc
	
html:
	ln -s doc/robodoc.rc .
	robodoc  --html --multidoc --index --doc ./doc/html
	rm robodoc.rc
	

.PHONY: clean
clean:
	$(MAKE) -C src clean
	$(MAKE) -C doc/tex clean
	-rm doc/html/*
	-rm -f $(SRCTARGETS)
