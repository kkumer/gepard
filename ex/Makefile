# Makefile for gepard examples

# Interesting targets:
# 	*0.dat    -  produce data files *0.dat, *1.dat, ...
# 	*.eps     -  produce figs


# targets
export DATTARGETS = aux0.dat auxns0.dat scalesNS0.dat scales0.dat scalesNNLO0.dat
export FIGTARGETS = radNLONS.eps radNLO.eps evolutNS.eps \
                   evolut.eps radQ.eps radNNLONS.eps radNNLO.eps 

# Long-running targets
export LONG = anatomyNS.eps anatomy.eps

all: $(DATTARGETS) $(FIGTARGETS)

everything: all $(LONG)

# Don't delete intermediate files
.PRECIOUS: %0.dat

# rule for doing data files
%0.dat: %
	$^

# rule for doing figs
%.eps: %0.dat
	./plt -b $*


.PHONY: datclean
datclean:
	-rm -f *.dat

.PHONY: clean
clean:
	-rm -f *.dat *.eps

