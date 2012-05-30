allpdfs: gepard.pdf gepard-api.pdf

%.pdf: %.dvi
	dvipdf $^

%.ps: %.dvi
	dvips -o $@ $^

# FIXME: write/find more sophisticated rule
%.dvi: %.tex
	latex $^
	latex $^
	-makeindex $*.idx
	latex $^

.PHONY: clean
clean:
	-rm -f gepard-api.* gepard.log gepard.aux gepard.dvi \
	   	gepard.pdf gepard.ps gepard.toc gepard.blg gepard.out

realclean:
	-rm -f gepard-api.* gepard.log gepard.aux gepard.dvi \
	   	gepard.pdf gepard.ps gepard.blg gepard.out gepard.bbl
