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
	-rm -f gepard-api.* gepard.{log,aux,dvi,pdf,ps,toc,blg}

realclean:
	-rm -f gepard-api.* gepard.{log,aux,dvi,pdf,ps,blg,bbl}
