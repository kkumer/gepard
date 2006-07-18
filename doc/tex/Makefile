pdf: dvi
	dvipdf gepard.dvi

dvi: gepard.tex
	latex gepard.tex

apipdf: apidvi
	dvipdf gepard-api.dvi

apips: apidvi
	dvips -o gepard-api.ps gepard-api.dvi

apidvi: gepard-api.tex
	latex gepard-api.tex
	latex gepard-api.tex
	makeindex gepard-api.idx
	latex gepard-api.tex

.PHONY: clean
clean:
	-rm -f gepard-api.* gepard.{log,aux,dvi,pdf,ps,blg,bbl}
