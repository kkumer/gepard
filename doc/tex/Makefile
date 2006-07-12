dvi: gepard.tex
	latex gepard.tex

pdf: apidvi
	dvipdf gepard-api.dvi

apips: gepard-api.dvi
	dvips -o gepard-api.ps gepard-api.dvi

apidvi: gepard-api.tex
	latex gepard-api.tex
	latex gepard-api.tex
	makeindex gepard-api.idx
	latex gepard-api.tex

.PHONY: clean
clean:
	-rm -f gepard-api.* gepard.{log,aux,dvi,pdf,ps,blg,bbl}
