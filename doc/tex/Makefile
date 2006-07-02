pdf: dvi
	dvipdf gepard.dvi

ps: gepard.dvi
	dvips -o gepard.ps gepard.dvi

dvi: gepard.tex
	latex gepard.tex
	latex gepard.tex
	makeindex gepard.idx
	latex gepard.tex

.PHONY: clean
clean:
	-rm -f gepard.*
