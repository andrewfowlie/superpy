docdir = ${datadir}/doc/${PACKAGE}

doc_DATA = nmssmManual.pdf \
           rpcManual.pdf \
           rpvManual.pdf \
           neutManual.pdf

EXTRA_DIST = \
	$(doc_DATA) \
	rpcManual.tex \
	rpvManual.tex \
	neutManual.tex \
	nmssmManual.tex

doc: $(doc_DATA)

rpcManual.pdf: rpcManual.tex
	latex rpcManual.tex
	latex rpcManual.tex
	dvips -f rpcManual > rpcManual.ps
	ps2pdf rpcManual.ps

rpvManual.pdf: rpvManual.tex
	latex rpvManual.tex
	latex rpvManual.tex
	dvips -f rpvManual > rpvManual.ps
	ps2pdf rpvManual.ps

neutManual.pdf: neutManual.tex
	latex neutManual.tex
	latex neutManual.tex
	dvips -f neutManual > neutManual.ps
	ps2pdf neutManual.ps

nmssmManual.pdf: nmssmManual.tex
	latex nmssmManual.tex
	latex nmssmManual.tex
	dvips -f nmssmManual > nmssmManual.ps
	ps2pdf nmssmManual.ps

clean-local:
	-rm -f *.log *.spl *.aux *.dvi *.ps *.toc *.out

distclean-local: clean-local
	-rm -f $(doc_DATA)
