
all: manual3.pdf 

manual3.pdf:manual3.tex
	latex  manual3
	latex  manual3
	dvipdf  manual3.dvi	
clean: 
	rm -f  *.dvi *.log *.aux manual3.pdf	
