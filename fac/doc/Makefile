doc: manual.tex
	latex manual
	bibtex manual
	makeindex manual
	latex manual
	latex manual
	dvips -o manual.ps manual

clean:
	rm -rf *~ 

