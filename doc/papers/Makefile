all: structure.ps excitation.ps pirr.ps aidr.ps ionization.ps

.SUFFIXES: .tex .ps

.tex.ps:
	latex $*
	bibtex $*
	sed -f fixref.sed $*.bbl > $*_tmp.bbl
	mv $*_tmp.bbl $*.bbl
	latex $*
	latex $*
	dvips -o $@ $*

single: 
	sed '1,1s/documentclass.*{/documentclass{/' \
	structure.tex > structure.tmp
	mv structure.tmp structure.tex
	sed '1,1s/documentclass.*{/documentclass{/' \
	excitation.tex > excitation.tmp
	mv excitation.tmp excitation.tex
	sed '1,1s/documentclass.*{/documentclass{/' \
	pirr.tex > pirr.tmp
	mv pirr.tmp pirr.tex
	sed '1,1s/documentclass.*{/documentclass{/' \
	aidr.tex > aidr.tmp
	mv aidr.tmp aidr.tex
	sed '1,1s/documentclass.*{/documentclass{/' \
	ionization.tex > ionization.tmp
	mv ionization.tmp ionization.tex


double: 
	sed '1,1s/documentclass.*{/documentclass[reviewcopy]{/' \
	structure.tex > structure.tmp
	mv structure.tmp structure.tex
	sed '1,1s/documentclass.*{/documentclass[reviewcopy]{/' \
	excitation.tex > excitation.tmp
	mv excitation.tmp excitation.tex
	sed '1,1s/documentclass.*{/documentclass[reviewcopy]{/' \
	pirr.tex > pirr.tmp
	mv pirr.tmp pirr.tex
	sed '1,1s/documentclass.*{/documentclass[reviewcopy]{/' \
	aidr.tex > aidr.tmp
	mv aidr.tmp aidr.tex
	sed '1,1s/documentclass.*{/documentclass[reviewcopy]{/' \
	ionization.tex > ionization.tmp
	mv ionization.tmp ionization.tex