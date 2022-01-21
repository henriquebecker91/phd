output/thesis.pdf: clean thesis.tex thesis.bib infufrgs/inputs/*
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis.tex
	bibtex output/thesis
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis.tex
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis.tex

.PHONY: clean
clean:
	rm -rf ./output/*.{aux,bbl,blg,lof,log,lot,out,toc}

