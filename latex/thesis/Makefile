all: output/thesis.pdf output/table_workbench.pdf

output/thesis.pdf: thesis.tex thesis.bib infufrgs/inputs/*
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis.tex
	bibtex output/thesis
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis.tex
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis.tex

output/table_workbench.pdf: table_workbench.tex infufrgs/inputs/*
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output table_workbench.tex

.PHONY: clean
clean:
	rm -rf ./output/*.{aux,bbl,blg,lof,log,lot,out,toc}

