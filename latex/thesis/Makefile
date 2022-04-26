all: output/thesis.pdf output/table_workbench.pdf

output/thesis.pdf: Makefile thesis.tex related_problems.tex future_works.tex thesis.bib infufrgs/inputs/* plots/* diagrams/*
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis.tex
	bibtex output/thesis
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis.tex
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis.tex

output/table_workbench.pdf: Makefile table_workbench.tex infufrgs/inputs/* plots/* diagrams/*
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output table_workbench.tex
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output table_workbench.tex

.PHONY: clean
clean:
	rm -rf ./output/*.{aux,bbl,blg,lof,log,lot,out,toc}

