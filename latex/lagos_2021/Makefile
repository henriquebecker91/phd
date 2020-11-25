output/thesis_proposal.pdf: clean lagos.tex lagos.bib entcs_format/*
	TEXINPUTS=".:./entcs_format:" pdflatex --output-directory output lagos.tex
	BSTINPUTS=".:./entcs_format:" bibtex output/lagos
	TEXINPUTS=".:./entcs_format:" pdflatex --output-directory output lagos.tex
	TEXINPUTS=".:./entcs_format:" pdflatex --output-directory output lagos.tex

.PHONY: clean
clean:
	rm -rf ./output/*.{aux,bbl,blg,lof,log,lot,out,toc}

