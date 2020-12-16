ALL: clean output/thesis_proposal.pdf #output/appendix.pdf

#output/appendix.pdf: appendix.tex
#	pdflatex --output-directory output appendix.tex
#	pdflatex --output-directory output appendix.tex

output/thesis_proposal.pdf: lagos.tex lagos.bib entcs_format/*
	TEXINPUTS=".:./entcs_format:" pdflatex --output-directory output lagos.tex
	BSTINPUTS=".:./entcs_format:" bibtex output/lagos
	TEXINPUTS=".:./entcs_format:" pdflatex --output-directory output lagos.tex
	TEXINPUTS=".:./entcs_format:" pdflatex --output-directory output lagos.tex

.PHONY: clean
clean:
	rm -rf ./output/*.{aux,bbl,blg,lof,log,lot,out,toc}

