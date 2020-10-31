output/thesis_proposal.pdf: clean thesis_proposal.tex thesis_proposal.bib infufrgs/inputs/*
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis_proposal.tex
	bibtex output/thesis_proposal
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis_proposal.tex
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis_proposal.tex

.PHONY: clean
clean:
	rm -rf ./output/*.{aux,bbl,blg,lof,log,lot,out,toc}

