output/thesis_proposal.pdf: thesis_proposal.tex thesis_proposal.bib infufrgs/inputs/*
	rm -rf ./output
	mkdir output
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis_proposal.tex
	bibtex output/thesis_proposal
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis_proposal.tex
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis_proposal.tex

