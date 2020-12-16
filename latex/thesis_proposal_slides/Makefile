output/thesis_proposal_slides.pdf: clean thesis_proposal_slides.tex thesis_proposal_slides.bib
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis_proposal_slides.tex
	bibtex output/thesis_proposal_slides
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis_proposal_slides.tex
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output thesis_proposal_slides.tex

.PHONY: clean
clean:
	rm -rf ./output/*.{aux,bbl,blg,lof,log,lot,out,toc,nav,pdf,vrb,snm}

