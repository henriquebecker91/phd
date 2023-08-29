output/thesis_slides.pdf: thesis_slides.tex thesis_slides.bib
	pdflatex thesis_slides.tex
	bibtex thesis_slides
	pdflatex thesis_slides.tex
	pdflatex thesis_slides.tex

.PHONY: clean
clean:
	rm -rf ./*.{aux,bbl,blg,lof,log,lot,out,toc,nav,vrb,snm}

