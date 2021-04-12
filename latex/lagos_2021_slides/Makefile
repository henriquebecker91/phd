output/lagos_2021_slides.pdf: clean lagos_2021_slides.tex lagos_2021_slides.bib
	TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output lagos_2021_slides.tex
	#bibtex output/lagos_2021_slides
	#TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output lagos_2021_slides.tex
	#TEXINPUTS=".:./infufrgs/inputs:" pdflatex --output-directory output lagos_2021_slides.tex

.PHONY: clean
clean:
	rm -rf ./output/*.{aux,bbl,blg,lof,log,lot,out,toc,nav,pdf,vrb,snm}

