ALL: latex.out/manuscript.pdf latex.out/cover_letter.pdf

latex.out:
	mkdir latex.out

SV_FILES= sv/spmpsci.bst sv/svglov3.clo sv/svjour3.cls

#EXTERNAL_FIGURES= T1a.pdf

latex.out/manuscript.pdf: Makefile $(SV_FILES) latex.out revised_PPG2KP.tex mybib.bib
	TEXINPUTS=".:./itor:" pdflatex --shell-escape --output-directory latex.out revised_PPG2KP.tex
	BSTINPUTS=".:./itor:" bibtex latex.out/revised_PPG2KP
	TEXINPUTS=".:./itor:" pdflatex --shell-escape --output-directory latex.out revised_PPG2KP.tex
	TEXINPUTS=".:./itor:" pdflatex --shell-escape --output-directory latex.out revised_PPG2KP.tex

latex.out/cover_letter.pdf: Makefile $(SV_FILES) cover_letter.tex
	TEXINPUTS=".:./sv:" pdflatex --output-directory latex.out cover_letter.tex
	BSTINPUTS=".:./sv:" bibtex latex.out/cover_letter
	TEXINPUTS=".:./sv:" pdflatex --output-directory latex.out cover_letter.tex
	TEXINPUTS=".:./sv:" pdflatex --output-directory latex.out cover_letter.tex

.PHONY: clean
clean:
	rm -f ./latex.out/*.{aux,bbl,blg,lof,log,lot,out,toc,spl}

