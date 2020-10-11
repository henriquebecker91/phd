#!/bin/sh

pdflatex --shell-escape revised_PPG2KP.tex
bibtex revised_PPG2KP
pdflatex --shell-escape revised_PPG2KP.tex
pdflatex --shell-escape revised_PPG2KP.tex

