#!/bin/sh

pdflatex revised_PPG2KP.tex
bibtex revised_PPG2KP
pdflatex revised_PPG2KP.tex
pdflatex revised_PPG2KP.tex

