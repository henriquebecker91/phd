#!/bin/sh

pdflatex revised_furini.tex
bibtex revised_furini
pdflatex revised_furini.tex
pdflatex revised_furini.tex

