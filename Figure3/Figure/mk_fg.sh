#!/bin/bash

Rscript median_error.r
pdflatex Figure3.tex > /dev/null 2>&1
open Figure3.pdf
