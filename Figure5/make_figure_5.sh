#!/bin/bash

echo 'Causality test of VCR and environmental variables'
declare -a sample2=('divergence,temperature' 'divergence,Waves')

for j in "${sample2[@]}"
do
	echo '... Running'
	Rscript main.R $j > /dev/null
	echo '... done'
done

echo 'Now make the figure'
cd outputTables/Figure/
pdflatex Figure5.tex > /dev/null
rm *.log *.aux
open Figure5.pdf

