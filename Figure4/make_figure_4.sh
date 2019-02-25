#!/bin/bash

echo 'This code reproduce Figure 4 in the main text'
echo 'Data are take from: Species fluctuations sustained by a cyclic succession at the edge of chaos'
echo 'Beninca et al. PNAS (2015)'
echo 'This code will take a very long time ...'
echo 'Steps:'
echo '1) Compute the Jacobian matrix with the regularized S-map from the multivariate time series'
echo '   using 40 different ratios of L1 to L2 norms (this is the slow step)'
echo '2) Construct an average Jacobian matrix from the ensemble using only those coefficients'
echo '   that minimize training and test error'
echo '3) Compute the VCR (Tr[J]) and its confidence intervals'

cd Code/
Rscript main.R ### Create the ensemble for the computation of the VCR
cd outputFiles/
Rscript weighted_ensemble_method.r ### Make the weighted jacobian matrix
Rscript final_computation_vcr.r ### Compute the VCR with confidence intervals
cd ../../
## Now make the Figure
cd Figure/
pdflatex Fig4.tex
open Fig4.pdf


