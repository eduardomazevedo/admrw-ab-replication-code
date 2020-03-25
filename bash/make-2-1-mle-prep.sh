#!/bin/bash
#$ -N make-2-1-mle-prep
#$ -hold_jid make-1-r
#$ -q short.q
#$ -o ./log-bash/
#$ -e ./log-bash/
#$ -cwd

############################
## MLE in MATLAB.
run="matlab -nodisplay < "
echo "symbolic differentiation of the t distribution"
$run ./matlab/scripts/create_tpdf_gradient.m
