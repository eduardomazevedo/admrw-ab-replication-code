#!/bin/bash
#$ -N make-2-3-mle-collect
#$ -hold_jid make-2-2-mle-array-job,make-2-2-5-mle-array-job-disaggregated
#$ -q short.q
#$ -o ./log-bash/
#$ -e ./log-bash/
#$ -cwd

############################
## MLE in MATLAB.
run="matlab -nodisplay < "
echo "symbolic differentiation of the t distribution"
$run ./matlab/scripts/mle_array_collect.m
