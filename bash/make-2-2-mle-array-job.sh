#!/bin/bash
#$ -N make-2-2-mle-array-job
#$ -hold_jid make-2-1-mle-prep
#$ -q all.q
#$ -t 1-54
#$ -o ./log-bash/
#$ -e ./log-bash/
#$ -cwd

############################
## MLE in MATLAB.
# Prep

# Run MLE
# matlab -nodesktop -r 'addpath ./matlab/scripts; mle_array_job;'
matlab -nodisplay < ./matlab/scripts/mle_array_job.m
