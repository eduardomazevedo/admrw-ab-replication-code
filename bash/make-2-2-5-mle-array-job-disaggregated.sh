#!/bin/bash
#$ -N make-2-2-5-mle-array-job-disaggregated
#$ -hold_jid make-2-1-mle-prep
#$ -q short.q
#$ -t 1-12
#$ -o ./log-bash/
#$ -e ./log-bash/
#$ -cwd

############################
## MLE in MATLAB.
# Prep

# Run MLE
# matlab -nodesktop -r 'addpath ./matlab/scripts; mle_array_job;'
matlab -nodisplay < ./matlab/scripts/mle_array_job_disaggregated.m
