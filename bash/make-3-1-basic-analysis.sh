#!/bin/bash
#$ -N make-3-1-basic-analysis
#$ -hold_jid make-2-3-mle-collect
#$ -q short.q
#$ -t 1-9
#$ -o ./log-bash/
#$ -e ./log-bash/
#$ -cwd

# Run
matlab -nodesktop -r 'addpath ./matlab/scripts; basic_analysis_publish;'
