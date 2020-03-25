#!/bin/bash
#$ -N f-tabulation-collect
#$ -hold_jid f-tabulation-array-job
#$ -q short.q
#$ -o ./log-bash/
#$ -e ./log-bash/
#$ -cwd

############################
## MLE in MATLAB.
run="matlab -nodisplay < "

$run ./matlab/exploratory/f_tabulation_collect.m
