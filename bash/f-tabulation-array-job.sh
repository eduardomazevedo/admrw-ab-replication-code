#!/bin/bash
#$ -N f_tabulation_array
#$ -q short.q
#$ -t 1-400
#$ -o ./log-bash/
#$ -e ./log-bash/
#$ -cwd

############################
# Tabulate f
matlab -nodisplay < ./matlab/exploratory/f_tabulation_array_job.m
