#!/bin/bash
#$ -N make-1-r
#$ -q short.q
#$ -o ./log-bash/
#$ -e ./log-bash/
#$ -cwd

############################
## Data cleaning steps in R.

#1
Rscript r/clean-upstream-foray-dataset.R

#2
Rscript r/extend-probability-legit.R

#3
Rscript r/summary-stats.R
rm Rplots.pdf

#4
Rscript r/create-deconvolution-data.R
