#!/bin/bash
#$ -N make-3-3-triage
#$ -hold_jid make-2-3-mle-collect
#$ -q short.q
#$ -o ./log-bash/
#$ -e ./log-bash/
#$ -cwd

############################
## Triage analysis.

#1
Rscript r/triage-rmarkdown.R
