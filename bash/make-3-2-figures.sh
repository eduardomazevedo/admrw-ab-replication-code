#!/bin/bash
#$ -N make-3-2-figures
#$ -hold_jid make-2-3-mle-collect
#$ -q short.q
#$ -o ./log-bash/
#$ -e ./log-bash/
#$ -cwd

############################
run="matlab -nodisplay < "
echo "Produce MATLAB figures."
$run ./matlab/scripts/make_figures.m
$run ./matlab/scripts/relevant_range.m
$run ./matlab/scripts/placebo_metrics.m
$run ./matlab/scripts/make_table_mleresults.m
$run ./matlab/scripts/make_table_mleresults_disaggregated.m
$run ./matlab/scripts/make_table_mle_disaggregated2.m
