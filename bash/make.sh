#!/bin/bash

# Clean data and summary stats on r.
qsub ./bash/make-1-r.sh

# Maximum likelihood estimation.
qsub ./bash/make-2-1-mle-prep.sh
qsub ./bash/make-2-2-mle-array-job.sh
qsub ./bash/make-2-2-5-mle-array-job-disaggregated.sh
qsub ./bash/make-2-3-mle-collect.sh

# Analysis
qsub ./bash/make-3-1-basic-analysis.sh
qsub ./bash/make-3-2-figures.sh
qsub ./bash/make-3-3-triage.sh
