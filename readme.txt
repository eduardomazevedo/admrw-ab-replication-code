---------
# Replication code for "A/B Testing with Fat Tails"
---------
This repository replicates the results of "A/B Testing with Fat Tails" by Azevedo, Deng, Montiel-Olea, Rao, and Weyl. The data used is proprietary. To obtain permission to use the data contact Microsoft Corporation at 1 Microsoft Way, Redmond, WA 98052.


---------
# General Comments about this repository:
---------
The initial data queries from Cosmos and data cleaning were done "upstream"
at MSR / Microsoft laptops. The initial analysis produces non-sensitive data
used in this "downstream" part of the analysis.

All scripts in this project (.m, .R, .sh ) should be run from the current
directory (i.e., the `.` directory).

Last Revision: March 2020.


---------
# Input from upstream analysis at MSR used in this repository:
---------
There is an initial phase of data cleaning that is done in the Microsoft Cosmos and in a Microsoft computer. This is output to `data/upstream-inputs`

There are two key files that will be required to reproduce our analysis.

The first file is our main dataset of A/B tests:

- `data/upstream-inputs/foray-dataset-non-sensitive.Rdata`

The second file is a small dataset of innovations from the relevance team
for which we have pre-flight metrics

- `data/upstream-inputs/triage-dataset-non-sensitive.Rdata`

*Please note that the .R scripts will not run without these two files.

There are two other files that were generated from the "upstream" data cleaning.

The first file is a summary of statistics that will be used to provide tables
and figures in the paper:

- `data/upstream-inputs/clean-cosmos-output.csv`

The second file is the result of the second iteration of our data audit:

- `data/upstream-inputs/surya-final-audit.csv`


-----------
# Required Software
-----------
R with data.table, tidyverse, glmnet, openxlsx, plotmo, gridExtra, openxlsx, rmarkdown.
MATLAB.


-----------
# How to run the files in this repository? (Pipeline)
-----------
To reproduce the entire analysis, run `bash/make.sh`. This will run the following steps.

1) Data cleaning and summary stats in R.
  - `bash/make-1-r.sh`
  - `r/clean-upstream-foray-dataset.R`
  - Filters observations for NTreatments = 1.
  - Filters observations for no triggers / filters.
  - Saves as `data/intermediate/data-for-probability-legit-fitting.Rdata`.
  - `r/extend-probability-legit`
  - Creates intermediate data with legit probabilities.
  - Currently this is fit with a LASSO regression using Delta and absolute value of Delta.
  - Currently the best fit is just a constant mean probability, and this seems pretty robust.
  - Saves to `data/intermediate/foray-dataset-with-fitted-probability.Rdata`.
  - `r/summary-stats.R`
  - Produces summary stats.
  - `r/create-deconvolution-data.R`
  - creates file `data/intermediate/convolution-data.csv` for estimation in MATLAB.
  - creates file `data/intermediate/convolution-data-disaggregated.csv` for disaggregated estimation.

2) Maximum likelihood estimation.
  - `bash/make-2-1-mle-prep.sh` and `matlab/scripts/create_tpdf_gradient.m`
  - Prepares for estimation by writing a function with derivatives of the t pdf.
  - `bash/make-2-2-mle-array-job` and `matlab/scripts/mle_array_job.m`
  - Runs several specifications of mle.
  - `bash/make-2-2-5-mle-array-job-disaggregated` and `matlab/scripts/mle_array_job_disaggregated.m`
  - run disaggregated mle.
  - `bash/make-3-3-mle-collect.sh` and `matlab/scripts/mle_array_collect.m`
  - Collect mle results into `matlab/mat/mle.mat`.

3) Analysis.
  - `bash/make-3-1-basic-analysis.sh` and `matlab/scripts/basic_analysis.m`
  - Creates reports for internal use.
  - `bash/make-3-2-figures.sh` and `matlab/scripts/make_figures.m`
  - Creates figures for the paper.
  - Produces output from the triage analysis.

4) Convert tables to latex (This is not included in the makefile because of cluster software requirements).
  - Kinda buggy and only works on Ed's mac.
  - First open the tables in .xlsx in Excel to get the numbers to update.
  - `bash/make-4-tables.sh`. Needs gnumeric and ssconvert.

Most numerical methods are implemented in the MATLAB Tweedie class. See MATLAB readme.


-----------
# List of output
-----------
- Figures:
  - Production function t distribution illustrative (`matlab/scripts/make_figures.m`).
  - Production function normal illustrative (`matlab/scripts/make_figures.m`).
  - Histogram (`r/summary-stats.R`).
  - Log log plots (`r/summary-stats.R`).
  - MLE estimates mean and scale (`matlab/scripts/make_figures.m`).
  - MLE estimates tail (`matlab/scripts/make_figures.m`).
  - Posterior mean (`matlab/scripts/make_figures.m`).
  - Estimated production function (`matlab/scripts/make_figures.m`).
  - Counterfactual from lean experimentation (`matlab/scripts/make_figures.m`).
- Tables:
  - Summary stats (`r/summary-stats.R`)

# Other files
- `r/start.R`: Loads paths and libraries for this project.
