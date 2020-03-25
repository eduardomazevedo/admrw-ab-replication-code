# Overview
  - All scripts should be run from the `.` directory.
  - For the pipeline, see the readme file for the entire project in the root folder.


# scripts
  -`create_t_pdf_gradient.m`: Creates functions for the gradient and hessian of the t distribution.
  -`mle_array_job.m` and `mle_array_collect.m`: performs maximum likelihood estimates assuming delta has a t distribution. Saves the result to `./mat/mle.mat`
  -`basic_analysis.m`: Plots basic results: goodness of fit, Bayesian correction, decomposition of past gains, and the production function. Expects a variable `metric_i` which is the number 1-11 of metric to look at.
  -`basic_analysis_publish.m`: Publishes the `basic_analysis.m` script. It is usually called by `basic_analysis.sh`, which runs it as an array job.
  - `make_figures.sh`: Produces the figures for the paper.


# classes
  -`Twee.m` has the functions that do all of the heavy lifting of the computations. See the help file.
