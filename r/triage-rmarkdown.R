##===============================================================##
## Triage Analysis                                               ##
## RA A/B Testing                                                ##
## Author: Felipe Flores Golfin                        June 2018 ##
##===============================================================##

library(rmarkdown)
library(ggplot2)
library(gridExtra)

# Clear working space
rm(list = ls())

path = getwd()

render("./r/triage-analysis.rmd", output_dir = "./log", knit_root_dir = path)