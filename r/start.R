# Check if we are in the project folder.
stopifnot(file.exists("./data/upstream-inputs/foray-dataset-non-sensitive.Rdata"))

# Clear memory
rm(list = ls())

# Libraries
library(stats)
library(data.table)
library(tidyverse)
library(lubridate)
library(readxl)
library(magrittr)

# Custom function
save_graph <- function (file_path) {
  ggsave(file_path, width = 16.1/2, height = 10/2, units = "in")
}