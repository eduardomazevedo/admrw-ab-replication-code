################################################
# Start
sink("./log/clean.R.log", split = TRUE)
source("./r/start.R")


# Load upstream dataset and load Surya's audit data
load("./data/upstream-inputs/foray-dataset-non-sensitive.Rdata")


# TODO: produce summary stats from data cleaning for the paper.

# Filter sample
foray_dataset_non_sensitive %<>%
  filter(NTreatments == 1,
         ExperimentHasTriggerFilterFlag == FALSE)

# Throw out AuditResult categorical variable because it does not use the latest version of the survey (see Codebook for explanation)
foray_dataset_non_sensitive %<>% select(-AuditResult)

# Save
foray_dataset_clean <- foray_dataset_non_sensitive
save(foray_dataset_clean, file = "./data/intermediate/data-for-probability-legit-fitting.Rdata")


################################################
# End
sink()
