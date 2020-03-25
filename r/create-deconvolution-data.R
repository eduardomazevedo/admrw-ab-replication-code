################################################
# Start
sink("./log/create-deconvolution-data.R.log", split = TRUE)
source("./r/start.R")
load("./data/intermediate/foray-dataset-with-fitted-probability.Rdata")


## Part 1: aggregated data
df <- foray_dataset_with_fitted_probability

# Filter
df %<>%
  filter(ProbabilityLegitFitted > 0) %>%
  select(Metric,
         N, Delta, StdErrorDelta, sigma,
         ProbabilityLegitFitted) %>%
  na.omit()


# Save
write_csv(df, "./data/intermediate/convolution-data.csv")


## Part 2: disaggregated data
df <- foray_dataset_with_fitted_probability %>%
  filter(ProbabilityLegitFitted > 0,
         Metric == "session success rate")

# Create data frame where we will stack all the specifications.
ddf <- tribble()


# Specs 1 to 3: date ranges
first_date = min(df$AnalysisStartDate)
last_date = max(df$AnalysisStartDate)
l <- floor((last_date - first_date) / 3) + 1  # About 1/3 of the date range (324 days)

for (ii in 1:3) {
  start_date = first_date + (ii - 1) * l
  end_date = start_date + l - ddays(1)
  if (ii == 3)
    end_date = end_date + ddays(1)
  spec_description = paste0("Date: ", as.character(start_date), " to ", as.character(end_date))
  
  ddf <- rbind(ddf,
               df %>%
                 filter(AnalysisStartDate >= start_date,
                        AnalysisStartDate <= end_date) %>%
                 select(Metric,
                        N, Delta, StdErrorDelta, sigma,
                        ProbabilityLegitFitted) %>%
                 na.omit() %>%
                 mutate(Spec = ii,
                        SpecDescription = spec_description))
}

# Check that we partitioned the data
stopifnot(nrow(ddf) == nrow(df))


# Spects 4 to 6: budget areas
df %<>%
  mutate(BudgetArea = as.character(BudgetArea)) %>%
  separate(BudgetArea, into = c("BudgetArea"), sep = " - ") %>%
  separate(BudgetArea, into = c("BudgetArea"), sep = " ") %>%
  separate(BudgetArea, into = c("BudgetArea"), sep = "-") %>%
  mutate(BudgetArea = if_else(BudgetArea == "coreux", "ux", BudgetArea))

list_of_areas <- unique(df$BudgetArea)
n_areas <- length(list_of_areas)
n_specs_for_time <- max(ddf$Spec)

for (ii in 1:n_areas) {
  spec_description <- paste0("Budget area: ", list_of_areas[[ii]])
  spec_number <- n_specs_for_time + ii
  
  ddf <- rbind(ddf,
               df %>%
                 filter(BudgetArea == list_of_areas[[ii]]) %>%
                 select(Metric,
                        N, Delta, StdErrorDelta, sigma,
                        ProbabilityLegitFitted) %>%
                 na.omit() %>%
                 mutate(Spec = spec_number,
                        SpecDescription = spec_description))
}

# Check that we have two partitions
stopifnot(nrow(ddf) == 2* nrow(df))


# Specs 7 and 8: one week experiment versus longer time periods.
spec_number = 7
spec_description = "Length: one week"
ddf %<>% rbind(
  df %>%
    filter(AnalysisRange == 7) %>%
    select(Metric,
           N, Delta, StdErrorDelta, sigma,
           ProbabilityLegitFitted) %>%
    na.omit() %>%
    mutate(Spec = spec_number,
           SpecDescription = spec_description))

spec_number <- 8
spec_description = "Length: over one week"
ddf %<>% rbind(
  df %>%
    filter(AnalysisRange > 7) %>%
    select(Metric,
           N, Delta, StdErrorDelta, sigma,
           ProbabilityLegitFitted) %>%
    na.omit() %>%
    mutate(Spec = spec_number,
           SpecDescription = spec_description))

# Check that we have three partitions
stopifnot(nrow(ddf) == 3 * nrow(df))


# Spec 9: exact multiples of one week
spec_number <- 9
spec_description = "Length: multiple of one week"
ddf %<>% rbind(
  df %>%
    filter((AnalysisRange %% 7) == 0) %>%
    select(Metric,
           N, Delta, StdErrorDelta, sigma,
           ProbabilityLegitFitted) %>%
    na.omit() %>%
    mutate(Spec = spec_number,
           SpecDescription = spec_description))

n_multiples_week <- df %>%
  filter((AnalysisRange %% 7) == 0) %>%
  nrow()

# Specs 10 to 12: quartiles of sample size
tercile_vector <- ntile(df$N, n = 3)

for(n_tercile in 1:4) {
  spec_number = 9 + n_tercile
  spec_description = paste0("Sample size tercile: ", n_tercile)
  ddf %<>% rbind(
    df %>%
      filter(tercile_vector == n_tercile) %>%
      select(Metric,
             N, Delta, StdErrorDelta, sigma,
             ProbabilityLegitFitted) %>%
      na.omit() %>%
      mutate(Spec = spec_number,
             SpecDescription = spec_description))
}

# Check that we have four partitions
stopifnot(nrow(ddf) == 4 * nrow(df) + n_multiples_week)


# Save disaggregated MLE data
write_csv(ddf, "./data/intermediate/convolution-data-disaggregated.csv")

################################################
# End
sink()

