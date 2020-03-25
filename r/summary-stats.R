############################################
# Start
############################################
sink("./log/summary-stats.R.log", split = TRUE)
source("./r/start.R")
library(forcats)
library(openxlsx)

# Load data
load("./data/intermediate/foray-dataset-with-fitted-probability.Rdata")

# Save number of experiments before chucking the ones that didn't pass the audit
n_experiments_after_cleaning <- 
  foray_dataset_with_fitted_probability %>%
  select(ExperimentStepId) %>%
  unique() %>%
  nrow()

# Save number of audited experiments
df <- foray_dataset_with_fitted_probability %>%
  select(ExperimentStepId, ProbabilityLegit) %>%
  na.omit() %>%
  unique()

n_audited_experiments <- nrow(df)

# Drop experiments that are flukes in the audit for the rest of the analysis.
foray_dataset_with_fitted_probability %<>%
  filter(ProbabilityLegitFitted > 0)

n_experiments_after_cleaning_and_audit <- 
  foray_dataset_with_fitted_probability %>%
  select(ExperimentStepId) %>%
  unique() %>%
  nrow()

# Keep only metrics for our analysis
# Rename metric names and order them so that the summary stats table come out right
# Only keep observations with positive legit probability.
foray_dataset_with_fitted_probability %<>%
  filter(Metric %in%
           c("page click rate", "queries per user", "quickback rate",
             "session success rate", "sessions per user", "time to success")) %>%
  mutate(Metric = recode(Metric,
                         `session success rate` = "session success rate",
                         `page click rate` = "short-run metric 1",
                         `quickback rate` = "short-run metric 2",
                         `time to success` = "short-run metric 3",
                         `queries per user` = "long-run metric 1",
                         `sessions per user` = "long-run metric 2"),
         Metric = parse_factor(Metric,
                               levels = c("session success rate",
                                          "short-run metric 1",
                                          "short-run metric 2",
                                          "short-run metric 3",
                                          "long-run metric 1",
                                          "long-run metric 2")))
  


############################################
# 1: summary stats table
############################################
# Prepare data for two tables: ii = 1 for all experiments in the data, and ii = 2 for guaranteed legit experiments.
for (ii in c(1, 2)) {
  if (ii == 1) {
    dataset_to_use <- foray_dataset_with_fitted_probability
  } else {
    dataset_to_use <- foray_dataset_with_fitted_probability %>%
      filter(ProbabilityLegitFitted == 1)
  }
  
  # Experiment-level variables
  df <- dataset_to_use %>%
    select(N, AnalysisRange, ProbabilityLegitFitted)
  
  table_experiment_level <- do.call(data.frame, 
                                    list(mean = apply(df, 2, mean),
                                         min = apply(df, 2, min),
                                         max = apply(df, 2, max),
                                         sd = apply(df, 2, sd)
                                    ))
  
  # Experiment-metric level variables
  table_experiment_metric_level_delta <-
    dataset_to_use %>%
    select(Metric, Delta) %>%
    group_by(Metric) %>%
    summarize(mean = mean(Delta),
              min = min(Delta),
              max = max(Delta),
              sd  = sd(Delta),
              iqr = IQR(Delta))
  
  table_experiment_metric_level_error <- 
    dataset_to_use %>%
    select(Metric, StdErrorDelta) %>%
    group_by(Metric) %>%
    summarize(mean = mean(StdErrorDelta),
              min = min(StdErrorDelta),
              max = max(StdErrorDelta),
              sd  = sd(StdErrorDelta))
  
  # Count observations
  n_observations <- 
    dataset_to_use %>%
    select(ExperimentStepId) %>%
    unique() %>%
    nrow() %>%
    format(nsmall = 0, big.mark = ",")

  
  # Save table input to excel
  if (ii == 1) {
    spreadsheet_path <- "./data/table-construction-data/input-summary-stats-table-all.xlsx"
  } else {
    spreadsheet_path <- "./data/table-construction-data/input-summary-stats-table-confirmed.xlsx"
  }
  write.xlsx(list("experiment-level" = table_experiment_level,
                  "delta" = table_experiment_metric_level_delta,
                  "se" = table_experiment_metric_level_error,
                  "n" = n_observations),
             file = spreadsheet_path,
             rowNames = TRUE)
}


############################################
#2 Histogram
############################################
df <- foray_dataset_with_fitted_probability

ggplot(data = df %>% filter(Metric == "session success rate"), aes(x = Delta)) +
  geom_histogram(binwidth = 0.01) +
  xlim(-0.3, 0.3) +
  xlab("Sample deltas (session success)") +
  ylab("Count")
save_graph("./output/figures/histogram-session-success.pdf")

ggplot(data = df, aes(x = Delta)) +
  geom_histogram(binwidth = 0.01) +
  xlim(-0.5, 0.5) +
  xlab("Sample deltas (session success)") +
  ylab("Count") +
  facet_wrap(~Metric)
save_graph("./output/figures/histogram-all-metrics.pdf")


############################################
#3 Log log plot
############################################
# Plot parameters
n_observations_to_calculate_tails <- 15
n_observations_to_plot <- 150

# Set seed because of the simulation
set.seed(1)

# Make plot for 1 left tail, 2 right tail, 3 both tails, and 4 simulated data from normal distribution's both tails.
for (spec in c("left", "right", "both", "simulated")) {
  # Prep data
  df <- foray_dataset_with_fitted_probability
  
  # Select appropriate data
  switch(spec,
         left = {
           df %<>% filter(Delta < 0)
           x_axis_title <- "-Delta (log scale)"
           file_name <- "pareto-tail-left.pdf"
           },
         right = {
           df %<>% filter(Delta > 0)
           x_axis_title <- "Delta (log scale)"
           file_name <- "pareto-tail-right.pdf"
           },
         both = {
           file_name <- "pareto-tail-absolute-value.pdf"
           x_axis_title <- "Absolute value of delta (log scale)"
         },
         simulated = {
           df$Delta <- rnorm(n = nrow(df), sd = df$StdErrorDelta)
           x_axis_title <- "Absolute Value of Simulated Delta (log scale)"
           file_name <- "pareto-tail-simulated-normal-absolute-value.pdf"
           })
  
  ddf <-
    df %>%
    group_by(Metric) %>%
    mutate(
      abs_rank = -abs(Delta) %>% dense_rank(),
      log_abs_rank = abs_rank %>% log10(),
      log_abs_delta = abs(Delta) %>% log10()
      )
  
  # Calculate slopes and add to data
  metrics <- levels(ddf$Metric)
  slope <- vector(mode = "numeric", length = length(metrics))
  
  for (ii in 1:length(metrics)) {
    fit <- lm(data = ddf %>% filter(abs_rank <= n_observations_to_calculate_tails, Metric == metrics[[ii]]),
              formula = log_abs_rank ~ log_abs_delta)
    slope[[ii]] <- -fit$coefficients[[2]]
  }
  
  aux <- data.frame(Metric = metrics, Slope = slope) %>%
    mutate(Metric = factor(Metric, levels(ddf$Metric)))
  
  ddf <- left_join(ddf, aux)
  
  # Prep pretty titles that report slopes
  ddf %<>%
    mutate(FacetTitle = paste0(
      Metric,
      " (slope = ",
      Slope %>% round(2),
      ")"
    ))
  
  ddf$FacetTitle <- as.factor(ddf$FacetTitle)
  
  # levels(ddf$FacetTitle) <- levels(ddf$FacetTitle)[c(3, 4, 5, 6, 1, 2)]
  
  # Plot
  ggplot(data = ddf %>% filter(abs_rank <= n_observations_to_plot),
         aes(x = abs(Delta), y = abs_rank)) +
    geom_point() +
    geom_smooth(data = ddf %>% filter(abs_rank <= n_observations_to_calculate_tails),
                method = "lm") +
    facet_wrap(~FacetTitle) +
    ylab("Rank (log scale)") +
    xlab(x_axis_title) +
    scale_x_log10() + scale_y_log10()
  
  save_graph(paste0("./output/figures/", file_name))
}

############################################
#4: Prepare constants for paper
############################################
save_constant <- function(x, file_name) {
  write(x, file = paste0("./output/constants/", file_name, ".txt"))
}

# n_experiments_after_cleaning
n_experiments_after_cleaning %>%
  format(nsmall = 0, big.mark = ",") %>%
  save_constant("n-experiments-after-cleaning")

# n_audited_experiments
n_audited_experiments %>%
  format(nsmall = 0, big.mark = ",") %>%
  save_constant("n-audited-experiments")

# mean_percentage_valid_audit
mean_percentage_valid_audit <- 100*median(foray_dataset_with_fitted_probability$ProbabilityLegitFitted)

mean_percentage_valid_audit %>%
  round() %>%
  format(nsmall = 0, big.mark = ",") %>%
  save_constant("mean-percentage-valid-audit")

# mean_std_error_ssr
foray_dataset_with_fitted_probability %>%
  filter(Metric == "session success rate") %>%
  .$StdErrorDelta %>%
  mean() %>%
  round(3) %>%
  save_constant("mean-std-error-ssr")

# std_dev_sample_ssr
foray_dataset_with_fitted_probability %>%
  filter(Metric == "session success rate") %>%
  .$Delta %>%
  sd() %>%
  round(3) %>%
  save_constant("std-dev-sample-ssr")


######################
# End
######################
sink()

