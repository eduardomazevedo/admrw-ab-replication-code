################################################
# Start
sink("./log/clean.R.log", split = TRUE)
library(glmnet)
source("./r/start.R")


# Load data
load("./data/intermediate/data-for-probability-legit-fitting.Rdata")


# Keep data for SSR and fetures for the model that fits legit probability.
df <- foray_dataset_clean %>%
  filter(Metric == "session success rate") %>%
  select(Delta, ProbabilityLegit) %>%
  arrange(Delta) %>%
  filter(abs(Delta) < 1) %>%
  na.omit() %>%
  mutate(AbsDelta = abs(Delta),
         DeltaSquared = Delta^2)


################################################
# Calculate best fit with a LASSO model. TODO: Pepe has some idea to improve this.
x <- model.matrix(ProbabilityLegit ~ ., data = df)
y <- df$ProbabilityLegit

fit <- glmnet(x, y, alpha = 1)
cvfit <- cv.glmnet(x, y, alpha = 1, nfolds = length(y))

pdf("./output/figures/legit-fit-lasso-path.pdf",
    width = 8.05,
    height = 5)
plotmo::plot_glmnet(fit)
dev.off()

pdf("./output/figures/legit-fit-cross-validation.pdf",
    width = 8.05,
    height = 5)
plot(cvfit)
dev.off()

beta <- coef(cvfit, s = "lambda.1se")


################################################
# Graphical representation of lack of relationship between Delta and being legit.
ggplot(data = df %>% na.omit(),
     aes(x = Delta, y = ProbabilityLegit)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Measured delta in success rate") + ylab("Probability valid") +
  theme_minimal()

save_graph("./output/figures/legit-is-uncorrelated-with-effect-size.pdf")


################################################
# Fit legit probabilities for each experiment
df <- foray_dataset_clean %>%
  filter(Metric == "session success rate") %>%
  select(ExperimentStepId, Delta, ProbabilityLegit) %>%
  mutate(AbsDelta = abs(Delta),
         DeltaSquared = Delta^2) %>%
  mutate(ProbabilityLegit = 0)

x <- model.matrix(ProbabilityLegit ~ ., data = df %>% select(-ExperimentStepId))

y <- predict(cvfit, newx = x, s = "lambda.1se")

df$ProbabilityLegitFitted <- y

df %<>%
  select(ExperimentStepId, ProbabilityLegitFitted)


################################################
# Create fitted probability in the main dataset, with the audit value when available
foray_dataset_with_fitted_probability <-
  foray_dataset_clean %>%
  left_join(df, by = "ExperimentStepId") %>%
  mutate(ProbabilityLegitFitted = ifelse(is.na(ProbabilityLegit), ProbabilityLegitFitted, ProbabilityLegit))


################################################
# Save legit probabilities
# Spread probability among duplicated experiments
save(foray_dataset_with_fitted_probability,
  file = "./data/intermediate/foray-dataset-with-fitted-probability.Rdata")


################################################
# End
sink()