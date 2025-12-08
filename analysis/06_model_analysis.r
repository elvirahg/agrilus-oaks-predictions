#### SET ENVIRONMENT ####
# Custom functions
source("R/analysis_functions.r")
source("R/plot_functions.r")

# Libraries
library(brms)
library(ggplot2)

# Load brms models from previous script
oak_models <- readRDS("data/results/oak_models.rds")


#### BASIC EXPLORATION OF SELECTED MODEL ####
# Read interaction data frame
interaction_data <- read.table("data/results/interaction_data.tsv",
                               header = TRUE,
                               sep = "\t")

# This is model 46 which is the same model that was selected in the original
# analyses (model 43 there)

# The model is a binary logistic regression (link: mu = logit), i.e., response
# = log-odds of being a host. To convert to probability: exp(y)/(1+exp(y))

# Summary stats
oak_models$mod046$formula
summary(oak_models$mod046)

# Compare the observed outcome variable y to simulated datasets yrep from the
# posterior predictive distribution; mc-stan.org/rstanarm/reference/
pp_check(oak_models$mod046)
pp_check(oak_models$mod046, "stat")

# nspect chains: trace and density plots (b: population-level)
plot(oak_models$mod046, N = 2, ask = FALSE)


#### PREDICT INTERACTION STATUS ####
predictions <- data.frame(
  quercus_sp = interaction_data$quercus_sp,
  agrilus_sp = interaction_data$agrilus_sp,
  host_status = factor(interaction_data$interaction),
  prediction_prob = fitted(oak_models$mod046, scale = "response")[, "Estimate"],
  prediction_lodds  = fitted(oak_models$mod046, scale = "linear")[, "Estimate"]
)

# Simple initial plot
par(mfrow = c(1, 2))
plot(prediction_prob ~ host_status, data = predictions,
     xlab = "Host status", ylab = "Predicted probability",
     col = c("steelblue", "coral"))
plot(prediction_lodds ~ host_status, data = predictions,
     xlab = "Host status", ylab = "Predicted log-odds",
     col = c("steelblue", "coral"))
par(mfrow = c(1, 1))

# Write predictions to file
# write.table(predictions,
#             "data/results/predictions.tsv",
#             row.names = FALSE,
#             col.names = TRUE,
#             sep = "\t",
#             quote = FALSE)


#### CHECK EFFECT OF DIFFERENT VARIABLES ####
comparisons <- data.frame(
  host_status = factor(interaction_data$interaction),
  oak = factor(interaction_data$quercus_sp),
  mod046 = predictions$prediction_lodds,
  mod003 = fitted(oak_models$mod003,
                  scale = "linear")[, "Estimate"],  # geo dist
  mod005 = fitted(oak_models$mod005,
                  scale = "linear")[, "Estimate"],  # phy dist
  mod008 = fitted(oak_odel$mod008,
                  scale = "linear")[, "Estimate"],  # (quercus | phylo)
  mod010 = fitted(oak_models$mod010,
                  scale = "linear")[, "Estimate"],  # geo dist + phy dist
  mod022 = fitted(oak_models$mod022,
                  scale = "linear")[, "Estimate"],  # geo dist + (quercus | phylo)
  mod029 = fitted(oak_models$mod029,
                  scale = "linear")[, "Estimate"]   # phy dist + (quercus | phylo)
)

# Plot comparisons of predictions for the full model vs the other models
plot_model_comparison(
  df = comparisons,
  full = "mod046",
  other = "mod003",
  xlab_text = "Geographic distance effect (model 3)"
)

plot_model_comparison(
  df = comparisons,
  full = "mod046",
  other = "mod005",
  xlab_text = "Phylogenetic distance fixed effect (model 5)"
)

plot_model_comparison(
  df = comparisons,
  full = "mod046",
  other = "mod008",
  xlab_text = "Phylogenetic random effect (model 8)"
)

plot_model_comparison(
  df = comparisons,
  full = "mod046",
  other = "mod010",
  xlab_text = "Geographic and fixed phylogenetic effect (model 10)"
)

plot_model_comparison(
  df = comparisons,
  full = "mod046",
  other = "mod022",
  xlab_text = "Geographic and random phylogenetic effect (model 22)"
)

plot_model_comparison(
  df = comparisons,
  full = "mod046",
  other = "mod029",
  xlab_text = "Fixed and random phylogenetic effects (model 29)"
)


#### FIND A PREDICTION THRESHOLD ####
# Compute prediction thresholds
# https://stackoverflow.com/questions/23240182/deciding-threshold-for-glm-
# logistic-regression-model-in-r
# https://statinfer.com/203-4-2-calculating-sensitivity-and-specificity-in-r/
thresholds_df <- compute_threshold_metrics(
  pred_values = predictions$prediction_lodds,
  obs_values = predictions$host_status
)

# NB, threshold using 'baseline' log-odds of sth being a host, given
# 7552 total observations & 116 'positive' observations:
# logit(prob) = log(prob/(1-prob)), i.e. log((116/7552)/(1-(116/7552)))
# Which is very close to the inferred threshold!
thr_intercepts <- compute_intercept(thr_df = thresholds_df)

# Plot sensitivity and specificity against threshold values with intersection
plot_sens_spec_threshold(thr_df = thresholds_df,
                         thr_ints = thr_intercepts)

# Plot ROC curve to assess the model's ability to differentiate between known
# hosts and alleged non-hosts under different probability thresholds
# To visualise ROC: plot(thresholds_df$sens~c(1-thresholds_df$spec))
ggplot(thresholds_df, aes(fpr, sensitivity)) +
  geom_line() +
  labs(x = "FPR", y = "TPR") +
  theme_minimal()

# Compute AUC value (0.97)
DescTools::AUC(thresholds_df$fpr, thresholds_df$sens)


#### PLOT RESULTS ####
# Genral violin plot
plot_predictions_violin(predictions = predictions,
                        threshold = thr_intercepts["threshold"])

median(predictions[predictions$host_status == 1, ]$prediction_lodds)
median(predictions[predictions$host_status == 0, ]$prediction_lodds)

# Add binary predictions to predictions dataframe
predictions$prediction_binary <- factor(
  ifelse(predictions$prediction_lodds > thr_intercepts["threshold"], 1, 0),
  levels = c(0, 1)
)

# Plot general binary predictions
plot_fourfold(observations = predictions$host_status,
              predictions = predictions$prediction_binary,
              print_metrics = TRUE)

# Plot by oak, all oaks
plot_predictions_by_species(data = predictions,
                            species_col = "quercus_sp",
                            threshold = thr_intercepts["threshold"],
                            pred_var = "prediction_lodds",
                            hosts_only = FALSE,
                            ncol = 16,
                            point_size = 1,
                            title = NULL,
                            italic_text = TRUE)

# Plot by oak, hosts only
plot_predictions_by_species(data = predictions,
                            species_col = "quercus_sp",
                            threshold = thr_intercepts["threshold"],
                            pred_var = "prediction_lodds",
                            hosts_only = TRUE,
                            ncol = 16,
                            point_size = 1,
                            title = NULL,
                            italic_text = TRUE)

# Plot by Agrilus species
plot_predictions_by_species(data = predictions,
                            species_col = "agrilus_sp",
                            threshold = thr_intercepts["threshold"],
                            pred_var = "prediction_lodds",
                            hosts_only = FALSE,
                            ncol = 8,
                            point_size = 1,
                            title = NULL,
                            italic_text = TRUE)

# Plot for a specific Agrilus species, labelling novel predicted hosts
plot_highlighted_predictions(
  data = subset(predictions, agrilus_sp == "Agrilus angustulus"),
  hosted_col = "agrilus_sp",
  host_col = "quercus_sp",
  threshold = thr_intercepts["threshold"]
)
