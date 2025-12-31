#### SET ENVIRONMENT ####
# Custom functions
source("R/model_functions.r")
source("R/analysis_functions.r")
source("R/plot_functions.r")

# Libraries
library(brms)
library(foreach)
library(ggplot2)

# Seed (not needed for actual analyses, only if exact reproducibility of
# results is desired)
set.seed(24601)

# Load brms models from previous script
oak_models <- readRDS("data/results/oak_models.rds")


#### LEAVE-ONE-OUT CHECKS: CODE INTERACTIONS AS '0', ONE AT A TIME ####
# Load original predictions
predictions <- read.table("data/results/predictions.tsv",
                          sep = "\t",
                          header = TRUE)

# Calculate threshold
thresholds_df <- compute_threshold_metrics(
  pred_values = predictions$prediction_lodds,
  obs_values = predictions$host_status
)

thr_intercepts <- compute_intercept(thr_df = thresholds_df)

# Compute one_as_zero analysis
# See www.r-bloggers.com/2015/12/calculate-leave-one-out-prediction-for-glm/
preds_one_as_zero <- predict_as_zero_loo(model = oak_models$mod046,
                                         interaction_col = "interaction",
                                         cores = 5,
                                         iter = 1000,
                                         warmup = 500,
                                         chains = 1,
                                         recompile = FALSE,
                                         scale = "linear")

# Write results to file
# write.table(preds_one_as_zero,
#             "data/results/predictions_loo_one_as_zero.tsv",
#             row.names = FALSE,
#             col.names = TRUE,
#             sep = "\t",
#             quote = FALSE)

# preds_one_as_zero <- read.table("data/results/predictions_loo_one_as_zero.tsv",
#                                 sep = "\t",
#                                 header = TRUE)

# Plot prediction changes, to check whether there are any cases where the
# orignal prediction value was > threshold but is now < threshold
preds_original_pos <- subset(predictions, host_status == 1)[["prediction_lodds"]]

plot_prediction_change(pred_original = preds_original_pos,
                       pred_new = preds_one_as_zero$pred,
                       threshold = thr_intercepts["threshold"],
                       print_change = TRUE,
                       x_lab = "Original predictions (log-odds)",
                       y_lab = "'One-as-zero' predictions (log-odds)",
                       title = "Quercus sp. - Agrilus sp.")

# Calculate new TPR (%)
one_as_zero_pred_hosts <- ifelse(preds_one_as_zero$pred > thr_intercepts["threshold"], 1, 0)
tpr(predictions = one_as_zero_pred_hosts)


#### LEAVE-ONE-OUT CHECKS: LEAVE ONE INTERACTION OUT AT A TIME ####
# Compute predictions: loo_predict() is the equivalent to predict(),
# and loo_linpred() to fitted(method = "linear")
preds_loo <- loo_linpred(oak_models$mod046, type = "mean")[, 1]
preds_loo_pos <- preds_loo[which(predictions$host_status == 1)]

# Compare against the previous method
plot_prediction_change(pred_original = preds_one_as_zero$pred,
                       pred_new = preds_loo_pos,
                       threshold = thr_intercepts["threshold"],
                       print_change = TRUE,
                       col_by = "pred_change",
                       x_lab = "'One-as-zero' predictions (log-odds)",
                       y_lab = "LOO predictions (log-odds)",
                       title = "Quercus sp. - Agrilus sp.")

# Compare against all original predictions, colouring by obs_vals
plot_prediction_change(pred_original = predictions$prediction_lodds,
                       pred_new = preds_loo,
                       threshold = thr_intercepts["threshold"],
                       obs_vals = as.numeric(predictions$host_status),
                       print_change = TRUE,
                       col_by = "obs_vals",
                       x_lab = "Original predictions (log-odds)",
                       y_lab = "LOO predictions (log-odds)",
                       title = "Quercus sp. - Agrilus sp.")

# Compare against positive original predictions, colouring by pred_change
plot_prediction_change(pred_original = preds_original_pos,
                       pred_new = preds_loo_pos,
                       threshold = thr_intercepts["threshold"],
                       print_change = TRUE,
                       col_by = "pred_change",
                       x_lab = "Original predictions (log-odds)",
                       y_lab = "LOO predictions (log-odds)",
                       title = "Quercus sp. - Agrilus sp.")

# Calculate new TPR (%)
loo_pred_hosts <- ifelse(preds_loo > thr_intercepts["threshold"], 1, 0)
tpr(predictions = loo_pred_hosts,
    observations = predictions$host_status)
