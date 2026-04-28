##########################################################################################
####                           AGRILUS OAK HOSTS ANALYSES                          ######
####                           ELVIRA HERNANDEZ GUTIERREZ                          ######
##########################################################################################

####        SET ENVIRONMENT                                                         #####
# Set directory and clean environment
# setwd("~/oak_analyses/01_original_models")         # Whatever your working directory is

# Load packages
library(ape)
library(ggtree)
library(RColorBrewer)
library(ggtreeExtra)
library(rgbif)
library(caper)
library(plyr)
library(dplyr)
library(CoordinateCleaner)
library(geodist)
library(mcpdiversity)
library(brms)
library(ggplot2)
# library(raster)
# library(geodata)
# library(ggbiplot)
# library(factoextra)
# library(doParallel)



#
##### XI)   LEAVE-ONE-OUT CHECKS                                                     #####
load("analysis/results/oak_models.RData", sep = "")

# A) FIRST WAY - CODE ONE INTERACTION AS A '0' AT A TIME
# This modified (using 'fitted' instead of 'predict') piece of code comes from:
# https://www.r-bloggers.com/2015/12/calculate-leave-one-out-prediction-for-glm/

# i) Run analysis
# NEEDS TO BE RUN WITH A JOBSCRIPT: 1 core, 15 GB, 116 array jobs; took ca. 2h
n <- as.numeric(commandArgs(trailingOnly = TRUE))
obs <- which(interaction_data$interaction == 1)

loo_fitted <- function(mod, obs, n) {
  ndat <- mod$data; i <- obs[n]
  print(interaction_data[i, 1:3])

  # If the current interaction == 1, then update the model by making it == 0,
  # and predict the value for the current observation
  ndat$interaction[i] <- 0
  return(data.frame(row.no = i,
                    prediction.logodds.loo = fitted(update(mod, newdata = ndat),
                                                    mod$data[i, ],
                                                    scale = "linear")[, 1],
                    row.names = NULL))
}

predictions_loo0 <- loo_fitted(oak_mod043, obs, n)

# write.table(predictions_loo0,
#             paste("analysis/results/loo/predictions_loo0", n, ".tsv", sep = ""),
#             row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# END OF JOBSCRIPT

# ii) Merge files
loo0_pred_loods <- lapply(paste0("analysis/results/loo/",
                                 list.files(path = "analysis/results/loo")),
                          read.table, header = TRUE, sep = "\t")
loo0_pred_loods <- do.call(rbind.data.frame, loo0_pred_loods)
loo0_pred_loods <- loo0_pred_loods[order(loo0_pred_loods$row.no),]

# write.table(loo0_pred_loods, "analysis/results/predictions_loo0.tsv",
#             row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


# ii) Compare results against original predictions
loo0_pred_loods <- read.table("analysis/results/predictions_loo0.tsv",
                              sep = "\t", header = T)$prediction.logodds.loo
predictions1 <- subset(predictions, host.status == 1)
interaction_data1 <- subset(interaction_data, interaction == 1)
predictions_loo1 <- data.frame(agrilus.sp = interaction_data1$agrilus.sp,
                               quercus.sp = interaction_data1$quercus.sp,
                               interaction = as.factor(interaction_data1$interaction),
                               prediction.logodds.orig = predictions1$prediction.logodds,
                               prediction.logodds.loo0 = loo0_pred_loods)

# Check the differences between the "loo0" and the "original" predictions
plot(loo0_pred_loods ~ predictions1$prediction.logodds)
abline(a = 0, b = 1, col = "steelblue")
# identify(y = loo0_pred_loods, x = predictions$prediction.logodds)


# iii) Plot results (only for known hosts, i.e., 1s)
# LOO0 vs. original predicitons
predictions_loo1$thr <- "no"
for (i in 1:nrow(predictions_loo1)) {
  if (predictions_loo1$prediction.logodds.orig[i] >= cutoff[2] &
      predictions_loo1$prediction.logodds.loo0[i] < cutoff[2]) {
    predictions_loo1$thr[i] <- "yes"
  }
}

ggplot(predictions_loo1,
       aes(x = prediction.logodds.orig, y = prediction.logodds.loo0,
           col = thr)) +
  geom_abline(linetype = "dashed", col = "darkgrey", linewidth = 1.5) +
  geom_point() +
  scale_color_manual(values = c("black", "coral")) +
  xlab("Original predictions") + ylab("'One 0 in turn' predictions") +
  ggtitle("Quercus sp. - Agrilus sp.") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

predictions_loo1[predictions_loo1$thr == "yes",]
predictions_loo1$thr <- NULL


# iv) Check new rate of true positives
# Create confidence matrix
conf_matrix <- table(ifelse(predictions_loo1$prediction.logodds.loo0 > cutoff[2], 1, 0),
                     predictions_loo1$interaction)

# Rate of true positives = 83.62%
round(conf_matrix[2]/(conf_matrix[1] + conf_matrix[2])*100, 2)
rm(conf_matrix, predictions1, interaction_data1)



# B) SECOND WAY - LEAVE ONE INTERACTION OUT AT A TIME
# i) Compute loo predictions
# I believe loo_predict() is the equivalent to predict(), and loo_linpred() to
# fitted(method = "linear")
# Each one takes <5 min to run
loo_pred_lodds <- loo_linpred(oak_mod043, type = "mean")
loo_pred_prob <- loo_predict(oak_mod043, type = "mean")


# ii) Compare against other method
predictions_loo1$prediction.logodds.loo <- loo_pred_lodds[which(predictions$host.status
                                                                == 1)]
predictions_loo1$thr <- "no"
for (i in 1:nrow(predictions_loo1)) {
  if (predictions_loo1$prediction.logodds.loo[i] >= cutoff[2] &
      predictions_loo1$prediction.logodds.loo0[i] < cutoff[2]) {
    predictions_loo1$thr[i] <- "loo"
  } else if (predictions_loo1$prediction.logodds.loo0[i] >= cutoff[2] &
        predictions_loo1$prediction.logodds.loo[i] < cutoff[2]) {
    predictions_loo1$thr[i] <- "loo0"
  }
}

ggplot(predictions_loo1,
       aes(x = prediction.logodds.loo0, y = prediction.logodds.loo,
           col = thr)) +
  geom_abline(linetype = "dashed", col = "darkgrey", linewidth = 1.5) +
  geom_point() +
  scale_color_manual(values = c("black", "coral")) +
  xlab("'One 0 in turn' predictions") + ylab("LOO predictions") +
  ggtitle("Quercus sp. - Agrilus sp.") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

predictions_loo1$thr <- NULL


# iii) Compare against original predictions
predictions_loo <- data.frame(agrilus.sp = interaction_data$agrilus.sp,
                              quercus.sp = interaction_data$quercus.sp,
                              interaction = as.factor(interaction_data$interaction),
                              prediction.logodds.orig = predictions$prediction.logodds,
                              prediction.logodds.loo = loo_pred_lodds)

# Check the differences between the "loo" and the "original" predictions
plot(loo_pred_lodds ~ predictions$prediction.logodds)
abline(a = 0, b = 1, col = "steelblue")


# iv) Plot results
# Colour by interaction
ggplot(predictions_loo,
       aes(x = prediction.logodds.orig, y = prediction.logodds.loo, col = interaction)) +
  geom_abline(linetype = "dashed", col = "darkgrey") +
  geom_point() +
  xlab("Original model") + ylab("Leave One Out ('automatic')") +
  scale_color_manual(values = c("black", "coral")) +
  ggtitle("Quercus sp. - Agrilus sp.") +
  theme(plot.title = element_text(hjust = 0.5))

# Only plot 1s, colouring those points that have moved under the threshold
predictions_loo$thr <- "no"
for (i in 1:nrow(predictions_loo)) {
  if (predictions_loo$prediction.logodds.orig[i] >= cutoff[2] &
      predictions_loo$prediction.logodds.loo[i] < cutoff[2]) {
    predictions_loo$thr[i] <- "yes"
  }
}

ggplot(predictions_loo[predictions_loo$interaction == 1, ],
       aes(x = prediction.logodds.orig, y = prediction.logodds.loo,
           col = thr)) +
  geom_abline(linetype = "dashed", col = "darkgrey", linewidth = 1.5) +
  geom_point() +
  scale_color_manual(values = c("black", "coral")) +
  xlab("Original model") + ylab("Leave One Out ('automatic')") +
  ggtitle("Quercus sp. - Agrilus sp.") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

predictions_loo[predictions_loo$interaction == 1 & predictions_loo$thr == "yes",]
predictions_loo$thr <- NULL


# v) Check new rate of true positives
# Create confidence matrix
conf_matrix <- table(ifelse(predictions_loo1$prediction.logodds.loo > cutoff[2], 1, 0),
                     predictions_loo1$interaction)

# Rate of true positives = 83.62%
round(conf_matrix[2]/(conf_matrix[1] + conf_matrix[2])*100, 2)

rm(conf_matrix, i)



# C) KFOLD ANALYSIS
# i) Compute kfold predictions
# Run 5-fold cross-validation, storing the fitted models (save_fits = TRUE) to
# be able to generate predictions
kfold <- brms::kfold(oak_mod043, K = 5, save_fits = TRUE)

# For each observation, use the model where that observation was held out, and
# return posterior draws of predictions (yrep)
kfold_preds <- brms::kfold_predict(kfold, method = "posterior_epred")

# colMean gives one predicted probability per observation (yrep = matrix of
# posterior predictions; rows = draws, columns = observations)
kfold_prob <- colMeans(kfold_preds$yrep)

# Transform to log-odds
kfold_lodds <- qlogis(kfold_prob)


# ii) Compare kfold vs loo (second-method) results (log-odds)
loo_lodds <- predictions_loo$prediction.logodds.loo

# Compare kfold vs loo results (log-odds)
length(loo_lodds) == length(kfold_lodds)

# Check correlation
cor(loo_lodds, kfold_lodds)
cor(loo_lodds[idx_rep_pos], kfold_lodds[idx_rep_pos])

# Accuracy (% 'True Positives' detected)
length(which(predictions$host.status == 1))

round(nrow(predictions[predictions$host.status == 1
                       & predictions$prediction.logodds > cutoff[2], ]) /
        nrow(predictions[predictions$host.status == 1, ]) * 100, 2)
round(nrow(predictions[predictions$host.status == 1
                       & loo_lodds > cutoff[2], ]) /
        nrow(predictions[predictions$host.status == 1, ]) * 100, 2)
round(nrow(predictions[predictions$host.status == 1
                       & kfold_lodds > cutoff[2], ]) /
        nrow(predictions[predictions$host.status == 1, ]) * 100, 2)

# % points over threshold
length(which(predictions$prediction.logodds > cutoff[2]))
round(length(which(predictions$prediction.logodds > cutoff[2])) /
        nrow(predictions) * 100, 2)

length(which(loo_lodds > cutoff[2]))
round(length(which(loo_lodds > cutoff[2])) / nrow(predictions) * 100, 2)

length(which(kfold_lodds > cutoff[2]))
round(length(which(kfold_lodds > cutoff[2])) / nrow(predictions) * 100, 2)

# No. predictions missing
nrow(predictions[predictions$prediction.logodds > cutoff[2]
                 & loo_lodds < cutoff[2], ])
nrow(predictions[predictions$prediction.logodds > cutoff[2]
                 & kfold_lodds < cutoff[2], ])


# iii) Plot comparison
# Index reported interactions
idx_rep_pos <- which(predictions$host.status == 1)
idx_rep_neg <- which(predictions$host.status != 1)

# Index predicted positive interactions
idx_pred_pos <- which(predictions$prediction.logodds > cutoff[2])

# Index intersections
idx_pos_all <- intersect(idx_rep_pos, idx_pred_pos)
idx_rep_pos_pred_neg <- setdiff(idx_rep_pos, idx_pred_pos)

idx_rep_neg_pred_pos <- intersect(idx_rep_neg, idx_pred_pos)
idx_neg_all <- setdiff(idx_rep_neg, idx_pred_pos)

# Plot negative interactions
plot(loo_lodds[idx_neg_all], kfold_lodds[idx_neg_all],
     xlim = c(min(loo_lodds), max(loo_lodds)),
     ylim = c(min(kfold_lodds), max(kfold_lodds)),
     col = rgb(0, 0, 0, 0.02),
     pch = 16,
     xlab = "LOO CV predicted log-odds",
     ylab = "5-fold CV predicted log-odds")

# Add negative interactions predicted as positive (triangle)
points(loo_lodds[idx_rep_neg_pred_pos], kfold_lodds[idx_rep_neg_pred_pos],
       col = rgb(0, 0, 0, 0.1),
       pch = 17)

# Add reported positive interactions (blue) but not predicted as such
points(loo_lodds[idx_rep_pos_pred_neg], kfold_lodds[idx_rep_pos_pred_neg],
       bg = rgb(70, 130, 180, 255, maxColorValue = 255),
       col = rgb(0, 0, 0, 0.5),
       pch = 21)

# Add reported positive interactions (blue) predicted as such (triangle)
points(loo_lodds[idx_pos_all], kfold_lodds[idx_pos_all],
       bg = rgb(70, 130, 180, 255, maxColorValue = 255),
       col = rgb(0, 0, 0, 0.5),
       pch = 24)

# Add 1:1 line, and thresholds
abline(0, 1, lty = 3)
abline(cutoff[2], 0, lty = 3)
abline(v = cutoff[2], lty = 3)

# save.image("analysis/results/oak_models.RData")



#