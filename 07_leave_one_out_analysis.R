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
load("results/oak_models.RData", sep = "")

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
#             paste("results/loo/predictions_loo0", n, ".tsv", sep = ""),
#             row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# END OF JOBSCRIPT

# ii) Merge files
loo0_pred_loods <- lapply(paste0("results/loo/", list.files(path = "results/loo")),
                          read.table, header = TRUE, sep = "\t")
loo0_pred_loods <- do.call(rbind.data.frame, loo0_pred_loods)
loo0_pred_loods <- loo0_pred_loods[order(loo0_pred_loods$row.no),]

# write.table(loo0_pred_loods, "results/predictions_loo0.tsv",
#             row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


# ii) Compare results against original predictions
loo0_pred_loods <- read.table("results/predictions_loo0.tsv",
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

# save.image("results/oak_models.RData")



#