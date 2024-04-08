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
#####      EXPLORE MODEL 43                                                         #####
load("results/oak_models.RData")

# This is a binary logistic regression (link: mu = logit), i.e., the response ==
# log-odds of something being a host. To convert it to probability: exp(y)/(1+exp(y))
oak_mod043$formula

# NB, to chech variable names: variables(oak_mod043)[1:5]

# A) SUMMARY AND DIAGNOSTICS
# i) First, have a look at the summary
summary(oak_mod043)


# ii) Diagnostics
# Compare the observed outcome variable y to simulated datasets yrep from the
# posterior predictive distribution (looks good).
# yrep: predicted scores from single replication; y: observed outcome
pp_check(oak_mod043)                # looks OK
pp_check(oak_mod043, "stat")        # mc-stan.org/rstanarm/reference/

# Diagnostic plots (Inspect chains - Trace and Density Plots)
# b: population-level
plot(oak_mod043, N = 2, ask = FALSE)


# B) PREDICT INTERACTION STATUS
# *fitted(, scale = ""response") (== predict()) uses the response scale, i.e., results
# after applying the inverse link function
# *fitted(, scale = "linear") return the value pre inverse link transformation

# Turn model predictions into a DF
predictions <- data.frame(cbind(interaction_data$quercus.sp,
                                interaction_data$agrilus.sp,
                                interaction_data$interaction,
                                fitted(oak_mod043, scale = "response")[,1],
                                fitted(oak_mod043, scale = "linear")[,1]))
colnames(predictions) <- c("quercus.sp", "agrilus.sp", "host.status",
                           "prediction.prob", "prediction.logodds")
predictions$host.status <- as.factor(predictions$host.status)
predictions$prediction.prob <- as.numeric(predictions$prediction.prob)
predictions$prediction.logodds <- as.numeric(predictions$prediction.logodds)

# Make a very simple plot
plot(predictions$prediction.prob ~ predictions$host.status)
plot(predictions$prediction.logodds ~ predictions$host.status)



# C) CHECK EFFECT OF DIFFERENT VARIABLES
# NB, fitted(oak_mod10, s = "l") == fitted(oak_mod043, re_formula = NA, s = "l")
comparisons <- data.frame(
  host.status = as.factor(interaction_data$interaction),
  oak = as.factor(interaction_data$quercus.sp),
  mod043 = fitted(oak_mod043, scale = "linear")[,1],  ## full
  mod002 = fitted(oak_mod002, scale = "linear")[,1],  ## geo
  mod004 = fitted(oak_mod004, scale = "linear")[,1],  ## phy1
  mod008 = fitted(oak_mod008, scale = "linear")[,1],  ## (1 | gr(oak | phy2))
  mod010 = fitted(oak_mod010, scale = "linear")[,1],  ## geo + phy1
  mod018 = fitted(oak_mod018, scale = "linear")[,1],  ## geo + (1 | gr(oak | phy2))
  mod025 = fitted(oak_mod025, scale = "linear")[,1]   ## phy1 + (1 | gr(oak | phy2))
  )


# i) Full model vs. geographic-distance-only model
summary(lm(comparisons$mod043 ~ comparisons$mod002))

ggplot(data = comparisons,
  aes(y = mod043, x = mod002, col = host.status)) +
  scale_colour_manual(values = c("lightgrey", "steelblue")) +
  geom_point() +
  geom_abline(intercept = -0.70436, slope = 1.24781,
              colour = "black", linewidth = 2, linetype = 'dashed') +
  ylab("Full model (model 43)") +
  xlab("Geographic distance model (model 2)") +
  theme_bw() +
  theme(legend.position = "none")

# hist(interaction_data$min.dist.mean.norm, breaks = 100, ylim = c(0, 500))


# ii) Full model vs. phylogenetic-distance-(fixed)-only model
summary(lm(comparisons$mod043 ~ comparisons$mod004))

ggplot(data = comparisons,
       aes(y = mod043, x = mod004, col = host.status)) +
  scale_colour_manual(values = c("lightgrey", "steelblue")) +
  geom_point() +
  geom_abline(intercept = -0.77562, slope = 1.44191,
              colour = "black", linewidth = 2, linetype = 'dashed') +
  ylab("Full model (model 43)") +
  xlab("Phylogenetic distance (fixed) model (model 4)") +
  theme_bw() +
  theme(legend.position = "none")

# hist(interaction_data$phylo.dist.min, breaks = 100, ylim = c(0, 500))


# iii) Full model vs. phylogenetic-random-only model
summary(lm(comparisons$mod043 ~ comparisons$mod008))

ggplot(data = comparisons,
       aes(y = mod043, x = mod008, col = host.status)) +
  scale_colour_manual(values = c("lightgrey", "steelblue")) +
  geom_point() +
  geom_abline(intercept = -1.30016, slope = 1.00775,
              colour = "black", linewidth = 2, linetype = 'dashed') +
  ylab("Full model (model 43)") +
  xlab("Phylogenetic random model (model 8)") +
  theme_bw() +
  theme(legend.position = "none")


# iv) Full model vs. phylogenetic-only (both fixed + random) model
summary(lm(comparisons$mod043 ~ comparisons$mod025))

ggplot(data = comparisons,
       aes(y = mod043, x = mod025, col = host.status)) +
  scale_colour_manual(values = c("lightgrey", "steelblue")) +
  geom_point() +
  geom_abline(intercept = -1.093212, slope = 0.924603,
              colour = "black", linewidth = 2, linetype = 'dashed') +
  ylab("Full model (model 43)") +
  xlab(paste("Phylogenetic-only model\n",
             "(model 25, both fixed + random phylogenetic effects)")) +
  theme_bw() +
  theme(legend.position = "none")


# v) Full model vs. geographic + phylogenetic (fixed) distance model
summary(lm(comparisons$mod043 ~ comparisons$mod010))

# Clustered by oak species
ggplot(data = comparisons,
       aes(y = mod043, x = mod010, col = host.status)) +
  scale_colour_manual(values = c("lightgrey", "steelblue")) +
  geom_point() +
  geom_abline(intercept = -0.24889, slope = 1.30971,
              colour = "black", linewidth = 2, linetype = 'dashed') +
  ylab("Full model (model 43)") +
  xlab(paste("Geographic and phylogenetic (fixed) distance model\n",
             "(model 10, no phylogenetic random effect)")) +
  theme_bw() +
  theme(legend.position = "none")


# vi) Full model vs. geographic distance + phylogenetic random effect model
summary(lm(comparisons$mod043 ~ comparisons$mod018))

ggplot(data = comparisons,
       aes(y = mod043, x = mod018, col = host.status)) +
  scale_colour_manual(values = c("lightgrey", "steelblue")) +
  geom_point() +
  geom_abline(intercept = -1.093212, slope = 0.924603,
              colour = "black", linewidth = 2, linetype = 'dashed') +
  ylab("Full model (model 43)") +
  xlab(paste("Geographic distance and phylogenetic-random-effect model\n",
             "(model 10, no phylogenetic distance fixed effect)")) +
  theme_bw() +
  theme(legend.position = "none")



#
#####      FIND A PREDICTION THRESHOLD FOR MODEL 43                                 #####
# A) COMPUTE THRESHOLD
# i) Compute range of thresholds to test
# First, find the minimum and maximum predicted values
# Minimum predicted value (rounded to closest lower .0 or 0.5 decimal)
min_v <- min(predictions$prediction.logodds); min_v

# If round == floor, leave as round(); else, leave as round() - 0.5
if ((round(min_v)) == floor(min_v)) {
  min_v <- round(min_v)
} else {
  min_v <- round(min_v) - 0.5
}; min_v

# Maximum predicted value (rounded to closest upper .0 or 0.5 decimal)
max_v <- max(predictions$prediction.logodds); max_v

# If round == floor, leave as round() + 0.5; else, leave as round()
if (round(max_v) == floor(max_v)) {
  max_v <- round(max_v) + 0.5
} else {
  max_v <- round(max_v)
}; max_v

# Next, get a sequence of the threshold values we're going to be looking at (i.e.,
# values from the min. val to the max. val, with 0.5 increases)
thrs <- seq(from = min_v, to = max_v, by = 0.5); rm(min_v, max_v)


# ii) Compute threshold value
# See https://stackoverflow.com/questions/23240182/deciding-threshold-for-glm-logis
# tic-regression-model-in-r
# And https://statinfer.com/203-4-2-calculating-sensitivity-and-specificity-in-r/

# Initialise an empty vector for:
sens <- c()            # Sensitivity (i.e., True Positive Rate = TP/(TP+FN))
spec <- c()            # Specificity (i.e., True Negative Rate = TN/(TN+FP))
fpr <- c()             # FPR (i.e., False Positive Rate = FP/(FP + TN))

# Now, for every threshold value
for (thr in thrs) {
  # First, compute the Confusion Matrix (i.e., "a summary of prediction results on a
  # classification problem") under current threshold
  # Rows = Predicted values classified as 1s or 0s according to the current threshold
  # Cols = Observed values
  # Matrix will look like:
  #            0 (Obs)   1 (Obs)
  # 0 (Pred)   TN        FN
  # 1 (pred)   FP        TP
  conf_matrix <- table(ifelse(predictions$prediction.logodds > thr, 1, 0),   # rows
                       predictions$host.status)                              # cols
  # print(conf_matrix)

  # Then, compute sensitivity (TPR), specificity (TNR) and FPR
  # Let's deal with the exceptions first (i.e., if all predicted values are either 1
  # or 0 under the current threshold)
  if(nrow(conf_matrix) == 1) {
    # If all pred vals are == 1, then TPR (sens) and FPR = 1, and TNR (spec) = 0
    if(dimnames(conf_matrix)[[1]] == "1") {
      sens <- c(sens, 1); spec <- c(spec, 0); fpr <- c(fpr, 1)
      # If all pred vals are == 0, then TPR (sens) and FPR = 0, and TNR (spec) = 1
    } else {
      sens <- c(sens, 0); spec <- c(spec, 1); fpr <- c(fpr, 0)
    }

    # Now, if the current threshold divides predicted values into "1"s and "0"s
  } else {
    # Sensitivity (i.e., True Positive Rate) = TP/(TP+FN)
    sens <- c(sens, conf_matrix[2,2]/(conf_matrix[2,2] + conf_matrix[1,2]))

    # Specificity (i.e., True Negative Rate) = TN/(TN+FP)
    spec <- c(spec, conf_matrix[1,1]/(conf_matrix[1,1] + conf_matrix[2,1]))

    # False Positive Rate (FPR) = FP/(FP + TN)
    fpr <- c(fpr, conf_matrix[2,1]/(conf_matrix[2,1] + conf_matrix[1,1]))
  }
}

# Create a DF with all this info
thrs_df <- as.data.frame(cbind(thrs, sens, spec, fpr))
rm(sens, spec, thrs, fpr, thr, conf_matrix)


# iv) Plot results and find intercept
# First, plot sensitivity and specificity vs. threshold, to get a visual idea of where
# the cutoff point is
thrs_df %>% ggplot(aes(thrs, sens)) +
  geom_line(colour = "steelblue", linewidth = 2) +
  geom_line(aes(thrs, spec, col = "coral"), linewidth = 2) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Specificity")) +
  # xlim(-4.5, -4) + ylim(0.875, 0.925) +
  labs(x = "Threshold", y = "Sensitivity") +
  theme(axis.title.y = element_text(colour = "steelblue"),
        axis.title.y.right = element_text(colour = "coral"),
        axis.title = element_text(size = 12), # face = "bold"
        legend.position = "none")

# Now, find the intersect
# We saw above that the cutoff point is between -4.5 and -4.0 (as this is where the
# lines intersect), so now we can compute linear models for both the sensitivity and
# the specificity lines between these two points
lm_sens <- lm(thrs_df[thrs_df$thrs == -4.5 | thrs_df$thrs == -4.0, c(1)] ~
                thrs_df[thrs_df$thrs == -4.5 | thrs_df$thrs == -4.0, c(2)])
lm_spec <- lm(thrs_df[thrs_df$thrs == -4.5 | thrs_df$thrs == -4.0, c(1)] ~
                thrs_df[thrs_df$thrs == -4.5 | thrs_df$thrs == -4.0, c(3)])

# Compute the intercept
# See https://stackoverflow.com/questions/7114703/finding-where-two-linear-fits-inter
# sect-in-r
cm <- rbind(coef(lm_sens),coef(lm_spec))    # Coefficient matrix
cutoff <- c(-solve(cbind(cm[,2],-1))
            %*% cm[,1])
cutoff                                      # y (sens|spec) = 0.905; x (thr) = -4.154


# And finally, plot properly with cutoff points
thrs_df %>% ggplot(aes(thrs, sens)) +
  geom_line(colour = "steelblue", linewidth = 2) +
  geom_line(aes(thrs, spec, col = "coral"), linewidth = 2) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Specificity")) +
  geom_hline(yintercept = cutoff[1], linetype="dashed", color = "darkgrey") +
  geom_vline(xintercept = cutoff[2], linetype="dashed", color = "darkgrey") +
  labs(x = "Threshold", y = "Sensitivity") +
  theme(axis.title.y = element_text(colour = "steelblue"),
        axis.title.y.right = element_text(colour = "coral"),
        axis.title = element_text(size = 12), # face = "bold"
        legend.position = "none")

# Also have a look at the AUC plot, which can be used to assess the ability of the model
# to differentiate between known hosts and alleged non-hosts under different probability
# thresholds.
# ROC: plot(thrs_df$sens~c(1-thrs_df$spec))
thrs_df %>% ggplot(aes(fpr, sens)) +
  geom_line() +
  labs(x = 'FPR', y = "TPR") +
  theme(axis.title.y.right = element_text(colour = "red"), legend.position="none")

# Compute AUC value
DescTools::AUC(thrs_df$fpr, thrs_df$sens)    # 0.963
rm(thrs_df, lm_spec, lm_sens, cm)



#
#####      PLOT RESULTS (GENERAL)                                                   #####
# A) VIOLIN PLOT
ggplot(data = predictions,
       aes(y = prediction.logodds, x = host.status, fill = host.status)) +

  # Threshold using 'baseline' log-odds of sth being a host
  # NB, 7552 = no. observations & 116 == no. 'positive' observations
  # logit(p) = log(p/(1-p)), with p = prob
  # geom_hline(yintercept = log((116/7552)/(1-(116/7552))),
  #            linetype = "dashed", col = "black") +
  # New inferred threshold (very close to the one above!)
  geom_hline(yintercept = cutoff[2], linewidth = 1.5,
             linetype = "dotted", col = "maroon") +

  geom_violin(width = 1) +
  geom_boxplot(width = 0.1, fill = "grey") +

  scale_fill_manual(values = c("aquamarine2", "coral")) +
  ggtitle("Predicted interactions") +
  xlab("Host status") + ylab("Log-odds predicted host status") +

  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))

median(predictions[predictions$host.status == 1, ]$prediction.logodds)
median(predictions[predictions$host.status == 0, ]$prediction.logodds)



# B) BINARY PLOTS USING THRESHOLD
# i) Prepare data
predictions <- data.frame(predictions,
                          ifelse(predictions$prediction.logodds > cutoff[2], 1, 0))
colnames(predictions) <- c(colnames(predictions)[1:5], "prediction.cutoff")
predictions$prediction.cutoff <- as.factor(predictions$prediction.cutoff)

# Check confidence matrix (predictions = rows, observations = columns)
conf_matrix <- table(ifelse(predictions$prediction.logodds > cutoff[2], 1, 0),
                     predictions$host.status)

# % true positives: 90.52%
round(conf_matrix[2,2]/(conf_matrix[1,2] + conf_matrix[2,2])*100, 2)

# % true negatives: 90.64%
round(conf_matrix[1,1]/(conf_matrix[1,1] + conf_matrix[2,1])*100, 2)


# ii) Plot confusion matrix as a fourfold plot
# See https://stackoverflow.com/questions/23891140/r-how-to-visualize-confusion-ma
# trix-using-the-caret-package
#            0 (Obs)   1 (Obs)
# 0 (Pred)   TN        FN
# 1 (pred)   FP        TP

fourfoldplot(conf_matrix, color = c("#CC6666", "#99CC99"),
             conf.level = 0, main = "Confusion Matrix",
             margin = c(2))       # standardizing the col



#
#####      PLOT PREDICTIONS BY OAK/AGRILUS SPECIES                                  #####
# A) PER OAK SPECIES
# i) All oak species
ggplot(predictions,
       aes(x = host.status, y = prediction.logodds)) +
  geom_point() +
  geom_hline(yintercept = cutoff[2], linewidth = 1.5,
             linetype = "dotted", col = "maroon") +
  facet_wrap(~ quercus.sp, ncol = 16) +
  theme(strip.text = element_text(face = "italic"))


# ii) Only hosts
hosts <- unique(subset(predictions, host.status == 1)$quercus.sp)
hosts <- subset(predictions, quercus.sp %in% hosts)

ggplot(hosts, aes(x = host.status, y = prediction.logodds)) +
  geom_point() +
  geom_hline(yintercept = cutoff[2], linewidth = 1.5,
             linetype = "dotted", col = "maroon") +
  facet_wrap(~ quercus.sp, ncol = 10)



# B) PER AGRILUS SPECIES
# i) For all individual Agrilus sp.
ggplot(predictions,
       aes(x = host.status, y = prediction.logodds)) +
  geom_hline(yintercept = cutoff[2], linewidth = 1.5,
             linetype = "dotted", col = "maroon") +
  geom_point() + facet_wrap(~ agrilus.sp, ncol = 8) +
  theme(strip.text = element_text(face = "italic"))


# ii) For a particular Agrilus sp.
# Continuous plot
q <- subset(predictions, agrilus.sp == "A. angustulus" & host.status == 0 &
              prediction.logodds >= cutoff[2])$quercus.sp

ggplot(subset(predictions, agrilus.sp == "A. angustulus"),
       aes(x = host.status, y = prediction.logodds)) + # , fill = host.status))
  # geom_violin(width = 1) +
  # geom_boxplot(width = 0.1, fill = "grey") +
  # scale_fill_manual(values = c("aquamarine2", "coral")) +
  geom_point(size = 5) +
  facet_wrap(~ agrilus.sp, ncol = 8) +
  geom_hline(yintercept = cutoff[2],
             linetype = "dashed", color = "red") +
  ggrepel::geom_label_repel(aes(label = ifelse(quercus.sp %in% q,
                                               as.character(quercus.sp), ''),),
                            hjust = 0, vjust = 0, max.overlaps = 100, size = 4)

# save.image("results/oak_models.RData")



#