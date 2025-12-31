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
library(choroplethr)
library(choroplethrMaps)
# library(raster)
# library(geodata)
# library(ggbiplot)
# library(factoextra)
# library(doParallel)



#
#### QUICK EXPLORATION OF INPUT HOST DATA ####
# A) PREPARE INPUT
# Read in file with extended host status information
host_info <- read.csv("data/input/host_info_references.csv")
host_info <- host_info[, c("agrilus_sp",
                           "native_region_insect",
                           "host_family",
                           "host_sp",
                           "record_used",
                           "native_region_host")]

# Remove Q. gambelii as it's not in the final dataset (missing from phylogeny)
host_info <- subset(host_info, host_sp != "Quercus gambelii")

# Add columns of interest
host_info <- host_info |>
  dplyr::mutate(host_genus = stringr::word(host_sp, 1),
                host_is_oak = host_genus == "Quercus",
                host_is_sp = stringr::str_detect(host_sp, "\\s"))


# B) EXPLORE AGRILUS
agrilus_summary <- host_info |>
  dplyr::group_by(agrilus_sp) |>
  dplyr::summarise(only_oaks = all(host_is_oak),
                   number_host_spp = dplyr::n_distinct(host_sp[host_is_sp]),
                   number_host_genera = dplyr::n_distinct(host_genus),
                   number_continents = dplyr::n_distinct(
                     unlist(stringr::str_split(native_region_host, ";\\s*"))
                   ),
                   .groups = "drop")
agrilus_summary <- host_info |>
  dplyr::group_by(agrilus_sp) |>
  dplyr::summarise(only_oaks = all(host_is_oak),
                   number_host_spp = dplyr::n_distinct(host_sp[host_is_sp]),
                   number_host_genera = dplyr::n_distinct(host_genus),
                   number_continents = dplyr::n_distinct(
                     unlist(stringr::str_split(native_region_host, ";\\s*"))
                   ),
                   .groups = "drop")

# % Agrilus with only oak hosts (56.25%, 18 spp)
table(agrilus_summary$only_oaks)
mean(agrilus_summary$only_oaks) * 100

# Number of host spp
hist(table(agrilus_summary$number_host_spp),
     xlab = "Number of host species",
     main = "Number of host species per Agrlius species")

summary(agrilus_summary$number_host_spp)

# Number of host genera
hist(table(agrilus_summary$number_host_genera),
     xlab = "Number of host genera",
     main = "Number of host genera per Agrlius species")
summary(agrilus_summary$number_host_genera)

# Number of continents
sort(table(unlist(stringr::str_split(host_info$native_region_insect, ";\\s*"))))
table(agrilus_summary$number_continents)
summary(agrilus_summary$number_continents)


# C) EXPLORE OAKS
quercus_summary <- host_info |>
  dplyr::filter(host_is_oak, host_is_sp) |>
  dplyr::group_by(host_sp) |>
  dplyr::summarise(number_agrilus_spp = dplyr::n_distinct(agrilus_sp),
                   number_continents = dplyr::n_distinct(
                     unlist(stringr::str_split(native_region_host, ";\\s*"))
                   ),
                   .groups = "drop")

# Number of Agrilus species
hist(table(quercus_summary$number_agrilus_spp),
     xlab = "Number of Agrilus species",
     main = "Number of Agrilus species hosted per oak species")

summary(quercus_summary$number_agrilus_spp)

# Number of continents
sort(table(unlist(stringr::str_split(host_info$native_region_host, ";\\s*"))))
table(agrilus_summary$number_continents)
summary(agrilus_summary$number_continents)



#
#####      EXPLORE MODEL 43                                                         #####
load("analysis/results/oak_models.RData")

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
  mod043 = fitted(oak_mod043, scale = "linear")[,1],  # full
  mod002 = fitted(oak_mod002, scale = "linear")[,1],  # geo
  mod004 = fitted(oak_mod004, scale = "linear")[,1],  # phy1
  mod008 = fitted(oak_mod008, scale = "linear")[,1],  # (1 | gr(oak | phy2))
  mod010 = fitted(oak_mod010, scale = "linear")[,1],  # geo + phy1
  mod018 = fitted(oak_mod018, scale = "linear")[,1],  # geo + (1 | gr(oak | phy2))
  mod025 = fitted(oak_mod025, scale = "linear")[,1]   # phy1 + (1 | gr(oak | phy2))
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

# save.image("analysis/results/oak_models.RData")



#
#####      PLOT DISITRIBUTION OF HOSTS                                              #####
# A) PHYLO
oak_tree <- read.tree("analysis/input/quercus_hipp19_singleton_crown_sp_level.tre")
plants <- gsub(oak_tree$tip.label, pattern = "\\|.*", replacement = "")
plants <- gsub(plants, pattern = "_", replacement = " ", fixed = TRUE)

phosts  <- data.frame(matrix(nrow = 238, ncol = 2))
colnames(phosts) <- c("quercus.sp", "val")
phosts$quercus.sp <- oak_tree$tip.label

vals <- unique(subset(predictions, prediction.cutoff == 1 & host.status == 0)$quercus.sp)
phosts$val <- phosts$quercus.sp %in% vals

x <- full_join(as_tibble(oak_tree), phosts, by = c("label" = "quercus.sp"))
oak_tree_hosts <- treeio::as.treedata(x)

ggtree(oak_tree_hosts, layout = "circular") +
  geom_tiplab(aes(subset = !val), offset = 5, size = 1.7,
              fontface = "italic", color = "gray50") +
  geom_tiplab(aes(subset = val), offset = 5, size = 1.7,
              fontface = "bold.italic", color = "darkgoldenrod3")


# B) GEO
distribution <- read.table(paste0(pw, "input/gbif_quercus_country.tsv"), header = F, sep = "\t")
colnames(distribution) <- c("quercus.sp", "region")
head(distribution)

distribution$region <- countrycode::countrycode(distribution$region, "iso2c", "country.name")

# Attach the country map used in the choropleth function to make sure that country names in
# our dataset are the same
data(country.map)
distribution$region <- tolower(distribution$region)

unique(distribution$region)[which(!(unique(distribution$region) %in% unique(country.map$region)))]

distribution$region[which(distribution$region == "united states")] <- "united states of america"
distribution$region[which(distribution$region == "myanmar (burma)")] <- "myanmar"
distribution$region[which(distribution$region == "bosnia & herzegovina")] <- "bosnia and herzegovina"
distribution$region[which(distribution$region == "czechia")] <- "czech republic"
distribution$region[which(distribution$region == "isle of man")] <- "united kingdom"
distribution$region[which(distribution$region == "north macedonia")] <- "macedonia"
distribution$region[which(distribution$region == "serbia")] <- "republic of serbia"
distribution$region[which(distribution$region == "guernsey")] <- "united kingdom"
distribution$region[which(distribution$region == "jersey")] <- "united kingdom"
distribution$region[which(distribution$region == "tanzania")] <- "united republic of tanzania"
distribution$region[which(distribution$region == "congo - kinshasa")] <- "democratic republic of the congo"

# Create a DS with the number of oak species each country has, and plot them
distribution_oaks <-  as.data.frame(table(distribution$region))
colnames(distribution_oaks) <- c("region", "value")
distribution_oaks <- distribution_oaks[distribution_oaks$region %in% unique(country.map$region),]

extra <- unique(country.map$region)[which(!(unique(country.map$region) %in% distribution_oaks$region))]
extra <- data.frame(cbind(extra, rep(0, length(extra))))
colnames(extra) <- colnames(distribution_oaks)
extra$value <- as.integer(extra$value)

distribution_oaks <- rbind(distribution_oaks, extra); rm(extra)

country_choropleth(distribution_oaks, num_colors = 9) +
  scale_fill_brewer(palette = 10)

# Create a DS with the number of known hosts species each country has, and plot them
distribution_hosts <-  as.data.frame(table(subset(distribution,
                                                  quercus.sp %in% unique(agrilus_hosts$plant.sp))$region))
colnames(distribution_hosts) <- c("region", "value")
distribution_hosts <- distribution_hosts[distribution_hosts$region %in% unique(country.map$region),]

extra <- unique(country.map$region)[which(!(unique(country.map$region) %in% distribution_hosts$region))]
extra <- data.frame(cbind(extra, rep(0, length(extra))))
colnames(extra) <- colnames(distribution_hosts)
extra$value <- as.integer(extra$value)

distribution_hosts <- rbind(distribution_hosts, extra); rm(extra)

country_choropleth(distribution_hosts, num_colors = 9) +
  scale_fill_brewer(palette = 10)

# Create a DS with the number of predicted hosts species each country has, and plot them
distribution_phosts <-  as.data.frame(
  table(subset(distribution, quercus.sp %in%
                 gsub("Q\\.", "Quercus",
                      unique(subset(predictions, prediction.cutoff == 1)$quercus.sp)))$region))
colnames(distribution_phosts) <- c("region", "value")
distribution_phosts <- distribution_phosts[distribution_phosts$region %in% unique(country.map$region),]

extra <- unique(country.map$region)[which(!(unique(country.map$region) %in% distribution_phosts$region))]
extra <- data.frame(cbind(extra, rep(0, length(extra))))
colnames(extra) <- colnames(distribution_phosts)
extra$value <- as.integer(extra$value)

distribution_phosts <- rbind(distribution_phosts, extra); rm(extra)

country_choropleth(distribution_phosts, num_colors = 9) +
  scale_fill_brewer(palette = 10)



#
#####      EXPLORE AGRILUS SPECIES WHOSE REAL HOSTS ARE NOT PREDICTED AS SUCH       #####
# It seems to happen when:
# The beetle is only hosted by one oak species
# The hosts are far in the phylogeny ("cold clades") or geo
subset(predictions, prediction.cutoff == 0 & host.status == 1)

# A) SPECIES WITH ONLY ONE KNOWN HOST
# i) A. samai
unique(subset(predictions, agrilus.sp == "A. samai" & host.status == 1)$quercus.sp)

# ii) A. niveoguttatus
unique(subset(predictions, agrilus.sp == "A. niveoguttatus" & host.status == 1)$quercus.sp)

# iii) A. chiricahuae
unique(subset(predictions, agrilus.sp == "A. chiricahuae" & host.status == 1)$quercus.sp)

# iv) A. relegatoides
unique(subset(predictions, agrilus.sp == "A. relegatoides" & host.status == 1)$quercus.sp)

# v) A. acutipennis
unique(subset(predictions, agrilus.sp == "A. acutipennis" & host.status == 1)$quercus.sp)


# B) SPECIES WITH MORE THAN ONE KNOWN HOST
# i) A. coxalis (0s are far in geo and phylo)
# Check hosts
subset(predictions, agrilus.sp == "A. coxalis" & host.status == 1)

# Plot phylogeny
hosts  <- data.frame(matrix(nrow = 238, ncol = 2))
colnames(hosts) <- c("quercus.sp", "val")
hosts$quercus.sp <- oak_tree$tip.label

vals <- unique(subset(predictions, agrilus.sp == "A. coxalis" & host.status == 1)$quercus.sp)
hosts$val <- hosts$quercus.sp %in% vals

x <- full_join(as_tibble(oak_tree), hosts, by = c("label" = "quercus.sp"))
oak_tree2 <- treeio::as.treedata(x)
rm(x, hosts)

ggtree(oak_tree2, layout = "circular") +
  geom_tiplab(aes(color = val), offset = 5, size = 1.7, fontface = "italic") +
  scale_color_manual(values = c("black", "coral"))

# Check phylo distance
pdists <- subset(interaction_data, agrilus.sp == "A. coxalis" &
                       quercus.sp %in% c("Q. conzattii", "Q. peduncularis",
                                         "Q. chrysolepis", "Q. kelloggii",
                                         "Q. agrifolia"))[,c(1,3,12)]

barplot(pdists$phylo.dist.min, names.arg = pdists$quercus.sp,
        ylab = "min phylo mean",
        las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(rep(1, 6))])

# Check ranef values (high values make it more likely for sth to be a host)
ranef <- data.frame(rbind(ranef(oak_mod043)$quercus.sp[,,1]['Q. conzattii',],
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. peduncularis',],
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. chrysolepis',],
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. kelloggii',],
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. agrifolia',]
                          ))
ranef$quercus.sp <- c("Q. conzattii", "Q. peduncularis", "Q. chrysolepis", "Q. kelloggii",
                      "Q. agrifolia")

barplot(ranef$Estimate, names.arg = ranef$quercus.sp,
        ylab = "ranef estimate", las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(1,1,1,1,1,2)])

# Geographic distances
subset(interaction_data, agrilus.sp == "A. coxalis" & interaction == 1)[, c(1,2,4,5)]

gdists <- subset(interaction_data, agrilus.sp == "A. coxalis" &
                                       quercus.sp %in% c("Q. conzattii", "Q. peduncularis",
                                                         "Q. chrysolepis", "Q. kelloggii",
                                                         "Q. agrifolia"))[,c(1,3,5)]

barplot(gdists$min.dist.mean.norm, names.arg = gdists$quercus.sp,
        ylab = "min.dist.mean (geo)",
        las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(rep(1, 5), 2)])


# ii) A. hemiphanes (0s are far in geo and phylo)
# Check hosts
subset(predictions, agrilus.sp == "A. hemiphanes" & host.status == 1)

# Plot phylogeny
hosts  <- data.frame(matrix(nrow = 238, ncol = 2))
colnames(hosts) <- c("quercus.sp", "val")
hosts$quercus.sp <- oak_tree$tip.label

vals <- unique(subset(predictions, agrilus.sp == "A. hemiphanes" & host.status == 1)$quercus.sp)
hosts$val <- hosts$quercus.sp %in% vals

x <- full_join(as_tibble(oak_tree), hosts, by = c("label" = "quercus.sp"))
oak_tree2 <- treeio::as.treedata(x)
rm(x, hosts)

ggtree(oak_tree2, layout = "circular") +
  geom_tiplab(aes(color = val), offset = 5, size = 1.7, fontface = "italic") +
  scale_color_manual(values = c("black", "coral"))

# Check phylo distance
pdists <- subset(interaction_data, agrilus.sp == "A. hemiphanes" &
                       quercus.sp %in% c("Q. coccinea", "Q. coccifera",
                                         "Q. ilex"))[,c(1,3,12)]

barplot(pdists$phylo.dist.min, names.arg = pdists$quercus.sp,
        ylab = "min phylo mean",
        las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(rep(1, 6))])

# Check ranef values (high values make it more likely for sth to be a host)
ranef <- data.frame(rbind(ranef(oak_mod043)$quercus.sp[,,1]['Q. coccinea',],
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. coccifera',],
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. ilex',]
                          ))
ranef$quercus.sp <- c("Q. coccinea", "Q. coccifera",  "Q. ilex")

barplot(ranef$Estimate, names.arg = ranef$quercus.sp,
        ylab = "ranef estimate", las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(1,1,1,1,1,2)])

# Geographic distances
subset(interaction_data, agrilus.sp == "A. hemiphanes" & interaction == 1)[, c(1,2,4,5)]

gdists <- subset(interaction_data, agrilus.sp == "A. hemiphanes" &
                                       quercus.sp %in% c("Q. coccinea", "Q. coccifera",
                                                         "Q. ilex"))[,c(1,3,5)]

barplot(gdists$min.dist.mean.norm, names.arg = gdists$quercus.sp,
        ylab = "min.dist.mean (geo)",
        las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(rep(1, 5), 2)])


# iii) A. defectus (far in phylo)
subset(predictions, agrilus.sp == "A. defectus" & host.status == 1)

# Plot phylogeny
hosts  <- data.frame(matrix(nrow = 238, ncol = 2))
colnames(hosts) <- c("quercus.sp", "val")
hosts$quercus.sp <- oak_tree$tip.label

vals <- unique(subset(predictions, agrilus.sp == "A. defectus" & host.status == 1)$quercus.sp)
hosts$val <- hosts$quercus.sp %in% vals

x <- full_join(as_tibble(oak_tree), hosts, by = c("label" = "quercus.sp"))
oak_tree2 <- treeio::as.treedata(x)
rm(x, hosts)

ggtree(oak_tree2, layout = "circular") +
  geom_tiplab(aes(color = val), offset = 5, size = 1.7, fontface = "italic") +
  scale_color_manual(values = c("black", "coral"))

# Check phylo distance
pdists <- subset(interaction_data, agrilus.sp == "A. defectus" &
                       quercus.sp %in% c("Q. muehlenbergii", "Q. stellata",
                                         "Q. alba"))[,c(1,3,12)]

barplot(pdists$phylo.dist.min, names.arg = pdists$quercus.sp,
        ylab = "min phylo mean",
        las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(rep(1, 6))])

# Check ranef values (high values make it more likely for sth to be a host)
ranef <- data.frame(rbind(ranef(oak_mod043)$quercus.sp[,,1]['Q. muehlenbergii',],
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. stellata',],
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. alba',]
                          ))
ranef$quercus.sp <- c("Q. muehlenbergii", "Q. stellata", "Q. alba")

barplot(ranef$Estimate, names.arg = ranef$quercus.sp,
        ylab = "ranef estimate", las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(1,1,1,1,1,2)])

# Geographic distances
subset(interaction_data, agrilus.sp == "A. defectus" & interaction == 1)[, c(1,2,4,5)]

gdists <- subset(interaction_data, agrilus.sp == "A. defectus" &
                                       quercus.sp %in% c(c("Q. muehlenbergii", "Q. stellata",
                                                           "Q. alba")))[,c(1,3,5)]

barplot(gdists$min.dist.mean.norm, names.arg = gdists$quercus.sp,
        ylab = "min.dist.mean (geo)",
        las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(rep(1, 5), 2)])


# iv) A. albocomus (far in phylo, somewhat geo)
subset(predictions, agrilus.sp == "A. albocomus" & host.status == 1)

# Plot phylogeny
hosts  <- data.frame(matrix(nrow = 238, ncol = 2))
colnames(hosts) <- c("quercus.sp", "val")
hosts$quercus.sp <- oak_tree$tip.label

vals <- unique(subset(predictions, agrilus.sp == "A. albocomus" & host.status == 1)$quercus.sp)
hosts$val <- hosts$quercus.sp %in% vals

x <- full_join(as_tibble(oak_tree), hosts, by = c("label" = "quercus.sp"))
oak_tree2 <- treeio::as.treedata(x)
rm(x, hosts)

ggtree(oak_tree2, layout = "circular") +
  geom_tiplab(aes(color = val), offset = 5, size = 1.7, fontface = "italic") +
  scale_color_manual(values = c("black", "coral"))

# Check phylo distance
pdists <- subset(interaction_data, agrilus.sp == "A. albocomus" &
                       quercus.sp %in% c("Q. emoryi", "Q. grisea"))[,c(1,3,12)]

barplot(pdists$phylo.dist.min, names.arg = pdists$quercus.sp,
        ylab = "min phylo mean",
        las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(rep(1, 6))])

# Check ranef values (high values make it more likely for sth to be a host)
ranef <- data.frame(rbind(ranef(oak_mod043)$quercus.sp[,,1]['Q. emoryi',],
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. grisea',]
                          ))
ranef$quercus.sp <- c("Q. emoryi", "Q. grisea")

barplot(ranef$Estimate, names.arg = ranef$quercus.sp,
        ylab = "ranef estimate", las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(1,1,1,1,1,2)])

# Geographic distances
subset(interaction_data, agrilus.sp == "A. albocomus" & interaction == 1)[, c(1,2,4,5)]

gdists <- subset(interaction_data, agrilus.sp == "A. albocomus" &
                                       quercus.sp %in% c("Q. emoryi", "Q. grisea"))[,c(1,3,5)]

barplot(gdists$min.dist.mean.norm, names.arg = gdists$quercus.sp,
        ylab = "min.dist.mean (geo)",
        las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(rep(1, 5), 2)])



#