#### SET ENVIRONMENT ####
# Custom functions
source("R/phylo_functions.r")
source("R/geo_functions.r")
source("R/analysis_functions.r")
source("R/plot_functions.r")

# Libraries
library(brms)
library(dplyr)
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

# Add binary predictions to predictions dataframe
predictions$prediction_binary <- factor(
  ifelse(predictions$prediction_lodds > thr_intercepts["threshold"], 1, 0),
  levels = c(0, 1)
)

# Write predictions to file
# write.table(predictions,
#             "data/results/predictions.tsv",
#             row.names = FALSE,
#             col.names = TRUE,
#             sep = "\t",
#             quote = FALSE)


#### PLOT PREDICTION RESULTS (GRAPHS) ####
# Genral violin plot
plot_predictions_violin(predictions = predictions,
                        threshold = thr_intercepts["threshold"])

median(predictions[predictions$host_status == 1, ]$prediction_lodds)
median(predictions[predictions$host_status == 0, ]$prediction_lodds)

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


#### PLOT PREDICTION RESULTS (PHYLO TREE) ####
# Read in tree
# See Hipp et al., 2020; https://github.com/andrew-hipp/global-oaks-2019
oak_phylo <- ape::read.tree("data/input/tr.singletons.correlated.1.taxaGrepCrown_accepted_names.tre")

# Clean tree
old_names <- c("Quercus litoralis",
               "Quercus new",
               "Quercus sp")
new_names <- c("Quercus litoralis (Atuna excelsa)",
               "Quercus sp. nov. QUE000227",
               "Quercus sp. nov. QUE001568")

oak_phylo <- standardise_phylo(oak_phylo,
                               old_labels = old_names,
                               new_labels = new_names,
                               remove_duplicates = TRUE,
                               clean_labels = TRUE,
                               pattern = "^([A-Z][a-z]+)_([×|x]*)_*([a-z-]+).*",
                               replacement = "\\1 \\2\\3")

# Identify known and predicted host species
pred_hosts <- unique(subset(predictions,
                            prediction_binary == 1)$quercus_sp)
known_hosts <- unique(subset(predictions,
                             host_status == 1)$quercus_sp)

predicted_hosts <- data.frame(
  quercus_sp = oak_phylo$tip.label,
  is_known = oak_phylo$tip.label %in% known_hosts,
  is_pred = oak_phylo$tip.label %in% pred_hosts
)

# Plot oak phylogeny with host info
plot_phylo_hosts(phylo = oak_phylo,
                 pred_hosts = predicted_hosts,
                 label = "quercus_sp",
                 layout = "circular",
                 tiplab_offset = 5,
                 tiplab_size = 1.7,
                 tippoint_shape = 18,
                 tippoint_size = 2,
                 pred_color = "darkgoldenrod3",
                 nonpred_color = "gray50",
                 known_color = "#cf3530",
                 show_legend = FALSE)


#### PLOT PREDICTION RESULTS (GEO CHOROPLETHS) ####
# Read WCVP data
wcvp_v10_names <- read.table("data/input/wcvp_v10/wcvp_names.csv",
                             sep = "|",
                             header = TRUE,
                             quote = "",
                             comment.char = "")
wcvp_v10_distributions <- read.table("data/input/wcvp_v10/wcvp_distribution.csv",
                                     sep = "|",
                                     header = TRUE,
                                     quote = "",
                                     comment.char = "")

# Oak species in dataset
plant_names <- unique(predictions$quercus_sp)

# Replace names in dataset that need to be looked up as someting else
lookup <- c("Quercus frainetto" = "Quercus conferta",
            "Quercus litoralis (Atuna excelsa)" = "Atuna excelsa",
            "Quercus margarettae" = "Quercus margaretiae",
            "Quercus ×crenata" = "Quercus × crenata")

taxon_queries <- ifelse(plant_names %in% names(lookup),
                        lookup[plant_names], plant_names)

# Retrieve WCVP distribution information for each oak species in dataset
oak_distributions_list <- wcvp_distribution_list(
  taxon_queries = taxon_queries,
  plant_names = plant_names,
  known_hosts = known_hosts,
  pred_hosts = pred_hosts,
  wcvp_names = wcvp_v10_names,
  wcvp_distributions = wcvp_v10_distributions
)

# Combine into a single data frame
oak_distributions <- do.call(rbind, oak_distributions_list)

# Prepare palettes
palette_vals_l3 <- colorRampPalette(RColorBrewer::brewer.pal(9, "PuBuGn"))(12)
names(palette_vals_l3) <- paste0(seq(1, 56, 5), "–", seq(5, 60, 5))

palette_vals_country <- c("white",
                          colorRampPalette(RColorBrewer::brewer.pal(9, "PuBuGn"))(13))
names(palette_vals_country) <- c("0",
                                 paste0(seq(1, 130, 10), "–", seq(10, 130, 10)))

# Plot total oak species (LEVEL3_NAM)
plot_distribution_region(df = oak_distributions,
                         group_level = "LEVEL3_NAM",
                         counts_type = "species",
                         breaks = c(-Inf, seq(0, 60, 5)),
                         labels = c("0", names(palette_vals_l3)),
                         world_sf = NULL,
                         palette_vals = palette_vals_l3,
                         na_col = "gray90")

plot_distribution_region(df = oak_distributions,
                         group_level = "LEVEL3_NAM",
                         counts_type = "species_native",
                         breaks = c(-Inf, seq(0, 60, 5)),
                         labels = c("0", names(palette_vals_l3)),
                         world_sf = NULL,
                         palette_vals = palette_vals_l3,
                         na_col = "gray90")

# Plot known hosts (LEVEL3_NAM)
plot_distribution_region(df = oak_distributions,
                         group_level = "LEVEL3_NAM",
                         counts_type = "hosts",
                         breaks = c(-Inf, seq(0, 60, 5)),
                         labels = c("0", names(palette_vals_l3)),
                         world_sf = NULL,
                         palette_vals = palette_vals_l3,
                         na_col = "gray90")

plot_distribution_region(oak_distributions,
                         group_level = "LEVEL3_NAM",
                         counts_type = "hosts_native",
                         breaks = c(-Inf, seq(0, 60, 5)),
                         labels = c("0", names(palette_vals_l3)),
                         world_sf = NULL,
                         palette_vals = palette_vals_l3,
                         na_col = "gray90")

# Plot predicted hosts (LEVEL3_NAM)
plot_distribution_region(df = oak_distributions,
                         group_level = "LEVEL3_NAM",
                         counts_type = "pred_hosts",
                         breaks = c(-Inf, seq(0, 60, 5)),
                         labels = c("0", names(palette_vals_l3)),
                         world_sf = NULL,
                         palette_vals = palette_vals_l3,
                         na_col = "gray90")

plot_distribution_region(df = oak_distributions,
                         group_level = "LEVEL3_NAM",
                         counts_type = "pred_hosts_native",
                         breaks = c(-Inf, seq(0, 60, 5)),
                         labels = c("0", names(palette_vals_l3)),
                         world_sf = NULL,
                         palette_vals = palette_vals_l3,
                         na_col = "gray90")

# Plot total oak species (country)
plot_distribution_country(df = oak_distributions,
                          counts_type = "species",
                          breaks = c(0, 1, seq(10, 130, 10)),
                          labels = names(palette_vals_country),
                          palette_vals = palette_vals_country,
                          na_col = "gray90")
plot_distribution_country(df = oak_distributions,
                          counts_type = "species_native",
                          breaks = c(0, 1, seq(10, 130, 10)),
                          labels = names(palette_vals_country),
                          palette_vals = palette_vals_country,
                          na_col = "gray90")

# Plot known hosts (country)
plot_distribution_country(df = oak_distributions,
                          counts_type = "hosts",
                          breaks = c(0, 1, seq(10, 130, 10)),
                          labels = names(palette_vals_country),
                          palette_vals = palette_vals_country,
                          na_col = "gray90")
plot_distribution_country(df = oak_distributions,
                          counts_type = "hosts_native",
                          breaks = c(0, 1, seq(10, 130, 10)),
                          labels = names(palette_vals_country),
                          palette_vals = palette_vals_country,
                          na_col = "gray90")

# Plot predicted hosts (country)
plot_distribution_country(df = oak_distributions,
                          counts_type = "pred_hosts",
                          breaks = c(0, 1, seq(10, 130, 10)),
                          labels = names(palette_vals_country),
                          palette_vals = palette_vals_country,
                          na_col = "gray90")
plot_distribution_country(df = oak_distributions,
                          counts_type = "pred_hosts_native",
                          breaks = c(0, 1, seq(10, 130, 10)),
                          labels = names(palette_vals_country),
                          palette_vals = palette_vals_country,
                          na_col = "gray90")


##### EXPLORE AGRILUS SPECIES WHOSE REAL HOSTS ARE NOT PREDICTED AS SUCH #####
# Seems to happen when:
# Beetle only has one host species
# Hosts are far in phylogeny ("cold clades") or geography
subset(predictions, prediction_binary == 0 & host_status == 1)

# i) Explore species with only one known host
unique(subset(interaction_data,
              agrilus_sp == "Agrilus samai"
              & interaction == 1)$quercus_sp)
unique(subset(interaction_data,
              agrilus_sp == "Agrilus niveoguttatus"
              & interaction == 1)$quercus_sp)
unique(subset(interaction_data,
              agrilus_sp == "Agrilus chiricahuae"
              & interaction == 1)$quercus_sp)
unique(subset(interaction_data,
              agrilus_sp == "Agrilus relegatoides"
              & interaction == 1)$quercus_sp)
unique(subset(interaction_data,
              agrilus_sp == "Agrilus acutipennis"
              & interaction == 1)$quercus_sp)

# ii) Explore all other species
# Colour legend:
# * Blue: oak species is a known host of the given beetle in its current range
# * Light blue: predicted value fell under the binary threshold.
# For Agrilus auroguttatus:
# * Orange: novel host in its invasive range.
# * Lighter orange: predicted value fell under the binary threshold.
# * Pink: known non-host in the beetle’s introduced range.

# ii.i) Agrilus auroguttatus
a_auroguttatus_map <- c(
  # known hosts in native range
  "Quercus conzattii" = "#e6c27e",
  "Quercus hypoleucoides" = "#ae6773",
  "Quercus peduncularis" = "#e6c27e",
  "Quercus emoryi" = "#2775b4",
  # known novel hosts in new range
  "Quercus engelmannii" = "#2775b4",
  "Quercus chrysolepis" = "#e6c27e",
  "Quercus kelloggii" = "#2775b4",
  "Quercus agrifolia" = "#f0a720",
  # known non-host
  "Quercus arizonica" = "#7c9db9"
)

# Plot known hosts phylogeny
plot_phylogeny(phylo = oak_phylo,
               interaction_df = interaction_data,
               species = "Agrilus auroguttatus",
               host_spp = names(a_auroguttatus_map),
               main = "Agrilus auroguttatus",
               layout = "circular",
               offset = 5,
               size = 1.7,
               tip_colours = c("black", "coral"),
               legend_pos = "none")

# Plot phylo_dist_min values for known hosts
barplot_host_distance(interaction_df = interaction_data,
                      species = "Agrilus auroguttatus",
                      species_col = "agrilus_sp",
                      host_col = "quercus_sp",
                      host_status_col = "interaction",
                      dist_col = "phylo_dist_min",
                      host_spp = names(a_auroguttatus_map),
                      colour = a_auroguttatus_map,
                      cex_names = 0.5,
                      ylab = paste("Minimum phylogenetic distance to",
                                   "another host species"))

# Plot geo_min_dist_mean_norm values for known hosts
barplot_host_distance(interaction_df = interaction_data,
                      species = "Agrilus auroguttatus",
                      species_col = "agrilus_sp",
                      host_col = "quercus_sp",
                      host_status_col = "interaction",
                      dist_col = "geo_min_dist_mean_norm",
                      host_spp = names(a_auroguttatus_map),
                      colour = a_auroguttatus_map,
                      cex_names = 0.5,
                      ylab = paste("Minimum normalised geographic distance",
                                   "to another host species"))

# Plot random effects for known hosts
barplot_ranef(model = oak_models$mod046,
              host_spp = names(a_auroguttatus_map),
              species = "Agrilus auroguttatus",
              species_col = "agrilus_sp",
              host_col = "quercus_sp",
              host_status_col = "interaction",
              colour = a_auroguttatus_map)

# ii.ii) Rest of the Agrilus species
agrilus_list <- list(
  "Agrilus albocomus" = c("#2775b4", "#7c9db9"),
  "Agrilus coxalis" = c("#2775b4", "#2775b4", "#7c9db9", "#2775b4", "#7c9db9"),
  "Agrilus defectus" = c("#2775b4", "#2775b4", "#7c9db9"),
  "Agrilus hemiphanes" = c("#7c9db9", "#2775b4", "#2775b4")
)

for (i in seq_along(agrilus_list)) {
  # Plot known hosts phylogeny
  plot_phylogeny(phylo = oak_phylo,
                 interaction_df = interaction_data,
                 species = names(agrilus_list)[i],
                 species_col = "agrilus_sp",
                 host_status_col = "interaction",
                 host_col = "quercus_sp",
                 main = names(agrilus_list)[i],
                 layout = "circular",
                 offset = 5,
                 size = 1.7,
                 tip_colours = c("black", "coral"),
                 legend_pos = "none")

  # Plot phylo_dist_min values for known hosts
  barplot_host_distance(interaction_df = interaction_data,
                        species = names(agrilus_list)[i],
                        species_col = "agrilus_sp",
                        host_col = "quercus_sp",
                        host_status_col = "interaction",
                        dist_col = "phylo_dist_min",
                        colour = agrilus_list[[i]],
                        cex_names = 0.5,
                        ylab = paste("Minimum phylogenetic distance to",
                                     "another host species"))

  # Plot geo_min_dist_mean_norm values for known hosts
  barplot_host_distance(interaction_df = interaction_data,
                        species = names(agrilus_list)[i],
                        species_col = "agrilus_sp",
                        host_col = "quercus_sp",
                        host_status_col = "interaction",
                        dist_col = "geo_min_dist_mean_norm",
                        colour = agrilus_list[[i]],
                        cex_names = 0.5,
                        ylab = paste("Minimum normalised geographic distance",
                                     "to another host species"))

  # Plot random effects for known hosts
  barplot_ranef(interaction_df = interaction_data,
                model = oak_models$mod046,
                species = names(agrilus_list)[i],
                species_col = "agrilus_sp",
                host_col = "quercus_sp",
                host_status_col = "interaction",
                colour = agrilus_list[[i]])
}
