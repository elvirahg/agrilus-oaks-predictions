#### SET ENVIRONMENT ####
# Custom functions
source("R/phylo_functions.r")
source("R/distance_functions.r")
source("R/plot_functions.r")

# Packages
library(dplyr)
library(ggplot2)

# Seed (not needed for actual analyses, only if exact reproducibility of
# results is desired)s
set.seed(24601)


#### INITIALISE DATAFRAME TO STORE INFO ON AGRILUS - OAK INTERACTIONS ####
# Read in cleaned-up and filtered GBIF plant geo data
plant_geo <- read.table("data/results/gbif_plants_clean.tsv",
                        header = TRUE,
                        sep = "\t")

# Extract info on Agrilus species that use oaks, and their hosts
agrilus_hosts <- read.table("data/input/quercus_hosts.tsv",
                            header = TRUE,
                            sep = "\t")
colnames(agrilus_hosts) <- c("plant_sp", "agrilus_sp")

agrilus_hosts <- rbind(agrilus_hosts,
                       read.table("data/input/non_quercus_hosts.tsv",
                                  header = TRUE,
                                  sep = "\t"))

# Remove hosts not in the phylogeny (Q. gambelii)
agrilus_hosts <- agrilus_hosts[agrilus_hosts$plant_sp != "Quercus gambelii", ]

# Remove non-native hosts for Agrilus species introduced into new areas
agrilus_hosts <- agrilus_hosts |>
  filter(!(plant_sp == "Quercus agrifolia"
           & agrilus_sp == "Agrilus auroguttatus"),
         !(plant_sp == "Quercus macrocarpa"
           & agrilus_sp == "Agrilus sulcicollis"))

# Remove larval hosts with a confidence index < 3 in Jendek & Polakova (2014)
agrilus_hosts <- agrilus_hosts |>
  filter(!(plant_sp == "Ficus carica" & agrilus_sp == "Agrius obscuricollis"),
         !(plant_sp == "Quercus robur" & agrilus_sp == "Agrilus relegatus"),
         !(plant_sp %in% c("Corylus avellana", "Ostrya carpinifolia",
                           "Euonymus europaeus", "Castanea sativa")
           & agrilus_sp == "Agrilus graminis"))


#### ADD OAK GEO INFO TO INTERACTION DATAFRAME ####
# Initialise interaction datafame
oak_spp <- grep("Quercus", unique(plant_geo$plant_sp), value = TRUE)
interaction_data <- create_interaction_df(known_interactions = agrilus_hosts,
                                          host_col = "plant_sp",
                                          hosted_col = "agrilus_sp",
                                          host_taxa_include = oak_spp,
                                          verbose = TRUE)

# For each plant sp. - Agrilus sp. pair, find the mean min. distance to the
# nearest host (from a different plant sp.) of the current Agrilus sp.
# NB: The distance calculated is very different for large dists when using the
# 'cheap' metric (vs., e.g., 'harvestine'). Also note that, with 'cheap', the
# distance values depend on the number of location points used. Still, as we're
# transforming these distances, it's not something to worry about with our data.
interaction_data <- generate_geo_metrics_df(
  expanded_interaction_df = interaction_data,
  coords_df = plant_geo,
  known_interactions_df = agrilus_hosts,
  host_taxon_col = "plant_sp",
  hosted_taxon_col = "agrilus_sp",
  subsample_coords = 10000,
  measure = "cheap",
  metric_type = c("mean", "median", "norm", "log"),
  max_dist_km = 30000,
)

# Explore results using preliminary plots
par(mfrow = c(1, 2))
hist(interaction_data$geo_min_dist_mean[interaction_data$interaction == 1],
     seq(0, 3e04, length.out = 30),
     col = hcl.colors(1, "Warm"),
     main = "Values for hosts - other hosts of a given species",
     xlab = "mean(min distance to any host of the given species)")

hist(interaction_data$geo_min_dist_mean[interaction_data$interaction == 0],
     seq(0, 3e04, length.out = 30),
     col = hcl.colors(1, "Dynamic"),
     main = "Values for non-hosts – hosts of a given species",
     xlab = "mean(min distance to any host of the given species)")
par(mfrow = c(1, 1))

# Plot all thre geo distance metrics to compare them (mean)
plot_geo_dists(interaction_data,
               dist_cols = c("geo_min_dist_mean",
                             "geo_min_dist_mean_norm",
                             "geo_min_dist_mean_log"))


#### ADD OAK OCCURRENCE INFO TO INTERACTION DATAFRAME ####
# Extract no. oak occurrences
oak_occurrences <- data.frame(table(plant_geo$plant_sp))
colnames(oak_occurrences) <- c("plant_sp", "gbif_entries")
oak_occurrences$plant_sp <- as.character(oak_occurrences$plant_sp)

# Append to interaction_data
interaction_data <- merge(interaction_data,
                          oak_occurrences[, c("plant_sp", "gbif_entries")],
                          by = "plant_sp",
                          all.x = TRUE)


#### ADD OAK PHYLO INFO TO INTERACTION DATAFRAME ####
# Read and clean oak phylogenetic tree
oak_phylo <- ape::read.tree("data/input/tr.singletons.correlated.1.taxaGrepCrown_accepted_names.tre")

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

# Generate phylogenetic metrics
interaction_data <- generate_phylo_metrics_df(
  expanded_interaction_df = interaction_data,
  known_interactions_df = agrilus_hosts,
  tree = oak_phylo,
  host_taxon_col = "plant_sp",
  hosted_taxon_col = "agrilus_sp"
)

# Compare metrics
plot(phylo_dist_mean ~ phylo_dist_min, data = interaction_data)

# Rename plant species column, and save tsv
names(interaction_data)[names(interaction_data) == "plant_sp"] <- "quercus_sp"

# write.table(x = interaction_data,
#             file = "data/results/interaction_data.tsv",
#             row.names = FALSE,
#             col.names = TRUE,
#             sep = "\t",
#             quote = FALSE)
