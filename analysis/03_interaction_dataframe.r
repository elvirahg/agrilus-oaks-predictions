#### SET ENVIRONMENT ####
# Custom functions
source("R/functions.R")

# Packages
library(dplyr)
library(ggplot2)

# Seed
set.seed(24601)


#### INITIALISE DATAFRAME TO STORE INFO ON AGRILUS - OAK INTERACTIONS ####
# Read in cleaned-up and filtered GBIF plant geo data
plant_geo <- read.table("data/tmp/gbif_plants_clean.tsv",
                        header = TRUE,
                        sep = "\t")

# Extract info on Agrilus species that use oaks, and their hosts
agrilus_hosts <- read.table("data/input/agrilus_quercus_hosts.txt",
                            header = FALSE,
                            sep = "\t")
colnames(agrilus_hosts) <- c("agrilus", "hosts")

agrilus_hosts <- rbind(agrilus_hosts,
                       read.table("data/input/non_oak_hosts.tsv",
                                  header = TRUE,
                                  sep = "\t"))

# Swap order and rename columns
agrilus_hosts <- agrilus_hosts[, c("hosts", "agrilus")]
colnames(agrilus_hosts) <- c("plant.sp", "agrilus.sp")

# Remove hosts not in the phylogeny (Q. gambelii)
agrilus_hosts <- agrilus_hosts[!grepl("gambelii", agrilus_hosts$plant.sp), ]

# Remove non-native hosts for Agrilus species introduced into new areas
agrilus_hosts <- agrilus_hosts |>
  filter(!(plant.sp == "Quercus agrifolia"
           & agrilus.sp == "Agrilus auroguttatus"),
         !(plant.sp == "Quercus macrocarpa"
           & agrilus.sp == "Agrilus sulcicollis"))

# Remove larval hosts with a confidence index < 3 in Jendek & Polakova (2014)
agrilus_hosts <- agrilus_hosts |>
  filter(!(plant.sp == "Ficus carica" & agrilus.sp == "Agrius obscuricollis"),
         !(plant.sp == "Quercus robur" & agrilus.sp == "Agrilus relegatus"),
         !(plant.sp %in% c("Corylus avellana", "Ostrya carpinifolia",
                           "Euonymus europaeus", "Castanea sativa")
           & agrilus.sp == "Agrilus graminis"))


#### ADD OAK GEO INFO TO INTERACTION DATAFRAME ####
# Initialise interaction datafame
oak_spp <- grep("Quercus", unique(plant_geo$plant.sp), value = TRUE)
interaction_data <- create_interaction_df(known_interactions = agrilus_hosts,
                                          host_col = "plant.sp",
                                          hosted_col = "agrilus.sp",
                                          host_taxa_include = oak_spp,
                                          verbose = TRUE)

# For each plant sp. - Agrilus sp. pair, find the mean min. distance to the
# nearest host (from a different plant sp.) of the current Agrilus sp.
# NB: The distance calculated is very different for large dists when using the
# 'cheap' metric (vs., e.g., 'harvestine'). Also note that, with 'cheap', the
# distance values depend on the number of location points used. Still, as we're
# transforming these distances, it's not something to worry about with our data.
interaction_data <- generate_dist_metrics_df(
  expanded_interaction_df = interaction_data,
  coords_df = plant_geo,
  known_interactions_df = agrilus_hosts,
  host_taxon_col = "plant.sp",
  hosted_taxon_col = "agrilus.sp",
  subsample_coords = 10000,
  measure = "cheap",
  metric_type = c("mean", "median", "norm", "log"),
  max_dist_km = 30000,
)

# write.table(x = interaction_data,
#             file = "data/tmp/interaction_data.tsv",
#             row.names = FALSE,
#             col.names = TRUE,
#             sep = "\t",
#             quote = FALSE)


#### EXPLORE RESULTS ####
# interaction_data <- read.table("data/tmp/interaction_data.tsv",
#                                header = TRUE,
#                                sep = "\t")

# Preliminary plots
par(mfrow = c(1, 2))
hist(interaction_data$min.dist.mean[interaction_data$interaction == 1],
     seq(0, 3e04, length.out = 30),
     col = hcl.colors(1, "Warm"),
     main = "Values for hosts - other hosts of a given species",
     xlab = "mean(min distance to any host of the given species)")

hist(interaction_data$min.dist.mean[interaction_data$interaction == 0],
     seq(0, 3e04, length.out = 30),
     col = hcl.colors(1, "Dynamic"),
     main = "Values for non-hosts – hosts of a given species",
     xlab = "mean(min distance to any host of the given species)")
par(mfrow = c(1, 1))

# Plot all thre geo distance metrics to compare them (mean)
plot_geo_dists(interaction_data,
               dist_cols = c("min.dist.mean",
                             "min.dist.mean.norm",
                             "min.dist.mean.log"))
