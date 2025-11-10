#### SET ENVIRONMENT ####
# Custom functions
source("R/functions.R")

# Packages
library(ggplot2)
library(ggtree)
library(dplyr)


#### COMPUTE FRITZ AND PURVIS' D SINAL ####
## READ IN AND CLEAN PHYLOGENETIC TREE
# See Hipp et al., 2020; https://github.com/andrew-hipp/global-oaks-2019
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

## READ IN HOST DATA AND GENERATE PRESENCE-ABSENCE HOST STATUS DATAFRAME
oak_host_observations <- read.table("data/input/quercus_hosts_number_agrilus_hosted.txt",
                                    sep = "\t")
colnames(oak_host_observations) <- c("quercus_sp", "no_agrilus_spp")

oak_hosts_df <- generate_pres_abs_df(oak_host_observations$quercus_sp,
                                     oak_phylo$tip.label)

## COMPUTE SINGAL
# Prepare comparative object
oak_phylo$node.label <- NULL
oak_comparative <- caper::comparative.data(phy = oak_phylo,
                                           data = oak_hosts_df,
                                           # Name of the column with spp. names
                                           names.col = species,
                                           # Store variance covariance matrix
                                           vcv = TRUE,
                                           # Don't remove spp. w/o data for some variables
                                           na.omit = FALSE,
                                           # Shows dropped spp.
                                           warn.dropped = TRUE)

# Compute strength of the phylogenetic signal (Fritz & Purvis' D metric)
d_signal_oaks <- caper::phylo.d(data = oak_comparative,
                                binvar = host.status, permut = 1000)

# D = 0.55496, which suggests clumping is moderately strong; phylogenetic
# pattern differs significantly from Brownian & random expectations
d_signal_oaks

# Plot
plot(d_signal_oaks, lwd = 5, cex.axis = 2.5, cex.lab = 2)
axis(side = 2, lwd = 3, lab = FALSE)
axis(side = 1, lwd = 3, lab = FALSE)
box(lwd = 3)


#### EXTRACT GBIF OAK AND HOST DATA ####
## EXTRACT GBIF KEYS FOR OAK SPECIES IN THE PHYLOGENY
# Get GBIF keys for oak species in phylogeny
oak_gbif_keys <- get_gbif_keys(oak_phylo$tip.label, higher_taxon = "220")

# Add keys manually for spp. with missing keys
# Species names: oak_gbif_keys$species[which(is.na(oak_gbif_keys[, 2]))]
# Note that:
# * Q. frainetto: synonym according to GBIF (but accepted by WCVP)
# * Q. faginea, Q. petraea, Q. pyrenaica, Q. imbricaria, Q. rubra, Q. palustris,
# and Q. acuta: function above recognises another spp. first instead
# * Q. ×crenata: name does not match
# * Q. litoralis: Atuna excelsa subsp. excelsa. in GBIF
manual_oak_keys <- c(
  "Quercus sp. nov. QUE000227"        = NA,
  "Quercus faginea"                   = 2881480,
  "Quercus frainetto"                 = 2879504,
  "Quercus petraea"                   = 2880130,
  "Quercus pyrenaica"                 = 2878826,
  "Quercus imbricaria"                = 2881106,
  "Quercus rubra"                     = 2880539,
  "Quercus palustris"                 = 8313153,
  "Quercus sp. nov. QUE001568"        = NA,
  "Quercus ×crenata"                  = 2880003,
  "Quercus litoralis (Atuna excelsa)" = 10839442,
  "Quercus acuta"                     = 2876863
)

oak_gbif_keys$key <- ifelse(
  is.na(oak_gbif_keys$key) & oak_gbif_keys$species %in% names(manual_oak_keys),
  manual_oak_keys[oak_gbif_keys$species],
  oak_gbif_keys$key
)

## EXTRACT NUMBER OF GBIF ENTRIES
# Extract number of GBIF entries (takes a couple minutes)
oak_gbif_counts <- vapply(
  oak_gbif_keys$key,
  function(key) {
    rgbif::occ_count(
      taxonKey = key,
      hasCoordinate = TRUE,
      year = "1950, 2022"
    )
  },
  FUN.VALUE = numeric(1)
)

oak_gbif_entries <- data.frame(species = oak_gbif_keys$species,
                               no.entries = oak_gbif_counts)

nov_spp <- c("Quercus sp. nov. QUE000227", "Quercus sp. nov. QUE001568")
oak_gbif_entries$no.entries[oak_gbif_entries$species %in% nov_spp] <- 0


#### PLOT HOT CLADES ####
## PLOT HOT CLADES
# Load nexus tree created with phylocom nodesig (first tree: host status info)
oak_nodesig <- ape::read.nexus("data/input/quercus_nodesig_result.nex")[[1]]
oak_nodesig$tip.label <- gsub(pattern = "_",
                              replacement = " ",
                              oak_nodesig$tip.label)

# Check hot and cold clades using Mesquite Software, then import resulting tree
oak_nodesig$node.label <- gsub(pattern = "'",
                               replacement = "",
                               oak_nodesig$node.label)

# Group and colour by hot/cold nodes
nodes <- c("N246", "N247",
           "N293", "N294", "N295", "N296", "N297", "N298",
           "N316", "N317",
           "N342", "N343", "N344",
           "N381",
           "N404",
           "N419",
           "N431",
           "N443", "N444", "N445", "N446")

oak_nodesig <- ggtree::groupClade(.data = oak_nodesig,
                                  .node = nodes)

cols <- c("grey10",
          RColorBrewer::brewer.pal(n = 9, name = "Blues")[4:5],
          RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")[4:9],
          RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")[5:6],
          RColorBrewer::brewer.pal(n = 9, name = "Blues")[4:6],
          RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")[4],
          RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")[4],
          RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")[4],
          RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")[4],
          RColorBrewer::brewer.pal(n = 9, name = "Blues")[4:7])

# Plot tree with hot and cold clades
oak_ggtree <- ggtree::ggtree(oak_nodesig,
                             ggplot2::aes(color = group),
                             size = 1.3,
                             layout = "circular") +
  ggplot2::scale_color_manual(values = cols) +
  ggtree::geom_tiplab(color = "black",
                      offset = 15,
                      linewidth = 1.7,
                      fontface = "italic") +
  ggtree::geom_point2(ggplot2::aes(subset = (label %in% nodes)),
                      shape = 20, size = 5)

## PLOT HOT CLADES WITH INFO ON NO. AGRILUS SPP. HOSTED
# Add info on no. Agrilus spp. hosted
oak_ggtree_hosts <- oak_ggtree + ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(data = oak_host_observations,
                          geom = geom_bar,
                          mapping = aes(y = quercus_sp,
                                        x =  no_agrilus_spp),
                          offset = 0.15,
                          pwidth = 0.1,
                          # skip aggregation
                          stat = "identity",
                          # axis orientation
                          orientation = "y")

## PLOT HOT CLADES WITH INFO ON NO. AGRILUS SPP. HOSTED AND NO. GBIF ENTRIES
# Add no. GBIF entries
oak_ggtree_gbif <- oak_ggtree_hosts %<+% oak_gbif_entries +
  ggnewscale::new_scale_fill() +
  geom_tippoint(ggplot2::aes(x = x + 4, size = no.entries),
                shape = 21, fill = alpha("#86929C", .5),
                colour = "grey") +
  scale_size_continuous(range = c(0.1, 10),
                        limits = c(1, 500000),
                        breaks = c(1, 1000, 10000,
                                   100000, 500000)) +
  guides(color = "none") +
  theme(legend.position = c(0.9, 0.95))

oak_ggtree_gbif
