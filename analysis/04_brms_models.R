##########################################################################################
####                           AGRILUS OAK HOSTS ANALYSES                          ######
####                           ELVIRA HERNANDEZ GUTIERREZ                          ######
##########################################################################################

####       SET ENVIRONMENT                                                          #####
# Set directory and clean environment
# setwd("~/oak_analyses/01_original_models")          # Whatever your working directory is

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
#####      PREPARE DATA                                                             #####
# i) Read in interaction data
interaction_data <- read.table("tmp/interaction_data_oak_hosts.txt",
                               header = T, sep = "\t")
nrow(interaction_data) == 236*32  # 236 oak spp.*32 Agrilus spp.
interaction_data$quercus.sp <- gsub("-", "", interaction_data$quercus.sp, fixed = T)


# ii) Add info on total no. of total GBIF entries for each oak species
# First, read (cleaned) GBIF geo info DF and only retain info for oaks
geo_info <- read.table("tmp/gbif_plants_clean.tsv",
                       header = T, sep = "\t")
geo_info <- geo_info[, c("plant.sp", "lon", "lat")]
geo_info <- geo_info[grepl("Quercus", geo_info$plant.sp),]

# Extract no. occurrences
oak_occurrences <- data.frame(table(geo_info$plant.sp))
colnames(oak_occurrences) <- c("quercus.sp", "gbif.entries")
oak_occurrences$quercus.sp <- as.character(oak_occurrences$quercus.sp)
oak_occurrences$quercus.sp <- gsub("-", "", oak_occurrences$quercus.sp, fixed = T)
oak_occurrences$quercus.sp <-  gsub("Quercus", "Q.", oak_occurrences$quercus.sp,
                                    fixed = TRUE)

# Append this info to interaction_data
interaction_data$ID <- 1:nrow(interaction_data)
interaction_data <- merge(interaction_data, oak_occurrences, by = "quercus.sp",
                          all = T)
# interaction_data[is.na(interaction_data$gbif.entries),]$gbif.entries <- 0
interaction_data <- interaction_data[order(interaction_data$ID), ]
interaction_data$ID <- NULL

rm(geo_info, oak_occurrences)

# iii) Retrieve info on Agrilus hosts
agrilus_hosts <- read.table("input/agrilus_quercus_hosts.txt",
                            header = F, sep = "\t")
colnames(agrilus_hosts) <- c("agrilus", "hosts")
agrilus_hosts <- rbind(agrilus_hosts,
                       read.table("tmp/non_oak_hosts.tsv",
                                  header = T, sep = "\t"))
agrilus_hosts <- agrilus_hosts[,c(2,1)]
colnames(agrilus_hosts) <- c("plant.sp", "agrilus.sp")
agrilus_hosts$agrilus.sp <- gsub("Agrilus", "A.", agrilus_hosts$agrilus.sp)

# Remove any hosts that are not in the phylogeny (only Q. gambelii) from agrilus_hosts
agrilus_hosts <- agrilus_hosts[!grepl("gambelii", agrilus_hosts$plant.sp),]

# IMPORTANT: I have also decided to remove any reported hosts that are non-native
# for those Agrilus species that have been introduced into new areas
# A. auroguttatus: introduced from AZ (USA) to CA (USA). In CA, reported novel larval
# hosts are Q. kelloggii (not in DS), Q. chrysolepis (not in DS) & Q. agrifolia (in DS)
# A. bilineatus: introduced from N America to Turkey (no new reported hosts so far)
# A. sulcicollis: introduced from Europe to N America, where Q. macrocarpa (in my DS)
# is a new reported larval host
# head(agrilus_hosts[agrilus_hosts$plant.sp == "Quercus macrocarpa" &
#                    agrilus_hosts$agrilus.sp == "A. sulcicollis",])

agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Quercus agrifolia"
                                 & agrilus_hosts$agrilus.sp == "A. auroguttatus")), ]
agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Quercus macrocarpa"
                                 & agrilus_hosts$agrilus.sp == "A. sulcicollis")), ]

# IMPORTANT: Let's also remove an larval hosts according to J&P(2014) (i.e., there's at
# least one reference with '!') but with a  confidence index < 3
# A. graminis: Corylus avellana, Ostrya carpinifolia, Euonymus europaeus & Castanea sativa
# A. obscuricollis: Ficus carica
# A. relegatus (reported as A. dualis)	Quercus robur
agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Ficus carica"
                                 & agrilus_hosts$agrilus.sp == "A. obscuricollis")), ]

agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Quercus robur"
                                 & agrilus_hosts$agrilus.sp == "A. relegatus")), ]

agrilus_hosts <- agrilus_hosts[-(which((agrilus_hosts$plant.sp == "Corylus avellana" |
                                        agrilus_hosts$plant.sp == "Ostrya carpinifolia" |
                                        agrilus_hosts$plant.sp == "Euonymus europaeus" |
                                        agrilus_hosts$plant.sp == "Castanea sativa") &
                                        agrilus_hosts$agrilus.sp == "A. graminis")), ]


# iv) Create phylogenetic cov matrix
# Read nexus tree created with phylocom nodesig (first tree: host status info)
oak_phylo <- read.nexus("input/quercus_nodesig_result.nex")[[1]]
# plot(oak_phylo)

# Edit node label information
oak_phylo$node.label <- gsub("'", "", oak_phylo$node.label)
oak_phylo$node.label <- gsub("N", "", oak_phylo$node.label)

# Replace "_" for " "
oak_phylo$tip.label <- gsub(pattern = "_", replacement = " ", oak_phylo$tip.label)
oak_phylo$tip.label <- gsub(pattern = "Quercus", replacement = "Q.",
                            oak_phylo$tip.label, fixed = TRUE)

# Create a phylo object with info on the relationship between species. This can
# be used to construct a cov matrix of species (Hadfield & Nakagawa, 2010)
oak_phylo_cov <- vcv.phylo(oak_phylo)
rm(oak_phylo)


# v) Add phylo distance info to interaction_data
# Create a new column for phylo distance info
interaction_data$phylo.dist.mean <- NA
interaction_data$phylo.dist.min <- NA

# For every possible Quercus - Agrilus interaction, find the **oak** host species of
# the current Agrilus species, and then find the mean phylogenetic distance of the
# focal Quercus species to the hosts, and normalise it. If the focal Quercus species
# happens to be the only host, then set the distance to the max distance
for(i in 1:nrow(interaction_data)) {
  ## Get the oak sp, the Agrilus sp.
  q <- interaction_data$quercus.sp[i]
  a <- interaction_data$agrilus.sp[i]

  ## Extract the hosts of the given Agrilus species
  h <- grep("Quercus", subset(agrilus_hosts, agrilus.sp == a)$plant.sp, value = TRUE)
  h <- gsub("Quercus", "Q.", h, fixed = TRUE)

  ## Remove the current oak species from this list if it happens to be a host
  h <- h[!(h %in% q)]

  ## If the current oak species is the only known host of this beetle, set the phylo
  ## distance to the max distance (using 'min' function because, before rescaling
  ## below, high numbers here mean things are close)
  if (length(h) == 0) {
    # print(paste(a, q))
    d.mean <- min(oak_phylo_cov)
    d.min <- min(oak_phylo_cov)
  }
  ## Otherwise, extract the oak_phylo_cov values (i.e. the phylo distance) for each
  ## q - h, and calculate the mean/min phylogenetic distance
  else {
    r <- which(rownames(oak_phylo_cov) == q)
    c <- which(colnames(oak_phylo_cov) %in% h)
    d.mean <- mean(oak_phylo_cov[r,c])
    d.min <- max(oak_phylo_cov[r,c])
  }
  interaction_data$phylo.dist.mean[i] <- d.mean
  interaction_data$phylo.dist.min[i] <- d.min
}

# Now normalise distances
# plot(phylo.dist.mean ~ phylo.dist.min, data = interaction_data)
interaction_data$phylo.dist.mean <- scales::rescale(-c(interaction_data$phylo.dist.mean))
interaction_data$phylo.dist.min <- scales::rescale(-c(interaction_data$phylo.dist.min))

# And compare
# plot(phylo.dist.mean ~ phylo.dist.min, data = interaction_data)



#
#####      RUN MODELS - PRELIMINARY MODELS (GEO DISTANCE ONLY)                      #####
# Which geographic distance metric performs best?

# Prepare model formulas (all possible variable combinations)
formulas <- list(
  oak_modA1 = "interaction ~ min.dist.mean",
  oak_modB1 = "interaction ~ min.dist.median",
  oak_modC1 = "interaction ~ min.dist.mean.norm",
  oak_modD1 = "interaction ~ min.dist.median.norm",
  oak_modE1 = "interaction ~ min.dist.mean.log",
  oak_modF1 = "interaction ~ min.dist.median.log")

## Run models
for (i in 1:length(formulas)){
  oak_mod <- brm(formulas[[i]],
                 data = interaction_data,
                 family = "bernoulli",
                 iter = 7000, cores = 4,
                 control = list(adapt_delta = 0.995,  max_treedepth = 15),
                 data2 = list(oak_phylo_cov = oak_phylo_cov),
                 save_pars = save_pars(all = TRUE)
  )
  assign(names(formulas[i]), oak_mod)
}

# Compare models and choose best-performing geographic distance metric
loo_geo_dist <- loo(oak_modA1, oak_modB1, oak_modC1,
                    oak_modD1, oak_modE1, oak_modF1,
                    save_psis = T, moment_match = T,
                    reloo = T, compare = T)
loo_geo_dist$diff

# Best-fitting models = min.dist.mean.norm & min.dist.median.norm
#             elpd_diff   se_diff
# oak_modC1    0.0        0.0
# oak_modD1   -4.7        6.0
# oak_modF1   -57.9       9.4
# oak_modE1   -65.3       9.7
# oak_modA1   -146.0      17.1
# oak_modB1   -147.7      17.1

# save.image("/results/oak_models.RData")



#
#####      RUN MODELS - ALL COMBINATIONS                                            #####
load("/results/oak_models.RData")

# A) PREPARE FORMULAS
# Vector with all variables that we're going to use (minus ""variations"")
variables <- c("min.dist.mean.norm", "phylo.dist.min", "gbif.entries",
               "(1 | gr(quercus.sp, cov = oak_phylo_cov))", "(1 | agrilus.sp)")

# All variable combinations
# NB, got this bit from https://stackoverflow.com/questions/40049313/generate-all-
# combinations-of-all-lengths-in-r-from-a-vector
formulas <- do.call("c",
                    lapply(seq_along(variables),
                           function(i) combn(variables, i, FUN = list)))
formulas <- lapply(formulas, paste, collapse = " + ")

# We established above that "min.dist.mean.norm" and "min.dist.median.norm" are the
# best fitting geographic variables - so let's only use there two for the subsequent
# models. Here, every time there is "min.dist.mean.norm" in our model formula, make
# another formula that's the exact same, but with "min.dist.median.norm" instead
# Note that with all 6 geographic distance metrics we'd have 308 models to run, as:
# sum(grepl("min.dist.mean.norm", formulas), na.rm = TRUE)*4 + length(formulas)
formulas.tmp <- list()
for (i in 1:length(formulas)) {
  formulas.tmp <-  c(formulas.tmp, formulas[[i]])
  if (grepl("min.dist.mean.norm", formulas[[i]])) {
    formulas.tmp <- c(formulas.tmp,
                      gsub("min.dist.mean.norm", "min.dist.median.norm",
                           formulas[[i]]))
  }
}
formulas <- formulas.tmp; rm(formulas.tmp)

# We want to compare the performance of the different Quercus sp. - Agrilus sp.
# "specific" phylogenetic distances, and check whether we should add a random slope
# effect to our model. Here, every time there is "phylo.dist.min" in our model formula,
# make another formula that's the exact same, but with "phylo.dist.mean" instead, and
# another one that's "phylo.dist.min + phylo.dist.mean"
formulas.tmp <- list()
for (i in 1:length(formulas)) {
  formulas.tmp <-  c(formulas.tmp, formulas[[i]])
  if (grepl("phylo.dist.min", formulas[[i]])) {
    formulas.tmp <-  c(formulas.tmp,
                       gsub("phylo.dist.min", "phylo.dist.mean",
                            formulas[[i]], fixed = TRUE))
    formulas.tmp <-  c(formulas.tmp,
                       gsub("phylo.dist.min", "phylo.dist.min + phylo.dist.mean",
                            formulas[[i]], fixed = TRUE))
  }
}
formulas <- formulas.tmp; rm(formulas.tmp)

# Here, every time there is "phylo.dist.min|phylo.dist.mean" plus an Agrilus random
# effect in our model formula, make another formula that's the exact same, but with
# a random slope instead
formulas.tmp <- list()
for (i in 1:length(formulas)) {
  formulas.tmp <-  c(formulas.tmp, formulas[[i]])
  if (grepl("phylo.dist.min|phylo.dist.mean", formulas[[i]]) &
      grepl("agrilus.sp", formulas[[i]])) {
    if(grepl("phylo.dist.min \\+ phylo.dist.mean", formulas[[i]])) {
      formulas.tmp <- c(formulas.tmp,
                        gsub("(1 | agrilus.sp)",
                             "(phylo.dist.min + phylo.dist.mean | agrilus.sp)",
                             formulas[[i]], fixed = TRUE))
    } else {
      v <- gsub(".*(phylo.dist.min|phylo.dist.mean).*", "\\1", formulas[[i]])
      formulas.tmp <- c(formulas.tmp,
                        gsub("(1 | agrilus.sp)",
                             paste("(", v, " | agrilus.sp)", sep = ""),
                             formulas[[i]], fixed = TRUE))
    }
  }
}
formulas <- formulas.tmp; rm(formulas.tmp)

# Add response variable to our formulas
formulas <- gsub("^", "interaction ~ ", formulas)

# Add intercept null model as the first one
formulas <- c(list("interaction ~ 1"), formulas)

# Rename our formulas with useful model names
names(formulas) <- paste("oak_mod", 1:length(formulas))
names(formulas) <- gsub(" ([0-9])$", "00\\1", names(formulas))
names(formulas) <- gsub(" ([0-9][0-9])$", "0\\1", names(formulas))
names(formulas) <- gsub(" ", "", names(formulas))



# B) RUN MODELS
# i) Run models one at a time
# THIS NEEDS TO BE INSIDE A JOBSCRIPT: ask for 133 array jobs, 2G/core, 4 cores, 240 h
# (the longest-running jobs took ca. 3 days to run)
i <- as.numeric(commandArgs(trailingOnly=TRUE))
oak_mod <- brm(formulas[[i]],
               data = interaction_data,
               family = "bernoulli",
               iter = 7000, cores = 4,
               control = list(adapt_delta = 0.995,  max_treedepth = 15),
               data2 = list(oak_phylo_cov = oak_phylo_cov),
               save_pars = save_pars(all = TRUE)
)

assign(names(formulas[i]), oak_mod)

# save.image(paste("/results/oak_models", names(formulas[i]), ".RData", sep = ""))
# END OF JOBSCRIPT

# ii) Save all models in one .RData file
files <- paste("results/", list.files(path = "results/"), sep = "")
files <- grep("[0-9].RData", files, value = TRUE)

lapply(files, load, .GlobalEnv)

load("/results/oak_models.RData")
# save.image("/results/oak_models.RData")



#