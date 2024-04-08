##########################################################################################
####                           AGRILUS OAK HOSTS ANALYSES                          ######
####                           ELVIRA HERNANDEZ GUTIERREZ                          ######
##########################################################################################

####       SET ENVIRONMENT                                                          #####
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
#####      CHECK LOOIC - 3 GROUPS                                                   ####
# See https://discourse.mc-stan.org/t/understanding-looic/13409
# and https://mc-stan.org/loo/reference/loo-glossary.html#se-diff
# Lower LOOIC values indicate a better fit
load("results/oak_models.RData")

# NB, wtih WAIC (waic())), I get a message suggesting to use loo instead:
# "p_waic estimates greater than 0.4. We recommend trying loo instead."

# A) MAKE THREE SETS FOR COMPARISON
# THIS NEEDS TO BE INSIDE A JOBSCRIPT USING 1 CORE (3 ARRAY JOBS)
# cat(ls(pattern = "oak_mod"), sep = ", ")
i <- as.numeric(commandArgs(trailingOnly=TRUE))

if (i == 1) {
  # Takes about 1 day to run, uses ca. 80 GB
  rm(list=setdiff(ls(),
                  c("oak_mod001", "oak_mod002", "oak_mod003", "oak_mod004", "oak_mod005",
                    "oak_mod006", "oak_mod007", "oak_mod008", "oak_mod009", "oak_mod010",
                    "oak_mod011", "oak_mod012", "oak_mod013", "oak_mod014", "oak_mod015",
                    "oak_mod016", "oak_mod017", "oak_mod018", "oak_mod019", "oak_mod020",
                    "oak_mod021", "oak_mod022", "oak_mod023", "oak_mod024", "oak_mod025",
                    "oak_mod026", "oak_mod027", "oak_mod028", "oak_mod029", "oak_mod030",
                    "oak_mod031", "oak_mod032", "oak_mod033", "oak_mod034", "oak_mod035",
                    "oak_mod036", "oak_mod037", "oak_mod038", "oak_mod039", "oak_mod040",
                    "oak_mod041", "oak_mod042", "oak_mod043", "oak_mod044", "i")))

  print("start")
  loo_oak_mods1 <- loo(oak_mod001, oak_mod002, oak_mod003, oak_mod004, oak_mod005,
                       oak_mod006, oak_mod007, oak_mod008, oak_mod009, oak_mod010,
                       oak_mod011, oak_mod012, oak_mod013, oak_mod014, oak_mod015,
                       oak_mod016, oak_mod017, oak_mod018, oak_mod019, oak_mod020,
                       oak_mod021, oak_mod022, oak_mod023, oak_mod024, oak_mod025,
                       oak_mod026, oak_mod027, oak_mod028, oak_mod029, oak_mod030,
                       oak_mod031, oak_mod032, oak_mod033, oak_mod034, oak_mod035,
                       oak_mod036, oak_mod037, oak_mod038, oak_mod039, oak_mod040,
                       oak_mod041, oak_mod042, oak_mod043, oak_mod044,
                       save_psis = TRUE, moment_match = TRUE,
                       reloo = TRUE, compare = TRUE)
} else if (i == 2) {
  # Takes about 2 days to run, uses ca. 80 GB
  rm(list=setdiff(ls(),
                  c("oak_mod045", "oak_mod046", "oak_mod047", "oak_mod048", "oak_mod049",
                    "oak_mod050", "oak_mod051", "oak_mod052", "oak_mod053", "oak_mod054",
                    "oak_mod055", "oak_mod056", "oak_mod057", "oak_mod058", "oak_mod059",
                    "oak_mod060", "oak_mod061", "oak_mod062", "oak_mod063", "oak_mod064",
                    "oak_mod065", "oak_mod066", "oak_mod067", "oak_mod068", "oak_mod069",
                    "oak_mod070", "oak_mod071", "oak_mod072", "oak_mod073", "oak_mod074",
                    "oak_mod075", "oak_mod076", "oak_mod077", "oak_mod078", "oak_mod079",
                    "oak_mod080", "oak_mod081", "oak_mod082", "oak_mod083", "oak_mod084",
                    "oak_mod085", "oak_mod086", "oak_mod087", "oak_mod088", "i")))

  print("start")
  loo_oak_mods2 <- loo(oak_mod045, oak_mod046, oak_mod047, oak_mod048, oak_mod049,
                       oak_mod050, oak_mod051, oak_mod052, oak_mod053, oak_mod054,
                       oak_mod055, oak_mod056, oak_mod057, oak_mod058, oak_mod059,
                       oak_mod060, oak_mod061, oak_mod062, oak_mod063, oak_mod064,
                       oak_mod065, oak_mod066, oak_mod067, oak_mod068, oak_mod069,
                       oak_mod070, oak_mod071, oak_mod072, oak_mod073, oak_mod074,
                       oak_mod075, oak_mod076, oak_mod077, oak_mod078, oak_mod079,
                       oak_mod080, oak_mod081, oak_mod082, oak_mod083, oak_mod084,
                       oak_mod085, oak_mod086, oak_mod087, oak_mod088,
                       save_psis = TRUE, moment_match = TRUE,
                       reloo = TRUE, compare = TRUE)
} else if (i == 3) {
  # Takes about 5 days to run, uses ca. 100 GB
  rm(list=setdiff(ls(),
                  c("oak_mod089", "oak_mod090", "oak_mod091", "oak_mod092", "oak_mod093",
                    "oak_mod094", "oak_mod095", "oak_mod096", "oak_mod097", "oak_mod098",
                    "oak_mod099", "oak_mod100", "oak_mod101", "oak_mod102", "oak_mod103",
                    "oak_mod104", "oak_mod105", "oak_mod106", "oak_mod107", "oak_mod108",
                    "oak_mod109", "oak_mod110", "oak_mod111", "oak_mod112", "oak_mod113",
                    "oak_mod114", "oak_mod115", "oak_mod116", "oak_mod117", "oak_mod118",
                    "oak_mod119", "oak_mod120", "oak_mod121", "oak_mod122", "oak_mod123",
                    "oak_mod124", "oak_mod125", "oak_mod126", "oak_mod127", "oak_mod128",
                    "oak_mod129", "oak_mod130", "oak_mod131", "oak_mod132", "i")))

  print("start")
  loo_oak_mods3 <- loo(oak_mod089, oak_mod090, oak_mod091, oak_mod092, oak_mod093,
                       oak_mod094, oak_mod095, oak_mod096, oak_mod097, oak_mod098,
                       oak_mod099, oak_mod100, oak_mod101, oak_mod102, oak_mod103,
                       oak_mod104, oak_mod105, oak_mod106, oak_mod107, oak_mod108,
                       oak_mod109, oak_mod110, oak_mod111, oak_mod112, oak_mod113,
                       oak_mod114, oak_mod115, oak_mod116, oak_mod117, oak_mod118,
                       oak_mod119, oak_mod120, oak_mod121, oak_mod122, oak_mod123,
                       oak_mod124, oak_mod125, oak_mod126, oak_mod127, oak_mod128,
                       oak_mod129, oak_mod130, oak_mod131, oak_mod132,
                       save_psis = TRUE, moment_match = TRUE,
                       reloo = TRUE, compare = TRUE)
}

# Save image
# save.image(paste0("results/oak_models_looic", i, ".RData"))
# print("saved")

# END OF JOBSCRIPT


# B) EXPLORE RESULTS
# i) Group 1
# Load data
load("results/oak_models_looic1.RData")

# Explore results
# loo_oak_mods1
#             elpd_diff  se_diff
# oak_mod043   0.0       0.0
# oak_mod044  -4.4       3.2     *
# oak_mod018  -9.9       4.9
# oak_mod019  -24.7      8.4
# oak_mod037  -35.7      9.3


# ii) Group 2
# Load data
load("results/oak_models_looic2.RData")

# Explore results
# loo_oak_mods2
#              elpd_diff  se_diff
# oak_mod045    0.0       0.0
# oak_mod083   -0.3       3.0
# oak_mod085   -1.1       3.0
# oak_mod084   -4.5       4.3
# oak_mod065   -7.9       4.6
# oak_mod061   -8.7       5.6     *
# oak_mod046   -11.1      5.5
# oak_mod086   -11.6      6.6
# oak_mod048   -11.9      5.5


# iii) Group 3
# Load data
load("results/oak_models_looic3.RData")

# Explore results
# loo_oak_mods3
#             elpd_diff  se_diff
# oak_mod126   0.0       0.0
# oak_mod106  -0.2       3.0
# oak_mod102  -0.3       3.1
# oak_mod122  -1.2       1.0
# oak_mod101  -3.0       4.6
# oak_mod105  -3.9       4.6
# oak_mod121  -4.5       4.0
# oak_mod104  -4.5       4.7
# oak_mod125  -4.9       3.9
# oak_mod103  -5.1       4.8
# oak_mod124  -5.4       3.7
# oak_mod123  -6.1       4.0
# oak_mod132  -7.0       6.0
# oak_mod112  -7.2       6.5
# oak_mod128  -8.6       5.9
# oak_mod108  -8.9       6.5
# oak_mod113  -10.6      5.6  *
# oak_mod127  -15.1      6.7
# oak_mod107  -15.2      7.0
# oak_mod131  -15.6      6.7



#
#####      CHECK LOOIC - BEST PERFOMING MODELS FROM EACH SET                        #####
# Ask for 100 GB, takes ca 1 day to run

# i) Using all best performing models from all 3 groups
load("results/oak_models.RData")

# The "best-performing" models I selected are models for which, in each group:
# *Their value overlaps with that of the lowest with the lowest LOOIC in their group,
# i.e. 0 not contained in [eldp_diff - se_diff, eldp_diff + se_diff]
# *If the value does not overlap, then eldp_diff/se_diff < 2??
# See https://discourse.mc-stan.org/t/understanding-looic/13409/3
# "'Say [eldp_diff/se_diff] was only twice as big??? would that mean there is a negligible
# difference between the models?' 'Not necessarily negligible difference, but certainly
# there is non-negligible uncertainty about the difference'"

# Compare models
loo_oak_mods <- loo(oak_mod018,             # out of curiosity
                    oak_mod043, oak_mod044,
                    oak_mod045, oak_mod061, oak_mod065,
                    oak_mod083, oak_mod084, oak_mod085,
                    oak_mod101, oak_mod102, oak_mod103, oak_mod104,
                    oak_mod105, oak_mod106, oak_mod108, oak_mod112,
                    oak_mod113, oak_mod121, oak_mod122, oak_mod123,
                    oak_mod124, oak_mod125, oak_mod126, oak_mod128,
                    oak_mod132,
                    save_psis = TRUE, moment_match = TRUE,
                    reloo = TRUE, compare = TRUE)
loo_oak_mods$diffs


# ii) Using only the best performing models from all 3 groups that only had variables
# with significant effects
# Ask for 80 GB, takes ca 1 day to run

# The reasoning behind this is that, so, for example - model 45 was one of the best
# performing models from group 3, but I am not including it here, because:
# Model 43: host.status ~ geo.mean + phylo.min + (1 | phylo)
# Model 45: host.status ~ geo.mean + phylo.min + phylo.mean + (1 | phylo)
# In model 43, all variables have a significant effect (i.e. they do not overlap with
# 0) but in model 45, phylo.mean does not. As it does not perform better than model 43,
# and when removing its non-significant variable it is == model 43, I decided to remove
# it from this table
# For all models in the final comparison above with non-significant  variables, their
# "respective" model with only significant variables is also included above (i.e.,
# there is a "model 43" for every "model 45"), so I decided to exclude the models with
# non-significant terms

# Compare models
loo_oak_mods_sig <- loo(oak_mod018, # out of curiosity
                        oak_mod043, oak_mod044,
                        oak_mod061, oak_mod065, oak_mod083, oak_mod084,
                        oak_mod101, oak_mod103, oak_mod108, oak_mod113,
                        oak_mod121, oak_mod123, oak_mod128,
                        save_psis = TRUE, moment_match = TRUE,
                        reloo = TRUE, compare = TRUE)
loo_oak_mods$diffs

# Save image
# save.image("results/oak_models.RData")



#