#### SET ENVIRONMENT ####
# Custom functions
source("R/functions.R")

# Seed
set.seed(24601)

# Libraries
library(brms)
library(foreach)

# Load data from previous script (brms models)
load("data/results/models.RData")

#### GENERATE LOOIC COMPARISONS DIVIDING MODELS INTO 3 GROUPS ####
# See https://discourse.mc-stan.org/t/understanding-looic/13409
# and https://mc-stan.org/loo/reference/loo-glossary.html#se-diff
# Lower LOOIC values indicate a better fit

# Compare loo results between models within each group
# NB, wiht waic(), I get a message suggesting to use loo instead: "p_waic
# estimates greater than 0.4. We recommend trying loo instead."
# Note that loo_compare_parallel_groups() has only been tested with a subset
# of the oak_models data, due to time constraints; behaviour and performance
# with large brmsfit objects, or many models per group have not yet been
# evaluated.
loo_comp_results <- loo_compare_parallel_groups(
  models_list = oak_models,
  n_groups = 2,
  save = FALSE,
  verbose = TRUE
)

# Explore results
loo_comp_results$grp1
loo_comp_results$grp2
loo_comp_results$grp3

# Save image
# save.image("data/results/models_loo_comp.RData")


#### GENERATE LOOIC COMPARISONS FOR BEST-PERFORMING MODELS ACROSS GROUPS ####
# Compare best-performing models from all groups
# These are models for which, in each group:
# * Their value overlaps with that of the model with the lowest LOOIC in their
# group, i.e. 0 is not contained in [eldp_diff - se_diff, eldp_diff + se_diff]
# * If the value does not overlap, then eldp_diff/se_diff < 2
# See https://discourse.mc-stan.org/t/understanding-looic/13409/3
# "'Say [eldp_diff/se_diff] was only twice as big??? would that mean there is a
# negligible difference between the models?' 'Not necessarily negligible
# difference, but certainly there is non-negligible uncertainty about the
# difference'"

# NB, models tested below correspond to the best performing models in the
# original analyses
keys <- paste0("mod",
               c(39, 41, 43, 46, 71, 72, 73, 75, 77, 78, 80, 83, 91, 105, 106,
                 107, 108, 109, 112, 116, 124, 125, 127, 131, 132))
oak_models_best <- oak_models[keys]

loo_oak_mods <- brms::loo(oak_models_best,
                          moment_match = TRUE,
                          reloo = TRUE,
                          compare = TRUE)

# Compare best-performing models from all groups which only contain variables
# with significant effects
# For all models in the final comparison above with non-significant variables,
# their respective model with only significant variables is also included above
keys <- paste0("mod",
               c(39, 41, 43, 46, 71, 72, 75, 77, 80, 91, 105, 108, 116))
oak_models_best_sig <- oak_models[keys]
loo_oak_mods <- brms::loo(oak_models_best_sig,
                          moment_match = TRUE,
                          reloo = TRUE,
                          compare = TRUE)

# Save image
# save.image("data/results/models_loo_comp.RData")
