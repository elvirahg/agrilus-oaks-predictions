#### SET ENVIRONMENT ####
# Custom functions
source("R/functions.R")

# Seed
set.seed(24601)


#### PREPARE DATA ####
interaction_data <- read.table("data/tmp/interaction_data.tsv",
                               header = TRUE,
                               sep = "\t")

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

# Create a phylo cov distance matrix (Hadfield & Nakagawa, 2010)
phylo_cov <- ape::vcv.phylo(oak_phylo)


#### FIT PRELIMINARY GEO DISTANCE MODELS ####
# Prepare model formulas
formulas_geodist <- list(
  oak_mod_geodist01 = as.formula("interaction ~ geo_min_dist_mean"),
  oak_mod_geodist02 = as.formula("interaction ~ geo_min_dist_median"),
  oak_mod_geodist03 = as.formula("interaction ~ geo_min_dist_mean_norm"),
  oak_mod_geodist04 = as.formula("interaction ~ geo_min_dist_median_norm"),
  oak_mod_geodist05 = as.formula("interaction ~ geo_min_dist_mean_log"),
  oak_mod_geodist06 = as.formula("interaction ~ geo_min_dist_median_log")
)

# Fit models
oak_models_geodist <- lapply(formulas_geodist, function(f) {
  brms::brm(
    formula = f,
    data = interaction_data,
    family = "bernoulli",
    iter = 7000,
    cores = 4,
    control = list(adapt_delta = 0.995, max_treedepth = 15),
    save_pars = brms::save_pars(all = TRUE)
  )
})
names(oak_models_geodist) <- names(formulas_geodist)

# Compute and compare LOO metrics to select best-performing geo distance metric
loo_list <- lapply(oak_models_geodist, brms::loo,
                   moment_match = TRUE,
                   reloo = TRUE)
loo_list_compared <- brms::loo_compare(loo_list)

# "geo_min_dist_median_norm" and "geo_min_dist_mean_norm" are the best fitting geo vars,
# and thus to be used in models below
loo_list_compared
#                   elpd_diff se_diff
# oak_mod_geodist04    0.0       0.0
# oak_mod_geodist03  -21.3       9.4
# oak_mod_geodist06  -29.0       6.9
# oak_mod_geodist05  -40.7       8.7
# oak_mod_geodist01 -152.6      17.0
# oak_mod_geodist02 -153.8      16.9


#### FIT FULL MODELS ####
# Vector with all variables to be included in model formulas
variables <- c("geo_min_dist_mean_norm",
               "phylo_dist_min",
               "gbif_entries",
               "(1 | gr(quercus_sp, cov = phylo_cov))",
               "(1 | agrilus_sp)")

# Generate formulas
formulas <- generate_formulas(
  variables = variables,
  response = "interaction",
  include_null = TRUE
)

# Substitute geo variables: to compare results using "geo_min_dist_mean_norm" vs
# "geo_min_dist_median_norm", duplicate every formula containing
# "geo_min_dist_mean_norm" and substitute it with "geo_min_dist_median_norm"
formulas <- add_formula_variations(formulas = formulas,
                                   pattern = "geo_min_dist_mean_norm",
                                   replacement = "geo_min_dist_median_norm")

# Substitute phylo variables: use same logic as above to compare performance of
# the different Quercus sp - Agrilus sp hosts "specific" phylogenetic distances,
# and check whether adding a random slope improves performance
formulas <- add_formula_variations(formulas = formulas,
                                   pattern = "phylo_dist_min",
                                   replacements = c("phylo_dist_mean",
                                                    "phylo_dist_min + phylo_dist_mean"))

# Substitute phylo_dist(s) + (1 | agrilus_sp) with a random slope, i.e
# phylo_dist(s) + (phylo_dist(s) | agrilus_sp)
new_formulas <- list()
for (formula in formulas) {
  # Turn to string
  formula_str <- paste(deparse(formula), collapse = "")
  formula_str <- gsub("\\s+", " ", formula_str)

  # Pattern to catch
  pattern <- "phylo_dist_min( \\+ phylo_dist_mean)?|phylo_dist_mean"

  # If "(1 | agrilus_sp)" and either "phylo_dist_min", "phylo_dist_mean",
  # or "phylo_dist_min + phylo_dist_mean" in formula
  if (grepl(pattern, formula_str) && grepl("agrilus_sp", formula_str)) {
    # Extract phylo.dist variables in formula
    phylo_dist_vars <- regmatches(formula_str, regexpr(pattern, formula_str))

    # Substitute (1 | agrilus_sp) with (phylo_dist(s) | agrilus_sp)
    new_formula_str <- sub("(1 | agrilus_sp)",
                           paste("(", phylo_dist_vars, " | agrilus_sp)",
                                 sep = ""),
                           formula_str,
                           fixed = TRUE)

    # Convert back to formula
    new_formula <- as.formula(new_formula_str)
    new_formulas <- c(new_formulas, list(new_formula))
  }
}
formulas <- c(formulas, new_formulas)

formulas <- sort_formulas(formulas = formulas,
                          order_by = c("number_variables",
                                       "has_random",
                                       "alphabetic"),
                          random_regex = c("\\([0-9a-zA-Z\\+ ]*\\| gr\\(quercus_sp, cov = phylo_cov\\)\\)",
                                           "\\([0-9a-zA-Z\\+ ]*\\| agrilus_sp\\)"),
                          f_names = "oak_mod_")
formulas

# Fit models
# Could be parallelised
oak_models <- lapply(formulas, function(f) {
  brms::brm(
    formula = f,
    data = interaction_data,
    family = "bernoulli",
    iter = 7000,
    cores = 4,
    control = list(adapt_delta = 0.995, max_treedepth = 15),
    data2 = list(phylo_cov = phylo_cov),
    save_pars = brms::save_pars(all = TRUE)
  )
})
names(oak_models) <- names(formulas)

# save.image("data/results/models.RData")
