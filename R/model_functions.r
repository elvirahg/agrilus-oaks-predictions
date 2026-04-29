#' Generate all possible model formulas from a set of variables
#'
#' This function creates all possible combinations of predictor variables
#' for use in model formulae, optionally including a null (intercept-only)
#' model. The resulting list of formula objects can be used directly in
#' modelling functions.
#'
#' @param variables A character vector of predictor variables (right-hand
#' side terms).
#' @param response A single character string specifying the response variable
#' (left-hand side term).
#' @param include_null Logical; if `TRUE`, a null model (\code{response ~ 1})
#' is included as the first formula.
#'
#' @return A list of \code{formula} objects representing all model combinations.
#'
#' @examples
#' \dontrun{
#' vars <- c("geo.dist",
#'           "phylo.dist",
#'           "(1 | gr(plant_sp, cov = phylo_cov))")
#'
#' models <- generate_formulas(
#'   variables = vars,
#'   response = "interaction",
#'   include_null = TRUE
#' )
#' }
#'
#' @export
generate_formulas <- function(variables,
                              response,
                              include_null = TRUE) {
  # Checks
  if (!is.character(variables)) {
    stop("'variables' must be a character vector")
  }
  if (!is.character(response)) {
    stop("'response' must be a character vector")
  }
  if (!is.logical(include_null)) {
    stop("'inlude_null' must be logical")
  }

  # Generate all combinations of variables
  formulas <- do.call(
    "c",
    lapply(seq_along(variables), function(i) combn(variables, i, FUN = list))
  )

  # Collapse each combination into a single string
  formulas <- lapply(formulas, paste, collapse = " + ")

  # Add response variable
  formulas <- paste(response, "~", formulas)

  # Add a null model (intercept-only)
  if (include_null) {
    formulas <- c(paste(response, "~ 1"), formulas)
  }

  # Convert to formula objects
  formulas <- lapply(formulas, as.formula)

  formulas
}


#' Add variations of formula strings by pattern replacement
#'
#' Given a list of \code{formula} objects, this function generates additional
#' variations by substituting a specified `pattern` with one or more
#' `replacements`. The original formulas are preserved, and duplicates are
#' removed.
#'
#' @param formulas A list of \code{formula} objects.
#' @param pattern Character string specifying the exact text to replace.
#' Matching is done using fixed string matching (not regular expressions).
#' @param replacements Character vector of one or more replacement strings.
#'
#' @return A list of \code{formula} objects containing the original formulae and
#' any newly generated variants, with duplicates removed.
#'
#' @examples
#' \dontrun{
#' formulas <- list(as.formula("y ~ x + (1 | group)"))
#' pattern <- "(1 | group)"
#' replacements <- c("(x | group)", "(x + z | group)")
#' add_formula_variations(formulas, pattern, replacements)
#' # [1] "y ~ x + (1 | group)"
#' # [2] "y ~ x + (x | group)"
#' # [3] "y ~ x + (x + z | group)"
#' }
#'
#' @export
add_formula_variations <- function(formulas,
                                   pattern,
                                   replacements) {
  # Checks
  if (!is.list(formulas) || !all(sapply(formulas, inherits, "formula"))) {
    stop("'formulas' must be a list of formula objects")
  }
  if (!is.character(pattern)
      || !is.character(replacements)) {
    stop("'formulas', 'pattern', and 'replacements' must be characters")
  }

  # Convert formulas to character for substitution
  formula_strings <- vapply(formulas,
                            function(f) paste(deparse(f), collapse = ""),
                            character(1))

  out <- formula_strings
  for (rep in replacements) {
    out <- c(out, gsub(pattern, rep, formula_strings, fixed = TRUE))
  }

  # Remove duplicates and convert back to formula objects
  out <- unique(out)
  out <- lapply(out, as.formula)

  out
}


#' Sort a list of model formula objects
#'
#' This function sorts a list of model formula objects according to one or more
#' criteria applied to the right-hand side of each formula, ignoring the
#' response variable. It is useful for organising programmatically generated
#' sets of candidate models in a consistent order.
#'
#' The available sorting criteria are:
#' \itemize{
#'   \item \code{"number_variables"}: Number of fixed and random effects in the
#'          model.
#'   \item \code{"has_random"}: Presence of random-effect terms.
#'   \item \code{"alphabetic"}: Alphabetical order of the RHS formula string.
#' }
#'
#' @param formulas A list of \code{formula} objects to be sorted.
#' @param order_by A character vector specifying the sorting priority. Allowed
#'   values are \code{"number_variables"}, \code{"has_random"}, and
#'   \code{"alphabetic"}. Defaults to
#'   \code{c("number_variables", "has_random", "alphabetic")}.
#' @param random_regex Optional character vector of regular expressions used to
#'   identify random-effect terms when counting variables when `order_by`
#'   contains "number_variables". If not provided, random-effect terms may be
#'   miscounted.
#' @param f_names Optional character prefix for naming the sorted formulas. If
#'   not \code{NULL}, each formula in the output list will be named sequentially
#'   using this prefix and a zero-padded index (e.g. \code{"mod_001"},
#'   \code{"mod_002"}, ...).
#'
#' @return A named list of \code{formula} objects, sorted and optionally renamed
#'   according to the specified criteria.
#'
#' @examples
#' \dontrun{
#' formulas <- list(
#'   as.formula("resp ~ b"),
#'   as.formula("resp ~ a + (1 | e)"),
#'   as.formula("resp ~ a"),
#'   as.formula("resp ~ b + a"),
#'   as.formula("resp ~ b + a + (d | e)"),
#'   as.formula("resp ~ 1")
#' )
#'
#' sort_formulas(formulas)
#' }
#'
#' @export
sort_formulas <- function(formulas,
                          order_by = c("number_variables",
                                       "has_random",
                                       "alphabetic"),
                          random_regex = NULL,
                          f_names = "mod_") {
  # Checks
  if (!is.list(formulas) || !all(sapply(formulas, inherits, "formula"))) {
    stop("'formulas' must be a list of formula objects")
  }
  if (!is.character(order_by)) {
    stop("'order_by' must be a character vector")
  }
  if (!is.character(random_regex) && !is.null(random_regex)) {
    stop("'order_by' must be a character vector")
  }
  if (!is.character(f_names) && !is.null(f_names)) {
    stop("'f_names' must be NULL or a character vector")
  }

  # Convert formulas to character and isolate the right-hand side
  formula_strings <- vapply(formulas, function(f) {
    f_str <- paste(deparse(f), collapse = "")
    sub(".*~", "", f_str) |> trimws()
  }, character(1))

  # Initialise keys list
  keys <- list()

  # Split by '+' outside parentheses
  if ("number_variables" %in% order_by) {
    formula_rd <- Reduce(
      function(x, pattern) {
        gsub(pattern,
             "random",
             x,
             perl = TRUE)
      },
      random_regex,
      init = formula_strings
    )
    n_vars <- sapply(formula_rd, function(f) {
      length(strsplit(f,
                      " \\+ ",
                      perl = TRUE)[[1]])
    })

    # Add to list
    keys$number_variables <- n_vars
  }

  # Alphabetical ordering
  if ("alphabetic" %in% order_by) {
    keys$alphabetic <- unlist(formula_strings)
  }

  # Does formula contain randon effects
  if ("has_random" %in% order_by) {
    keys$has_random <- grepl("\\(", formula_strings)
  }

  # Order formulas
  sort_args <- lapply(order_by, function(k) keys[[k]])
  out_formulas <- formulas[do.call(order, sort_args)]

  # Name and return fomulas
  if (!is.null(f_names)) {
    names(out_formulas) <- sprintf("%s%03d", f_names, seq_along(out_formulas))
  }
  out_formulas
}


#' Run grouped LOO model comparison with optional parallelisation
#'
#' Performs leave-one-out cross-validation (LOO) for a list of `brmsfit` models.
#' Models can be split into groups, processed in parallel, and optionally saved.
#'
#' @param models_list Named list of `brmsfit` objects to compare.
#' @param n_groups Integer specifying number of groups to split the models into.
#' Optional if `model_groups` is supplied. Defaults to `NULL`. If used, the
#' number of parallel cores spawned equals `n_groups`.
#' @param model_groups Optional named list of character vectors, each containing
#' model names to define custom groups. If provided, `n_groups` must be set to
#' NULL.
#' @param save Logical; if `TRUE`, each group's `loo_compare` result is saved as
#' an `.rds` file in `save_path`. Defaults to `FALSE`.
#' @param save_path Character scalar specifying the directory in which to save
#' results if `save = TRUE`. Will be created if it does not exist. Defaults to
#' the current working directory.
#' @param verbose Logical; if `TRUE`, prints group assignment. Defaults to
#' `TRUE`.
#' @param ... Additional arguments passed directly to `brms::loo()`. For
#' example, `moment_match`, or `reloo` can be overridden.
#'
#' @return A named list of length `n_groups`, where each element contains the
#' corresponding `loo_compare` object for that group of models.
#'
#' @details
#' Models can be divided into groups in one of two ways:
#' 1. By specifying `n_groups`, in which case models are split automatically.
#' 2. By providing a named list `model_groups` specifying the exact models in
#' each group. This allows custom groupings. When `model_groups` is provided,
#' `n_groups` is ignored.
#'
#' The number of parallel cores used internally is equal to the number of
#' groups. Each group is processed sequentially, computing `brms::loo()` for
#' each model and then comparing with `brms::loo_compare()`. If `save = TRUE`,
#' the `loo_compare` result for each group is saved as
#' `loo_comp_<group_name>.RData` in `save_path`.
#'
#' This function has so far been tested only with relatively small datasets.
#' Behaviour and performance with very large `brmsfit` objects, or many models
#' per group have not yet been evaluated.
#'
#' @examples
#' \dontrun{
#' results <- run_grouped_loo(
#'   models_list = models,
#'   n_groups = 2
#' )
#' }
#'
#' @import brms
#' @import parallel
#' @import doParallel
#' @import foreach
#' @export
loo_compare_parallel_groups <- function(models_list,
                                        n_groups = NULL,
                                        model_groups = NULL,
                                        save = FALSE,
                                        save_path = getwd(),
                                        verbose = TRUE,
                                        ...) {
  # Checks
  if (!is.list(models_list)
      || !all(vapply(models_list, inherits, logical(1), "brmsfit"))) {
    stop("'models_list' must be a list of brmsfit objects")
  }
  if (!is.null(n_groups)) {
    if (!is.numeric(n_groups)
        || n_groups != as.integer(n_groups)
        || length(n_groups) != 1
        || n_groups > length(models_list)) {
      stop(paste("'n_groups' must be a single integer"))
    }
    if (n_groups > length(models_list)) {
      stop(paste("'n_groups' must be => length(models_list)"))
    }
  }
  if (!is.null(model_groups)) {
    if (!is.list(model_groups) || !all(sapply(model_groups, is.character))) {
      stop("'model_groups' must be a named list of character vectors")
    }
    if (!setequal(unlist(model_groups), names(models_list))) {
      stop("All models in 'models_list' must appear once in 'model_groups'")
    }
  }
  if ((!is.null(model_groups) && !is.null(n_groups))
      || (is.null(model_groups) && is.null(n_groups))) {
    stop("Either 'n_groups' or 'model_groups' must be provided")
  }
  if (!is.character(save_path) || length(save_path) != 1) {
    stop("'save_path' must be a single character string")
  }

  # Check directories if saving
  if (save && !dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  # Determine groups (useful if too many models to compare at once)
  if (is.null(model_groups)) {
    model_groups <- split_vector_into_lists(
      x = names(models_list),
      n_groups = n_groups
    )
    if (verbose) cat("Models divided into", length(model_groups), "group(s):\n")
    if (verbose) print(model_groups)
  }

  # I think it may be == ceiling(length(model_list)/n_groups) %% length(model_list)
  if (any(sapply(model_groups, length) == 1)) {
    stop("Cannot have groups of length 1, please use a lower 'n_groups' value")
  }

  results <- list()

  # Prepare per-group model lists to avoid huge exports
  group_model_lists <- lapply(model_groups, function(names) models_list[names])

  # Set up parallel environment
  cl <- parallel::makeCluster(n_groups)
  doParallel:::registerDoParallel(cl)

  # Export brms
  my_libs <- .libPaths()
  parallel::clusterExport(cl, "my_libs", envir = environment())
  parallel::clusterEvalQ(cl, .libPaths(my_libs))
  parallel::clusterEvalQ(cl, library(brms))

  # Loop over groups
  results <- foreach::foreach(grp = names(group_model_lists),
                              .packages = "brms",
                              .verbose = verbose,
                              .export = c("group_model_lists",
                                          "save",
                                          "save_path")) %dopar% {

    group_models <- group_model_lists[[grp]]

    # Compute LOO sequentially within this group
    loo_objs <- lapply(names(group_models), function(nm) {
      brms::loo(group_models[[nm]],
                moment_match = TRUE,
                reloo = TRUE,
                ...)
    })
    names(loo_objs) <- names(group_models)

    # Compare LOO objects
    loo_comp <- brms::loo_compare(loo_objs)

    # Optional save
    if (save) {
      save_file <- file.path(save_path, paste0("loo_comp_", grp, ".rds"))
      saveRDS(loo_comp, file = save_file)
    }

    # # Return LOO comparison
    loo_comp
  }

  # Name the list
  names(results) <- names(model_groups)

  # Stop cluster
  parallel::stopCluster(cl)

  results
}


#' Split character vector into evenly sized groups
#'
#' Divides a character vector into a specified number of approximately
#' equal-sized groups. If the total number of models is not divisible by
#' the number of groups,
#' group sizes will differ by at most one.
#'
#' @param x A character vector containing the items to be split.
#' @param n_groups An integer giving the number of groups to create.
#'
#' @return A named list of length \code{n_groups}, where each element contains
#' the subset of \code{x} belonging to that group. Group names are assigned as
#' \code{"grp1"}, \code{"grp2"}, and so on.
#'
#' @keywords internal
split_vector_into_lists <- function(x,
                                    n_groups) {

  # Checks
  if (!is.character(x)) {
    stop("'x' must be a character vector")
  }
  if (!is.numeric(n_groups)
      || n_groups != as.integer(n_groups)
      || length(n_groups) != 1) {
    stop("'n_groups' must be a single integer")
  }
  if (length(x) < n_groups) {
    stop("Error: length(x) < n_groups")
  }

  # Calculate number of models
  n_x <- length(x)

  # Group membership
  if (n_groups == 1) {
    groups <- rep(1, n_x)
  } else {
    groups <- cut(seq_len(n_x), breaks = n_groups, labels = FALSE)
  }

  # Split into a named list
  x_groups <- split(x, groups)
  names(x_groups) <- paste0("grp", seq_len(n_groups))

  x_groups
}


#' "Leave-one-observation-out"-like predictions by zeroing each positive case in
#' a brms model
#'
#' This function performs the follow: for every row in the model data where
#' `interaction_col == 1`, the value of that row is temporarily set to 0, and
#' the model is re-updated. It then computes the predictor (on the specified
#' `scale`) for the modified observation. This can be used to quantify how
#' individual positive observations contribute to the fitted model.
#'
#' Computation is parallelised across the number of workers specified by
#' `cores`.
#'
#' @param model A `brmsfit` object containing an attached data list
#' (i.e., `model$data` must exist).
#' @param interaction_col A single character string giving the name of the
#' binary column to perturb when generating 1-to-0 predictions. Defaults to
#' 'interaction'.
#' @param cores Integer. Number of parallel workers to use. Defaults to 2.
#' @param iter Integer. Total number of sampling iterations for each refitted
#' model. Defaults to 1000.
#' @param warmup Integer. Number of warmup iterations for each refitted model.
#' Defaults to 500.
#' @param chains Integer. Number of chains to run for each refitted model.
#' Defaults to 1.
#' @param recompile Logical. Passed to `brms::update()` to control whether
#' the underlying Stan model is recompiled on each iteration.
#' @param scale A single character string passed to `brms::fitted()`,
#' indicating the scale of the response variable. Defaults to 'linear'.
#'
#' @return A data frame with one row per positive observation. Columns:
#' #' \itemize{
#'   \item \code{row_no}: the index of the row modified
#'   \item \code{pred}: the predicted estimate for the modified row
#' }
#'
#' @details
#' This approach is similar in spirit to leave-one-out cross-validation,
#' except that the row is not removed; instead, its value in `interaction_col`
#' is set to 0 and the model is refit to estimate the corresponding
#' counterfactual prediction.
#'
#' @examples
#' \dontrun{
#' results <- predict_as_zero_loo(
#'   model = my_brms_model,
#'   interaction_col = "interaction",
#'   cores = 4
#' )
#' }
#'
#' @import brms
#' @import foreach
#' @import parallel
#' @import doParallel
#' @export
predict_as_zero_loo <- function(model,
                                interaction_col = "interaction",
                                cores = 2,
                                iter = 1000,
                                warmup = 500,
                                chains = 1,
                                recompile = FALSE,
                                scale = "linear") {
  # Checks
  if (!inherits(model, "brmsfit") || is.null(model$data)) {
    stop("'model' must be a 'brmsfit' object with a data slot (model$data)")
  }
  if (!is.character(interaction_col) || length(interaction_col) != 1) {
    stop("'interaction_col' must be a single character string")
  }
  if (!(interaction_col %in% names(model$data))) {
    stop(paste0("'", interaction_col, "' is not a column in model$data"))
  }
  if (!is.numeric(cores) || length(cores) != 1 || cores != as.integer(cores)
      || cores < 1) {
    stop("'cores' must be a single integer >= 1")
  }
  if (!is.numeric(iter) || length(iter) != 1 || iter != as.integer(iter)
      || iter < 1) {
    stop("'iter' must be a single integer >= 1")
  }
  if (!is.numeric(warmup) || length(warmup) != 1 || warmup != as.integer(warmup)
      || iter < 0) {
    stop("'warmup' must be a single non-negative integer")
  }
  if (!is.numeric(chains) || length(chains) != 1 || chains != as.integer(chains)
      || iter < 0) {
    stop("'chains' must be a single integer >= 1.")
  }
  if (!is.logical(recompile) || length(recompile) != 1) {
    stop("'recompile' must be a single logical value")
  }
  if (!is.character(scale) || length(scale) != 1) {
    stop("'scale' must be a single character string")
  }

  # Identify positive observations
  pos_int <- which(model$data[[interaction_col]] == 1)
  cat(paste("Number of positive interactions to be tested:",
            length(pos_int), "\n"))

  # Set up parallel workers
  cl <- parallel::makeCluster(cores)
  doParallel:::registerDoParallel(cl)

  # Make sure workers know where packages are
  libs <- .libPaths()
  parallel::clusterExport(cl, "libs", envir = environment())
  parallel::clusterEvalQ(cl, .libPaths(libs))
  parallel::clusterEvalQ(cl, library(brms))

  # Run loop
  results <- foreach::foreach(obs = pos_int,
                              .combine = rbind,
                              .packages = "brms") %dopar% {

    # Modify data row
    new_data <- model$data
    new_data[[interaction_col]][obs] <- 0

    # Refit
    mod_obs <- update(model,
                      newdata = new_data,
                      iter = iter,
                      chains = chains,
                      warmup = warmup,
                      recompile = recompile,
                      refresh = 0)

    # Predict linear predictor for that row
    pred <- fitted(mod_obs,
                   newdata = new_data[obs, ],
                   scale = scale)[, "Estimate"]

    data.frame(row_no = obs,
               pred = pred,
               row.names = NULL)
  }

  parallel::stopCluster(cl)

  results
}


#' K-fold cross-validated predictions per observation from a brms model
#'
#' Computes out-of-fold predictions for each observation using K-fold
#' cross-validation via \code{brms::kfold()}, and summarises posterior draws
#' into a single value per observation.
#'
#' For each observation, predictions are obtained from the model fit that did
#' not include that observation in training. Posterior expected predictions are
#' averaged across draws to yield a mean predicted probability, which can
#' optionally be transformed to log-odds.
#'
#' @param model A fitted \code{brmsfit} model object.
#' @param k Integer. Number of folds for K-fold cross-validation. Must be a
#' positive integer. Defaults to 5.
#' @param method Character string passed to \code{brms::kfold_predict()}.
#' Defaults to \code{"posterior_epred"}.
#' @param type Character string indicating the output scale:
#' \itemize{
#'   \item \code{"prob"}: return predicted probabilities
#'   \item \code{"lodds"}: return log-odds (logit-transformed probabilities)
#' }
#'
#' @return A numeric vector of length equal to the number of observations in the
#' original dataset, containing either predicted probabilities or log-odds.
#'
#' @examples
#' \dontrun{
#' fit <- brms::brm(y ~ x1 + x2, data = dat, family = bernoulli())
#'
#' # Predicted probabilities
#' kfold_pred <- kfold_predict_observations(fit, type = "prob")
#'
#' # Log-odds
#' kfold_lodds <- kfold_predict_observations(fit, type = "lodds")
#' }
#'
#' @import brms
#' @export
kfold_predict_observations <- function(model,
                                       k = 5,
                                       method = "posterior_epred",
                                       type = c("prob", "lodds")) {

  # Checks
  if (!inherits(model, "brmsfit") || is.null(model$data)) {
    stop("'model' must be a 'brmsfit' object with a data slot (model$data)")
  }
  if (!is.numeric(k) || k < 1 || k != as.integer(k)) {
    stop("'k' must be a positive integer")
  }
  if (!is.character(method) || length(method) != 1) {
    stop("'method' must be a single character string")
  }
  if ((length(type) != 1) || !(type %in% c("prob", "lodds"))) {
    stop("'type' must be 'prob', or 'lodds'")
  }

  # Run k-fold CV, storing the fitted models (save_fits = TRUE) to be able to
  # generate predictions
  kfold <- brms::kfold(model, K = k, save_fits = TRUE)

  # For each observation, use the model where that observation was held out, and
  # return posterior draws of predictions (yrep)
  preds <- brms::kfold_predict(kfold, method = method)

  # Compute one predicted probability per observation
  prob <- colMeans(preds$yrep)

  # Return prob or log-odds result according to type
  if (type == "prob") {
    message("Returning predicted probabilities")
    prob
  } else if (type == "lodds") {
    message("Returning predicted log-odds")
    lodds <- qlogis(prob)
    lodds
  }

}


#' Summarise threshold-based metrics for model and cross-validation predictions
#'
#' Computes simple classification summaries for a set of prediction vectors
#' evaluated against a common threshold. Metrics are calculated for model
#' predictions, cross-validation method 1, and cross-validation method 2
#' }
#'
#' For each method, the following are returned:
#' \itemize{
#'   \item \code{true_positive_rate}: percentage of reported positive cases
#'   (where \code{reported_status == 1}) with predictions above the threshold
#'   \item \code{count_over_threshold}: number of predictions above the
#'   threshold
#'   \item \code{percent_over_threshold}: percentage of all predictions above
#'   the threshold
#' }
#'
#' @param cv_method_1 Numeric vector. Predictions from the first
#' cross-validation method.
#' @param cv_method_2 Numeric vector. Predictions from the second
#' cross-validation method. Must be the same length and of the same scale
#' as \code{cv_method_1}.
#' @param model Numeric vector of model predictions used for
#' threshold-based classification. Must be the same length and of the same
#' scale as \code{cv_method_1}.
#' @param reported_status Numeric, factor, or logical vector indicating observed
#' status for each observation. Values equal to 1 are treated as "reported
#' positive", values equal to 0 as "reported negative". Must be the same length
#' as \code{cv_method_1}.
#' @param threshold Numeric. Threshold used to classify predictions as positive
#' or negative. Must be on the same scale as \code{cv_method_1}.
#' @param names Named character vector of length 3 giving row names for the
#' output. Must have names \code{cv_method_1}, \code{cv_method_2}, and
#' \code{model}. Values are used as row labels in the output.
#'
#' @return A \code{data.frame} with three rows (one per method) and columns:
#' \code{true_positive_rate}, \code{percent_over_threshold}, and
#' \code{count_over_threshold}.
#'
#' @seealso \code{\link{compute_summary_metrics}}
#'
#' @examples
#' \dontrun{
#' cv_summary_table(
#'   cv_method_1 = loo_lodds,
#'   cv_method_2 = kfold_lodds,
#'   model_predictions = predictions$prediction_lodds,
#'   reported_status = predictions$host_status,
#'   threshold = thr_intercepts["threshold"]
#' )
#' }
#'
#' @export
cv_summary_table <- function(cv_method_1,
                             cv_method_2,
                             model,
                             reported_status,
                             threshold,
                             names = c(cv_method_1 = "cv_method_1",
                                       cv_method_2 = "cv_method_2",
                                       model = "model")) {

  # Checks
  if (!is.numeric(cv_method_1)) {
    stop("'cv_method_1' must be a numeric vector")
  }
  if (!is.numeric(cv_method_2)) {
    stop("'cv_method_2' must be a numeric vector")
  }
  if (length(cv_method_2) != length(cv_method_1)) {
    stop("'cv_method_1' and 'cv_method_2' must have the same length")
  }
  if (!(is.logical(reported_status)
        || is.factor(reported_status)
        || is.numeric(reported_status)
        || all(reported_status %in% c(0, 1)))) {
    stop("'reported_status' must be logical, factor, or binary (0/1) numeric")
  }
  if (length(reported_status) != length(cv_method_1)) {
    stop("'reported_status' must match length of 'cv_method_1")
  }
  if (!is.numeric(model)) {
    stop("'model_predictions' must be a numeric vector")
  }
  if (length(model) != length(cv_method_1)) {
    stop("'model_predictions' must match length of 'cv_method_1")
  }
  if (!is.numeric(threshold) || length(threshold) != 1) {
    stop("'model_predictions' must be a single numeric value")
  }
  if (!is.character(names)
      || length(names) != 3
      || is.null(names(names))
      || !setequal(names(names), c("cv_method_1", "cv_method_2", "model"))) {
    stop("'names' must be a named character vector of length 3 with names:",
         "'cv_method_1', 'cv_method_2', 'model'")
  }

  # Build table
  metrics_summary <- rbind(
    compute_summary_metrics(predictions = model,
                            reported_status = reported_status,
                            threshold = threshold),
    compute_summary_metrics(predictions = cv_method_1,
                            reported_status = reported_status,
                            threshold = threshold),
    compute_summary_metrics(predictions = cv_method_2,
                            reported_status = reported_status,
                            threshold = threshold)
  )
  rownames(metrics_summary) <- names

  metrics_summary
}


#' Compute threshold-based performance metrics for a numeric prediction vector
#'
#' Internal helper function to summarise classification performance given a
#' threshold on the log-odds scale.
#'
#' Computes:
#' \itemize{
#'   \item \code{true_positive_rate}: percentage of reported positive cases
#'   (where \code{reported_status == 1}) with predictions above the threshold
#'   \item \code{count_over_threshold}: number of predictions above the
#'   threshold
#'   \item \code{percent_over_threshold}: percentage of all predictions above
#'   the threshold
#' }
#'
#' @param predictions Numeric vector of predictions (on the same scale as
#' \code{threshold}).
#' @param threshold Numeric. Threshold on the log-odds scale used to classify
#' predictions as positive or negative.
#' @param reported_status Numeric, factor, or logical vector indicating observed
#' status for each observation. Values equal to 1 are treated as "reported
#' positive"; all other values are treated as "reported negative".
#' @return A one-row \code{data.frame} containing the computed metrics.
#'
#' @keywords internal
compute_summary_metrics <- function(predictions,
                                    reported_status,
                                    threshold) {
  metrics <- data.frame(
    true_positive_rate = mean(predictions[reported_status == 1]
                              > threshold) * 100,
    count_over_threshold = sum(predictions > threshold),
    percent_over_threshold = mean(predictions > threshold) * 100
  )

  metrics
}
