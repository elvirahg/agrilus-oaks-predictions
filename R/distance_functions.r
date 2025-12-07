#' Create an expanded host – hosted interaction data frame
#'
#' This function takes a data frame of known host – hosted "postive"
#' interactions and produces an expanded data frame containing all
#' possible combinations between specified host taxon members and
#' hosted taxa. Each combination is annotated with a binary column
#' (`interaction`) indicating whether the pair is present in the known
#' interactions.
#'
#' Optionally, users can specify subsets of host or hosted taxa to include;
#' otherwise, all unique taxa in the known interactions are used.
#' The function checks that taxa in `known_interactions` are present in the
#' supplied subsets and issues informative warnings or errors accordingly.
#'
#' @param known_interactions A data frame containing known host – hosted
#' interaction pairs.
#' @param host_col A string giving the column name in `known_interactions`
#' that contains host taxa.
#' @param hosted_col A string giving the column name in `known_interactions`
#' that contains hosted taxa.
#' @param host_taxa_include Optional vector of host taxon members to include
#' in the expanded data frame. If `NULL`, all unique host taxa in
#' `known_interactions[[host_col]]` are used.
#' @param hosted_taxa_include Optional vector of hosted taxa to include in
#' the expanded data frame. If `NULL`, all unique hosted taxa in
#' `known_interactions[[hosted_col]]` are used.
#' @param verbose Logical; if `TRUE`, prints warnings when taxa in
#' `known_interactions` are missing from the supplied include lists.
#' Defaults to `TRUE`.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{host}{Host taxon name}
#'   \item{hosted}{Hosted taxon name}
#'   \item{interaction}{Integer (0/1) indicating whether the pair occurs in
#'   `known_interactions`}
#' }
#'
#' @examples
#' \dontrun{
#' interaction_df <- create_interaction_df(
#'   known_interactions = agrilus_hosts,
#'   host_col = "plant_sp",
#'   hosted_col = "agrilus_sp"
#' )
#' }
#'
#' @export
create_interaction_df <- function(known_interactions,
                                  host_col,
                                  hosted_col,
                                  host_taxa_include = NULL,
                                  hosted_taxa_include = NULL,
                                  verbose = TRUE) {
  # Checks
  if (!is.data.frame(known_interactions)) {
    stop("'known_interactions' must be a data frame")
  }
  if (!all(c(host_col, hosted_col) %in% colnames(known_interactions))) {
    stop("Specified column names not found in 'known_interactions'")
  }
  if (!is.null(host_taxa_include) && !is.character(host_taxa_include)) {
    stop("'host_taxa_include' must be a character vector, or NULL")
  }
  if (!is.null(hosted_taxa_include) && !is.character(hosted_taxa_include)) {
    stop("'hosted_taxa_include' must be a character vector, or NULL")
  }

  # Determine host-taxon species to include
  if (is.null(host_taxa_include)) {
    host_taxa_include <- unique(known_interactions[[host_col]])
  } else {
    # Issue warning/error if known hosts species missing in host_taxa_include
    known_hosts <- unique(known_interactions[[host_col]])
    missing_hosts <- setdiff(known_hosts, host_taxa_include)

    # Issue warnings or errors if necessary
    if (length(missing_hosts) == length(known_hosts)) {
      stop("No hosts in known_interactions are present in host_taxa_include")
    } else if (length(missing_hosts) > 0 && verbose) {
      warning(paste("Some hosts in known_interactions are not present in",
                    "host_taxa_include: ",
                    paste(missing_hosts, collapse = ", ")))
    }
  }

  # Determine hosted species to include
  if (is.null(hosted_taxa_include)) {
    hosted_taxa_include <- unique(known_interactions[[hosted_col]])
  } else {
    # Issue warning/error if known hosts species missing in hosted_taxa_include
    known_hosted <- unique(known_interactions[[hosted_col]])
    missing_hosted <- setdiff(known_hosted, hosted_taxa_include)

    # Issue warnings or errors if necessary
    if (length(missing_hosted) == length(known_hosted)) {
      stop("No hosts in known_interactions are present in hosted_taxa_include")
    } else if (length(missing_hosted) > 0 && verbose) {
      warning(paste("Some hosted taxa in known_interactions are not present in",
                    "hosted_taxa_include: ",
                    paste(missing_hosted, collapse = ", ")))
    }
  }

  # Create all combinations
  expanded <- expand.grid(setNames(list(host_taxa_include,
                                        hosted_taxa_include),
                                   c(host_col, hosted_col)),
                          stringsAsFactors = FALSE)
  expanded <- as.data.frame(expanded, stringsAsFactors = FALSE)
  attr(expanded, "out.attrs") <- NULL

  # Mark which combinations are known interactions
  expanded$interaction <- as.integer(
    paste(expanded[[host_col]], expanded[[hosted_col]]) %in%
      paste(known_interactions[[host_col]], known_interactions[[hosted_col]])
  )

  expanded
}


#' Calculate mean and median minimum geographic distance metrics
#'
#' This function computes distance-based metrics for each host-taxon-member –
#' hosted species interaction in `expanded_interaction_df`. For each species in
#' the host taxon, it calculates the minimum distance to all (other) hosts of
#' the same hosted species, optionally subsampling the coordinate data to
#' improve performance. Distances are precomputed in a lookup table to speed up
#' repeated queries, and the resulting metrics can be further transformed
#' (e.g. normalised or log-scaled) for downstream modelling.
#'
#' @param expanded_interaction_df Data frame of all possible host-taxon-member
#' – hosted species interactions to evaluate. Must include columns for host and
#' hosted taxa.
#' @param coords_df Data frame of location occurrences for all host-taxon-member
#' species. Column names must be "lon", "lat". Must include columns for host
#' taxa, longitude (`lon`), and latitude (`lat`).
#' @param known_interactions_df Data frame of known "positive" interactions.
#' Must include columns for host and hosted taxa.
#' @param host_taxon_col Character. Column name indicating the host taxon in all
#' input data frames.
#' @param hosted_taxon_col Character. Column name indicating the hosted taxon in
#' `expanded_interaction_df` and `known_interactions_df`.
#' @param subsample_coords Integer or NULL. Number of geographic records to
#' sample per host taxon to reduce computational load. Defaults to `10000`. Set
#' to NULL to use all coordinates.
#' @param measure Character. Distance calculation method passed to
#' `geodist::geodist`. Options: `"cheap"`, `"haversine"`, `"vincenty"`, or
#' `"geodesic"`. Defaults to `"cheap"`.
#' @param metric_type Character vector specifying which distance summaries and
#' transformations to compute. Options: `"mean"`, `"median"`, `"norm"`, `"log"`.
#' `"mean"` and/or `"median"` control which base metrics are calculated, while
#' `"norm"` and `"log"` apply normalisation and log-transformation (via
#' [transform_distances()]) to those metrics. Defaults to
#' `c("mean", "median", "norm", "log")`. At least one of `"mean"` or `"median"`
#' must be included.
#' @param dnorm_args Named numeric vector giving the `mean` and `sd` parameters
#' used by transform_distances()] for normalisation. Defaults to
#' `c(mean = 0, sd = 50)`.
#' @param dnorm_from Numeric vector of two values giving the range (in km) used
#' by [transform_distances()] to anchor the scaling of the normalised scores.
#' Defaults to `c(0, 1e5)`.
#' @param round_val Numeric. Integer indicating the number of decimal places
#' used for rounding results. Detaults to `2`.
#' @param max_dist_km Numeric. Maximum distance (in km) to assign when no other
#' hosts exist for a hosted species other than the curent taxon. Defaults to
#' `30000`.
#' @param verbose Logical. Set to false to ignore warnings and messages.
#' Defaults to `TRUE`.
#' @param quiet Logical. Set to true to ignore `geodist::geodist` (recommended
#' if `measure = "cheap"`.). Defaults to `TRUE`.
#'
#' @return A data frame identical to `expanded_interaction_df` with additional
#' columns depending on the options specified in `metric_type`:
#' - `geo_min_dist_mean`: Mean of the minimum distances (km), included if
#' `"mean"` is selected.
#' - `geo_min_dist_median`: Median of the minimum distances (km), included if
#' `"median"` is selected.
#' - `geo_min_dist_mean_norm`, `geo_min_dist_median_norm`: Normalised versions,
#' included if `"norm"` (and `mean` or `median`, respectively) is selected.
#' - `geo_min_dist_mean_log`, `geo_min_dist_median_log`: Log-transformed
#' versions, included if `"log"` (and `mean` or `median`, respectively) is
#' selected.
#'
#' @details
#' - If a host or its associated hosts lack coordinate records, `NA` is returned
#'   for that interaction's distance metrics.
#' - Precomputation of distances for all host taxa improves performance for
#'   large datasets.
#' - Warnings are generated for taxa with missing geographic information or no
#'   other hosts.
#'
#' @examples
#' # Generate example dataset
#' host_spp <- paste0("Quercus_", LETTERS[1:5])
#' hosted_spp <- paste0("Agrilus_", LETTERS[1:3])
#'
#' interaction_df <- expand.grid(plant_sp = host_spp,
#'                               agrilus_sp = hosted_spp)
#'
#' interaction_df$interaction <- c(1, 0, 0, 1, 0,
#'                                 1, 0, 1, 0, 0,
#'                                 0, 0, 0, 0, 1)
#'
#' coords_df <- do.call(rbind, lapply(host_spp, function(sp) {
#'   n_ind <- sample(3:5, 1)
#'   data.frame(plant_sp = sp,
#'              lon = runif(n_ind, -10, 10),
#'              lat = runif(n_ind, 45, 55))
#' }))
#'
#' hosts_df <- subset(interaction_df,
#'                    interaction == 1,
#'                    select = c("plant_sp", "agrilus_sp"))
#'
#' # Calculate distance metrics for an example dataset
#' results <- generate_geo_metrics_df(expanded_interaction_df = interaction_df,
#'                                     coords_df = coords_df,
#'                                     known_interactions_df = hosts_df,
#'                                     host_taxon_col = "plant_sp",
#'                                     hosted_taxon_col = "agrilus_sp")
#' head(results)
#'
#' @import dplyr
#' @export
generate_geo_metrics_df <- function(expanded_interaction_df,
                                    coords_df,
                                    known_interactions_df,
                                    host_taxon_col,
                                    hosted_taxon_col,
                                    subsample_coords = 10000,
                                    measure = "cheap",
                                    metric_type = c("mean", "median",
                                                    "norm", "log"),
                                    dnorm_args = c(mean = 0, sd = 50),
                                    dnorm_from = c(min = 0, max = 1e5),
                                    round_val = 2,
                                    max_dist_km = 30000,
                                    verbose = TRUE,
                                    geodist_quiet = TRUE) {
  # Checks
  if (!is.data.frame(expanded_interaction_df)
      || !is.data.frame(known_interactions_df)
      || !is.data.frame(coords_df)) {
    stop(paste("'expanded_interaction_df', coords_df,",
               "and 'known_interactions_df' must be data frames"))
  }
  if (!(host_taxon_col %in% colnames(expanded_interaction_df))
      || !(host_taxon_col %in% colnames(known_interactions_df))
      || !(host_taxon_col %in% colnames(coords_df))) {
    stop(paste("'host_taxon_col' must indicate the host taxon column name in",
               "all three data frames"))
  }
  if (!(hosted_taxon_col %in% colnames(expanded_interaction_df))
      || !(hosted_taxon_col %in% colnames(known_interactions_df))) {
    stop(paste("'hosted_taxon_col' must indicate the host taxon column name",
               "in both 'expanded_interaction_df' and 'known_interactions_df'"))
  }
  valid_metrics <- c("mean", "median", "norm", "log")
  if (!is.character(metric_type) || !all(metric_type %in% valid_metrics)) {
    stop(paste(
      "'metric_type' must be a character vector with any of:",
      paste(valid_metrics, collapse = ", ")
    ))
  }
  if (!any(c("mean", "median") %in% metric_type)) {
    stop("'metric_type' must include at least one of 'mean' or 'median'")
  }
  if (!is.numeric(dnorm_args) || !all(c("mean", "sd") %in% names(dnorm_args))) {
    stop(paste("'dnorm_args' must be a named numeric vector with names",
               "'mean' and 'sd'"))
  }
  if (dnorm_args["sd"] <= 0) {
    stop("dnorm_args['sd'] must be positive")
  }
  if (!is.numeric(dnorm_from)) {
    stop(paste("'dnorm_from' must be a named numeric vector with names",
               "'min' and 'max'"))
  }
  if (dnorm_from[1] >= dnorm_from[2]) {
    stop("dnorm_from[1] must be less than `dnorm_from[2s]`")
  }
  if (!is.numeric(max_dist_km) || max_dist_km < 0) {
    stop("'max_dist_km' must be a non-negative numeric value")
  }
  if (!is.null(subsample_coords)) {
    if (!is.numeric(subsample_coords)
        || subsample_coords < 0
        || subsample_coords != as.integer(subsample_coords)) {
      stop("'subsample_coords' must be a non-negative integer or NULL")
    }
  }

  # Initialise results df
  results <- expanded_interaction_df
  if ("mean" %in% metric_type) results$geo_min_dist_mean <- NA
  if ("median" %in% metric_type) results$geo_min_dist_median <- NA

  # subsample geo data frame
  if (!is.null(subsample_coords)) {
    coords_df <- coords_df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(host_taxon_col))) |>
      dplyr::slice_sample(n = subsample_coords) |>
      dplyr::ungroup() |>
      as.data.frame()
  }

  # Extract species in expanded_interaction_df that belong to host taxon
  taxa <- as.character(unique(expanded_interaction_df[[host_taxon_col]]))

  # Extract only host species
  hosts <- as.character(unique(known_interactions_df[[host_taxon_col]]))

  # Precompute distances
  dist_lookup <- precompute_geo_distances(taxa = taxa,
                                          hosts = hosts,
                                          coords_df = coords_df,
                                          host_taxon_col = host_taxon_col,
                                          measure = "cheap",
                                          geodist_quiet = geodist_quiet,
                                          verbose = verbose)

  if (verbose) message("Calculating metrics per interaction...")
  for (taxon in taxa) {
    # Get index for given taxon in expanded_interaction_df
    idxs <- which(expanded_interaction_df[[host_taxon_col]] == taxon)

    for (i in idxs) {
      # Extract hosted species and its hosts
      hosted_sp <- expanded_interaction_df[[hosted_taxon_col]][i]

      host_spp <- subset(known_interactions_df,
                         get(hosted_taxon_col) == hosted_sp
                         & get(host_taxon_col) != taxon,
                         select = host_taxon_col)
      host_spp <- as.character(host_spp[[1]])

      if (length(host_spp) == 0) {
        # No hosts other than the current taxon, so set distance to max
        if (verbose) message(taxon, ": no other hosts of ", hosted_sp,
                             "; setting distance to ", max_dist_km)
        if ("mean" %in% metric_type) results$geo_min_dist_mean[i] <- max_dist_km
        if ("median" %in% metric_type) results$geo_min_dist_median[i] <- max_dist_km
      } else {
        # Extract all minimum distances for given taxon
        taxon_dists <- dist_lookup[[taxon]]

        if (is.null(taxon_dists)) {
          min_dists <- NA
        } else {
          # Extract min distances of each individual record of the given taxon
          # to each given host
          min_dists <- sapply(host_spp, function(h) {
            if (!is.null(taxon_dists[[h]])) {
              taxon_dists[[h]]
            } else {
              NA
            }
          })
          min_dists <- as.matrix(min_dists)

          # Get min distance for each individual record of the given taxon
          # across hosts
          if (ncol(min_dists) > 1) {
            min_dists <- apply(min_dists, 1, min, na.rm = TRUE)
          }
        }

        # Compute metrics
        if ("mean" %in% metric_type) {
          results$geo_min_dist_mean[i] <- round(mean(min_dists, na.rm = TRUE),
                                                round_val)
          if (is.nan(results$geo_min_dist_mean[i])) results$geo_min_dist_mean[i] <- NA
        }
        if ("median" %in% metric_type) {
          results$geo_min_dist_median[i] <- round(median(min_dists, na.rm = TRUE),
                                                  round_val)
        }
      }
    }
  }
  # Transform distances (norm, log)
  results <- transform_distances(dist_df = results,
                                 metric_type = metric_type,
                                 dnorm_args = dnorm_args,
                                 dnorm_from = dnorm_from,
                                 round_val = round_val)

  # Return updated interaction_df
  results
}


#' Precompute minimum geographic distances between taxa and host species
#'
#' This internal function computes and stores the minimum distances (in km)
#' between all occurrences of each taxon in `taxa` and all host species in
#' `hosts`. Distances are calculated using the specified method from the
#' `geodist` package. The results are stored in a nested list for fast lookup
#' in downstream calculations.
#'
#' @param taxa Character vector of taxa for which distances will be calculated.
#' @param hosts Character vector of host taxa to compare against each taxon.
#' @param coords_df Data frame of georeferenced occurrences. Must include a
#' column for host taxa (`host_taxon_col`) and numeric columns `lon` and `lat`.
#' @param host_taxon_col Character. Column name in `coords_df` that identifies
#' the host taxon.
#' @param measure Character. Distance calculation method passed to
#' `geodist::geodist`. Options: `"cheap"`, `"haversine"`, `"vincenty"`,
#' `"geodesic"`. Defaults to `"cheap"`.
#' @param geodist_quiet Logical. If `TRUE`, suppress messages from `geodist`.
#' Defaults to `TRUE`.
#'
#' @return A nested list of the form `distance_lookup[[taxon]][[host]]`, where
#' each element is a numeric vector giving the minimum distance (km) from each
#' occurrence of `taxon` to all occurrences of `host`.
#'
#' @details
#' - If a taxon or host lacks coordinate records, it is skipped with a warning.
#' - Distances are calculated individually for each occurrence of a taxon to all
#'   occurrences of a host, and only the minimum distance per taxon occurrence
#'   is stored.
#'
#' @import geodist
#' @keywords internal
precompute_geo_distances <- function(taxa,
                                     hosts,
                                     coords_df,
                                     host_taxon_col,
                                     measure = "cheap",
                                     verbose = TRUE,
                                     geodist_quiet = TRUE) {
  # Checks
  if (!is.data.frame(coords_df)) {
    stop("'coords_df' must be a data frame")
  }
  if (!is.character(taxa) || !is.character(hosts)) {
    stop("'taxa' and 'hosts' must be a character object")
  }
  if (!(host_taxon_col %in% colnames(coords_df))) {
    stop("'host_taxon_col' must indicate host taxon column name in 'coords_df'")
  }
  if (!(measure %in% c("haversine", "vincenty", "geodesic", "cheap"))) {
    stop("'measure' must be 'haversine', 'vincenty', 'geodesic', or 'cheap'")
  }
  if (!is.logical(geodist_quiet)) {
    stop("quiet must be logical")
  }

  distance_lookup <- list()

  if (verbose) message(paste(Sys.time(), "Precomputing distances..."))
  n <- 0
  for (taxon in taxa) {
    n <- n + 1

    taxon_coords <- subset(coords_df, get(host_taxon_col) == taxon)

    if (verbose) message(paste0("Taxon ", n, "/", length(taxa), ": ", taxon,
                                "; no. coords = ", nrow(taxon_coords)))

    # If no geo records, skip
    if (nrow(taxon_coords) == 0) {
      if (verbose) warning("No geo records for ", taxon, "; skipping")
      next
    }

    # Distance from this plant to all host species
    dist_list <- list()

    for (host in hosts) {
      host_coords <- subset(coords_df, get(host_taxon_col) == host)

      # If no geo records, skip
      if (nrow(host_coords) == 0) {
        if (verbose) warning("No geo records for ", host, "; skipping")
        next
      }

      # Compute the min distance to a host for each individual in taxon df
      dists_km <- geodist::geodist(taxon_coords[, c("lon", "lat")],
                                   host_coords[, c("lon", "lat")],
                                   measure = measure,
                                   quiet = geodist_quiet) / 1000
      min_dists <- apply(dists_km, 1, min)

      dist_list[[host]] <- min_dists
    }
    distance_lookup[[taxon]] <- dist_list
  }

  message(paste(Sys.time(), "Finished precomputing distances!"))
  distance_lookup
}


#' Transform distance metrics (normalisation and log-scaling)
#'
#' Internal helper function that applies optional transformations to distance
#' metrics (mean and/or median) produced by [generate_geo_metrics_df()]. It can
#' generate normalised and/or log-transformed versions of distance columns,
#' depending on the options specified in `metric_type`. Normalisation parameters
#' can be adjusted using `dnorm_args` and `dnorm_from` to control the shape and
#' bounds of the weighting curve.
#'
#' @param dist_df Data frame containing at least one of the columns
#' `geo_min_dist_mean` and/or `geo_min_dist_median`, as generated by
#' [generate_geo_metrics_df()].
#' @param metric_type Character vector specifying which distance metrics and
#' transformations to apply. Supports `"mean"`, `"median"`, `"norm"`, and
#' `"log"`.
#' - `"norm"` applies a Gaussian weighting (`dnorm`, mean = 0, SD = 50) so that
#'   nearby points have higher scores and distant ones much lower, then rescales
#'   using fixed `dnorm(0)` and `dnorm(1e5)` limits for comparability.
#' - `"log"` applies a log-transform (with a 0.1 offset) and rescales to a
#'   0 – 1 range.
#' @param dnorm_args Named numeric vector giving the `mean` and `sd` parameters
#'   for the Gaussian weighting used in normalisation. Defaults to
#'   `c(mean = 0, sd = 50)`.
#' @param dnorm_from Numeric vector of two values giving the range (in km)
#'   over which to anchor the scaling of the normalised scores. Defaults to
#'   `c(0, 1e5)`.
#' @param round_val Numeric. Integer indicating the number of decimal places
#' used for rounding results. Detaults to `2`.
#'
#' @return
#' A data frame identical to `dist_df` with additional columns corresponding to
#' the transformations requested:
#'   - `geo_min_dist_mean_norm` if `"mean"` and `"norm"` are included.
#'   - `geo_min_dist_median_norm` if `"median"` `"norm"` are included.
#'   - `geo_min_dist_mean_log` if `"mean"` and `"log"` are included.
#'   - `geo_min_dist_median_log` if `"median"` and `"log"` are included.
#'
#' @importFrom scales rescale
#' @keywords internal
transform_distances <- function(dist_df,
                                metric_type,
                                dnorm_args = c(mean = 0, sd = 50),
                                dnorm_from = c(min = 0, max = 1e5),
                                round_val = 2) {
  # Checks
  if (!is.data.frame(dist_df)) {
    stop("'dist_df' must be data frames")
  }
  if (!any(c("geo_min_dist_mean", "geo_min_dist_median") %in% colnames(dist_df))) {
    stop(paste("'dist_df' must contain at least one of the columns",
               "'geo_min_dist_mean' or 'geo_min_dist_median'"))
  }
  valid_metrics <- c("mean", "median", "norm", "log")
  if (!is.character(metric_type) || !all(metric_type %in% valid_metrics)) {
    stop(paste0(
      "'metric_type' must be a character vector with any of: ",
      paste(valid_metrics, collapse = ", ")
    ))
  }
  if (!any(c("mean", "median") %in% metric_type)) {
    stop("'metric_type' must include at least one of 'mean' or 'median'")
  }
  if (!is.numeric(dnorm_args) || !all(c("mean", "sd") %in% names(dnorm_args))) {
    stop(paste("'dnorm_args' must be a named numeric vector with names",
               "'mean' and 'sd'"))
  }
  if (dnorm_args["sd"] <= 0) {
    stop("dnorm_args['sd'] must be positive")
  }
  if (!is.numeric(dnorm_from)) {
    stop(paste("'dnorm_from' must be a named numeric vector with names",
               "'min' and 'max'"))
  }
  if (dnorm_from[1] >= dnorm_from[2]) {
    stop("dnorm_from[1] must be less than dnorm_from[2s]")
  }
  dist_tf <- dist_df

  # Weight distances using a normal curve (e.g.: mean = 0, SD = 50) so nearby
  # points get higher scores and distant ones (e.g.: >100 km) much lower.
  # Rescale using dnorm_from limits to keep results comparable across datasets.
  if ("norm" %in% metric_type) {
    mu <- dnorm_args[["mean"]]
    sigma <- dnorm_args[["sd"]]
    from_vals <- range(-dnorm(dnorm_from, mean = mu, sd = sigma))

    if ("mean" %in% metric_type) {
      min_dist_mean_norm <- scales::rescale(-dnorm(dist_tf$geo_min_dist_mean,
                                                   mean = mu, sd = sigma),
                                            from = from_vals)
      min_dist_mean_norm <- round(min_dist_mean_norm, round_val)
      dist_tf$geo_min_dist_mean_norm <- min_dist_mean_norm
    }
    if ("median" %in% metric_type) {
      min_dist_median_norm <- scales::rescale(-dnorm(dist_tf$geo_min_dist_median,
                                                     mean = mu, sd = sigma),
                                              from = from_vals)
      min_dist_median_norm <- round(min_dist_median_norm, round_val)
      dist_tf$geo_min_dist_median_norm <- min_dist_median_norm
    }
  }
  # Log-transform distances (with small offset to deal with 0s)
  if ("log" %in% metric_type) {
    if ("mean" %in% metric_type) {
      min_dist_mean_log <- scales::rescale(log(dist_tf$geo_min_dist_mean + 0.1))
      min_dist_mean_log <- round(min_dist_mean_log, round_val)
      dist_tf$geo_min_dist_mean_log <- min_dist_mean_log
    }
    if ("median" %in% metric_type) {
      min_dist_median_log <- scales::rescale(log(dist_tf$geo_min_dist_median + 0.1))
      min_dist_median_log <- round(min_dist_median_log, round_val)
      dist_tf$geo_min_dist_median_log <- min_dist_median_log
    }
  }
  dist_tf
}


#' Calculate mean and minimum phylogenetic distance metrics
#'
#' This function calculates mean and minimum phylogenetic distances between each
#' focal host-taxon species and the known hosts of its associated hosted
#' species. The phylogenetic distances are derived from the inverse of a
#' covariance matrix computed from a phylogeny. Distances are rescaled to the
#' range [0, 1].
#'
#' @param expanded_interaction_df A data frame containing all possible
#' host-taxon species – hosted combinations for which phylogenetic metrics
#' should be calculated.
#' @param known_interactions_df A data frame of known interactions, containing
#' at least the columns specified by `host_taxon_col` and `hosted_taxon_col`.
#' Used to identify known hosts for each hosted species.
#' @param phylo A phylogenetic tree of class \code{"phylo"} (from the \pkg{ape}
#' package). The tree should include all host taxa in its tip labels.
#' @param host_taxon_col A character string giving the name of the column in
#' both input data frames that contains the host taxon.
#' @param hosted_taxon_col A character string giving the name of the column in
#' both input data frames that contains the hosted taxon.
#'
#' @return A data frame identical to `expanded_interaction_df`, with two new
#'   columns:
#'   \describe{
#'     \item{phylo_dist_mean}{Mean phylogenetic distance to known hosts.}
#'     \item{phylo_dist_min}{Minimum phylogenetic distance to known hosts.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Example data
#' phy <- ape::rtree(5, tip.label = c("A", "B", "C", "D", "E"))
#' 
#' known_df <- data.frame(
#'   host = c("A", "B", "C"),
#'   hosted = c("X", "X", "Y")
#' )
#' expanded_df <- expand.grid(host = c("A", "B", "C", "D", "E"),
#'                            hosted = c("X", "Y"))
#'
#' # Compute phylo metrics
#' generate_phylo_metrics_df(expanded_df, known_df, phy,
#'                           host_taxon_col = "host",
#'                           hosted_taxon_col = "hosted")
#' }
#'
#' @importFrom ape vcv.phylo
#' @importFrom scales rescale
#' @export
generate_phylo_metrics_df <- function(expanded_interaction_df,
                                      known_interactions_df,
                                      tree,
                                      host_taxon_col,
                                      hosted_taxon_col) {
  # Checks
  if (!is.data.frame(expanded_interaction_df)
      || !is.data.frame(known_interactions_df)) {
    stop(paste("'expanded_interaction_df', and",
               "'known_interactions_df' must be data frames"))
  }
  if (!(host_taxon_col %in% colnames(expanded_interaction_df))
      || !(host_taxon_col %in% colnames(known_interactions_df))) {
    stop(paste("'host_taxon_col' must indicate the host taxon column name in",
               "both data frames"))
  }
  if (!(hosted_taxon_col %in% colnames(expanded_interaction_df))
      || !(hosted_taxon_col %in% colnames(known_interactions_df))) {
    stop(paste("'hosted_taxon_col' must indicate the host taxon column name",
               "in both 'expanded_interaction_df' and 'known_interactions_df'"))
  }
  if (!inherits(tree, "phylo") || !("tip.label" %in% names(tree))) {
    stop("Input tree must be of class 'phylo' and contain a 'tip.label' field")
  }

  # Create a phylo cov distance matrix (Hadfield & Nakagawa, 2010)
  phylo_cov <- ape::vcv.phylo(tree)
  phylo_cov <- max(phylo_cov) - phylo_cov

  results <- expanded_interaction_df
  results$phylo_dist_mean <- NA
  results$phylo_dist_min <- NA

  # For every row in the df, calculate the mean phylogenetic distance of
  # the focal species (taxon) to the hosts of the hosted species
  for (i in seq_len(nrow(expanded_interaction_df))) {
    # Extract current host-taxon species and hosted species
    taxon <- expanded_interaction_df[[host_taxon_col]][i]
    hosted_sp <- expanded_interaction_df[[hosted_taxon_col]][i]

    # Extract hosts of the given hosted species (excluded focus_sp if host)
    host_spp <- subset(known_interactions_df,
                       get(hosted_taxon_col) == hosted_sp
                       & get(host_taxon_col) != taxon,
                       select = host_taxon_col)
    host_spp <- as.character(host_spp[[1]])

    # If the current taxon is the only known host, set the distance to max
    if (length(host_spp) == 0) {
      mean_dist <- max(phylo_cov)
      min_dist <- max(phylo_cov)
    } else {
      # Calculate mean/min phylo distance for host-taxon species to the hosts
      # of the given hosted species
      cov_rows <- which(rownames(phylo_cov) == taxon)
      cov_cols <- which(colnames(phylo_cov) %in% host_spp)
      if (length(cov_rows) != 0 && length(cov_cols) != 0) {
        mean_dist <- mean(phylo_cov[cov_rows, cov_cols])
        min_dist <- min(phylo_cov[cov_rows, cov_cols])
      }
    }
    results$phylo_dist_mean[i] <- mean_dist
    results$phylo_dist_min[i] <- min_dist
  }
  # Normalise distances
  results$phylo_dist_mean <- scales::rescale(results$phylo_dist_mean)
  results$phylo_dist_min <- scales::rescale(results$phylo_dist_min)
  results
}
