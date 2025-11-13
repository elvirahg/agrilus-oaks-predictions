#' Standardise the tip labels of a phylogeny
#'
#' Applies standardisation steps to a \code{phylo} tree's tip labels, including:
#' \itemize{
#'   \item Cleaning tip label names (using regex).
#'   \item Renaming specified tip labels.
#'   \item Removing duplicate tip labels, keeping only one instance.
#' }
#' It returns a phylogeny.
#'
#' @param tree A phylogeny of class \code{phylo}.
#' @param clean_labels Logical; if \code{TRUE}, clean tip labels (for plotting purposes).
#' @param pattern Character; regex pattern for extracting species-level labels.
#' @param replacement Character; replacement pattern for use in \code{gsub}.
#' @param remove_duplicates Logical; if \code{TRUE}, remove duplicated tip labels.
#' @param old_labels Character vector of labels to rename (passed to internal \code{rename_labels}).
#' @param new_labels Character vector of new labels to use (passed to internal \code{rename_labels}).
#' @return A \code{phylo} tree with standardised tip labels.
#'
#' @examples
#' \dontrun{
#' tree <- ape::read.tree("path/to/tree.tre")
#' tree <- standardise_phylo(tree,
#'                           clean_labels = TRUE,
#'                           pattern = "^([A-Z])[a-z]+_*([a-z]+).*",
#'                           replacement = "\\1. \\2",
#'                           remove_duplicates = TRUE,
#'                           old_labels = c"Quercus maxima",
#'                           new_labels = "Quercus rubra"))
#' plot(tree)
#' }
#'
#' @export
standardise_phylo <- function(tree,
                              clean_labels = TRUE,
                              pattern = "^([A-Z])[a-z]+_[×|x]*_*([a-z]+).*",
                              replacement = "\\1. \\2",
                              remove_duplicates = TRUE,
                              old_labels = NULL,
                              new_labels = NULL) {
  # Check input
  if (!inherits(tree, "phylo") || !("tip.label" %in% names(tree))) {
    stop("Input tree must be of class 'phylo' and contain a 'tip.label' field")
  }
  if (!is.logical(clean_labels) || length(clean_labels) != 1) {
    stop("'clean_labels' must be logical")
  }
  if (!is.character(pattern) || length(pattern) != 1) {
    stop("'pattern' must be a single string (regex)")
  }
  if (!is.character(replacement) || length(replacement) != 1) {
    stop("'replacement' must be a single string")
  }
  if (!is.logical(remove_duplicates) || length(remove_duplicates) != 1) {
    stop("'remove_duplicates' must be logical")
  }

  # Clean tip labels (for plotting purposes)
  if (clean_labels) {
    tree$tip.label <- gsub(pattern = pattern,
                           replacement = replacement,
                           x = tree$tip.label)
  }

  # Rename labels
  if (!is.null(old_labels)) {
    tree <- rename_labels(tree, old_labels, new_labels)
  }

  # Remove duplicates (keep one randomly)
  if (remove_duplicates) {
    tree <- remove_duplicate_tips(tree)
  }

  tree
}


#' Rename tip labels in a phylogeny
#'
#' Replaces specified tip labels in a \code{phylo} tree with new labels.
#' Both \code{old_labels} and \code{new_labels} must be character vectors
#' of the same length.
#'
#' @param tree Phylogeny of class \code{phylo} containing \code{tip.label}.
#' @param old_labels Character vector of tip labels to be replaced.
#' @param new_labels Character vector of new tip labels to replace the old ones.
#'
#' @return A \code{phylo} tree with specified tip labels renamed.
#'
#' @keywords internal
rename_labels <- function(tree, old_labels, new_labels) {

  if (!is.character(old_labels)) {
    stop("'old_labels' must be a character vector")
  }
  if (!is.character(new_labels)) {
    stop("'new_labels' must be a character vector")
  }
  if (length(old_labels) != length(new_labels)) {
    stop("Length of 'old_labels' and 'new_labels' must be equal")
  }

  for (i in seq_along(old_labels)) {
    old_label <- old_labels[i]
    new_label <- new_labels[i]

    idx <- which(tree$tip.label == old_label)

    if (length(idx) == 0) {
      warning(sprintf("No match found for label: '%s'", old_label))
    } else {
      tree$tip.label[idx] <- new_label
    }
  }

  tree
}


#' Remove duplicate tip labels
#'
#' Removes duplicated tip labels, keeping only one randomly selected instance.
#'
#' @param labels A character vector of tip labels.
#'
#' @return A character vector of tip labels with duplicates removed.
#'
#' @import ape
#' @keywords internal
remove_duplicate_tips <- function(tree) {
  labels <- tree$tip.label
  duplicate_names <- unique(labels[duplicated(labels)])

  idxs_to_remove <- unlist(
    lapply(duplicate_names, function(name) {
      idx <- which(labels == name)
      sample(idx, length(idx) - 1)
    })
  )

  if (length(idxs_to_remove) > 0) {
    tree <- ape::drop.tip(tree, idxs_to_remove)
  }

  tree
}


#' Generate host presence/absence data frame for `caper::comparative.data()`
#'
#' This function generates a presence/absence data frame for host species.
#' It takes a character vector of observed hosts and a full list of all
#' species, and returns a data frame marking species as hosts (1) or non-hosts
#' (0). Can be used to prepare data for input to `caper::comparative.data()`.
#'
#' @param observed_hosts A character vector of observed host species.
#' @param species_list A character vector of all species of interest. Each
#'   species in this list will appear in the output, with host.status 1 if
#'   present in `observed_hosts` and 0 otherwise.
#'
#' @return A data.frame with two columns:
#'   \describe{
#'     \item{species}{Species names from `species_list`.}
 #'    \item{host.status}{Integer (0/1) indicating whether the species is a known host.}
#'   }
#'
#' @examples
#' oak_hosts <- c("Quercus_robur", "Quercus_petraea")
#' all_oaks <- c("Quercus_rubra", "Quercus_robur",
#'               "Quercus_alba", "Quercus_petraea")
#' generate_pres_abs_df(oak_hosts, all_oaks)
#' #   species host.status
#' # 2 Quercus_rubra 0
#' # 3 Quercus_robur 1
#' # 4 Quercus_petraea 1
#' # 5 Quercus_alba 0
#'
#' @export
generate_pres_abs_df <- function(observed_hosts, species_list) {
  # Check input
  if (!is.character(observed_hosts)) {
    stop("'observed_hosts' must be a character vector")
  }
  if (!is.character(species_list)) {
    stop("'species_list' must be a character vector")
  }

  # Generate host status presence/absence data frame
  host_status <- data.frame(species = species_list,
                            host.status = 0)
  host_status$host.status[host_status[, 1] %in% observed_hosts] <- 1

  host_status
}


#' Retrieve GBIF taxon keys for a list of species
#'
#' This function is a lightweight wrapper around `rgbif::name_lookup()`. It
#' takes a character vector of species names and queries GBIF to retrieve their
#' GBIF taxon keys, and returns a data frame indicating the GBIF key for each
#' species. Some species may not be found and will require manual lookup.
#' Function allows manual overriding of any parameters in
#' `name_lookup()` via `...`.
#'
#' @param species Character vector of species names to look up.
#' @param higher_taxon Character string specifying a higher taxon key to
#' restrict the search.
#' @param rank Character string specifying the taxonomic rank to query
#' (default "SPECIES").
#' @param ... Additional arguments passed directly to `rgbif::name_lookup()`.
#'
#' @return A data.frame with two columns:
#'   \describe{
#'     \item{species}{Original species name.}
#'     \item{key}{Corresponding GBIF taxon key (NA if not found).}
#'   }
#'
#' @examples
#' # Simple example using dummy species
#' \dontrun{
#' species_list <- c("Quercus robur", "Quercus alba")
#' get_gbif_keys(species_list, higher_taxon = "220")
#' # Override `rgbif::name_lookup` `status` parameter
#' get_gbif_keys(species_list, higher_taxon = "220", status = "Synonym")
#' }
#'
#' @import rgbif
#' @export
get_gbif_keys <- function(species,
                          higher_taxon = "",
                          rank = "SPECIES",
                          ...) {
  # Check input
  if (!is.character(species)) {
    stop("'species' must be a character vector")
  }
  if (!is.character(higher_taxon)) {
    stop("'higher_taxon' must be a character (no default)")
  }
  if (!is.character(rank)) {
    stop("'rank' must be a character")
  }
  # Look up GBIF keys, forwarding any additional parameters to name_lookup()
  keys_info <- lapply(species, function(sp) {
    rgbif::name_lookup(query = sp,
                       rank = rank,
                       higherTaxonKey = higher_taxon,
                       status = "Accepted",
                       limit = 1,
                       ...)
  })

  # Fill keys_df with valid keys
  keys_df <- data.frame(species = species, key = NA)
  for (i in seq_along(keys_info)) {
    sp_info <- keys_info[[i]]
    if (!is.null(sp_info$hierarchies)) {
      sp <- sp_info$data$species
      if (sp %in% species) {
        keys_df$key[i] <- sp_info$data$key
      }
    }
  }

  keys_df
}


#' Clean Scientific Names
#'
#' This function takes a character vector of un-curated scientific names and
#' performs a series of cleaning steps to produce a standardized version
#' suitable for comparison or filtering. Built to clean `verbatimScientificName`
#' fields from an `rgbif` dataset.
#'
#' The cleaning includes:
#' - Removing abbreviations and everything after a capital letter, punctuation,
#'   or bracket (e.g., "Marsh.")
#' - Removing uncertainty markers, rank indicators, and common qualifiers
#'   (e.g., "subsp.")
#' - Removing hybrid markers ("x", "×", "X", "hybrid")
#' - Removing miscellaneous symbols (e.g., "_")
#' - Collapsing multiple spaces into a single space, replacing non-breaking
#'   spaces with regular spaces, and trimming leading/trailing whitespace
#'
#' @note
#' This approach is conservative: it removes common sources of noise and
#' non-standard elements in scientific names, but it is not guaranteed
#' to produce fully standardised or globally accurate names. Some valid
#' variations or subspecies may be altered or truncated.
#'
#' @param x A character vector of scientific names to be cleaned.
#'
#' @return A character vector of cleaned scientific names.
#'
#' @examples
#' names <- c("Quercus martensiana f. perplexans", "Quercus × robur L.")
#' clean_sci_name(names)
#'
#' @export
clean_taxon_name <- function(x) {
  if (!is.character(x)) {
    stop("'x' must be a character vector")
  }

  x |>
    # Remove non-breaking white-spaces
    gsub(" *\u00A0", " ", x = _) |>
    # remove everything after a capital letter, punctuation, or braket
    gsub("([a-z ])[A-Z].*|['‘,\\(\\[].*", "\\1", x = _) |>
    # remove abbreviations
    gsub(" [a-z]+\\.", " ", x = _) |>
    # remove uncertainty markers, rank indicators, and common qualifiers
    gsub("\\b(-|subsp|var|cf|f|forma|de|ex|vera|sensu)\\b\\.?", " ",
         x = _, ignore.case = TRUE) |>
    # remove hybrid markers
    gsub("\\b[xX]\\b|\\bhybrid\\b|×", "", x = _) |>
    # remove symbols
    gsub("[<>\\?\\&_\\.ß]|\\b-\\b", "", x = _) |>
    # collapse multiple spaces
    gsub("\\ +", " ", x = _) |>
    # Remove whitespaces
    trimws()
}


#' Filter GBIF occurrence data based on common quality criteria
#'
#' This function filters GBIF occurrence records retrieved via `rgbif`. All
#' filters are optional and have default settings.
#'
#' @param data A `data.frame` or `tibble` of GBIF occurrence records.
#'   Must contain the following columns depending on enabled filters:
#'   `decimalLatitude`, `decimalLongitude`, `coordinateUncertaintyInMeters`,
#'   `individualCount`, `occurrenceStatus`, `year`, and `issue`.
#' @param remove_no_coords Logical; if `TRUE`, removes records with missing
#'   latitude or longitude values. Defaults to `TRUE`.
#' @param coord_uncertainty_thr Numeric or `NULL`; maximum accepted value for
#'   `coordinateUncertaintyInMeters`. Records with uncertainty greater than this
#'   threshold are removed. If `NULL`, this filter is skipped. Defaults to
#'   `5000`.
#' @param remove_zero_indiv_count Logical; if `TRUE`, removes records with
#'   `individualCount` equal to 0. Defaults to `TRUE`.
#' @param max_indiv_count Numeric or `NULL`; removes records where
#'   `individualCount` exceeds this value. If `NULL`, this filter is skipped.
#'   Defaults to `10`.
#' @param remove_absent Logical; if `TRUE`, keeps only records with
#'   `occurrenceStatus == "PRESENT"`. Defaults to `TRUE`.
#' @param min_year Numeric or `NULL`; removes records with `year` earlier than
#'   this threshold. If `NULL`, this filter is skipped. Defaults to `1950`.
#' @param issues Character vector or `NULL`; list of GBIF issue codes to remove
#'   If `NULL`, this filter is skipped. Defaults to a set of common issues.
#'
#' @details
#' The function performs filtering operations in the following order:
#' 1. Remove records with missing coordinates (if enabled).
#' 2. Filter by coordinate uncertainty (if enabled).
#' 3. Filter by maximum and zero individual count (if enabled).
#' 4. Remove records where `occurrenceStatus != "PRESENT"` (if enabled).
#' 5. Filter by minimum collection year (if enabled).
#' 6. Remove records with specified GBIF data quality issues (if enabled).
#'
#' @return A filtered `data.frame` or `tibble` of occurrence records.
#'
#' @examples
#' \dontrun{
#' # Filter GBIF data with default thresholds
#' filtered_geo <- filter_gbif_data(geo_data)
#'
#' # Allow higher coordinate uncertainty
#' filtered_geo <- filter_gbif_data(or `tibble`, coord_uncertainty_thr = 10000)
#'
#' @export
filter_gbif_data <- function(data,
                             remove_no_coords = TRUE,
                             coord_uncertainty_thr = 5000,
                             remove_zero_indiv_count = TRUE,
                             max_indiv_count = 10,
                             remove_absent = TRUE,
                             min_year = 1950,
                             issues = c("COORDINATE_INVALID",
                                        "COORDINATE_OUT_OF_RANGE",
                                        "COORDINATE_REPROJECTION_SUSPICIOUS",
                                        "COORDINATE_UNCERTAINTY_METERS_INVALID",
                                        "COUNTRY_COORDINATE_MISMATCH",
                                        "PRESUMED_SWAPPED_COORDINATE",
                                        "PRESUMED_NEGATED_LONGITUDE",
                                        "PRESUMED_NEGATED_LATITUDE",
                                        "ZERO_COORDINATE",
                                        "FOSSIL_SPECIMEN")) {
  if (!inherits(data, "data.frame")) {
    stop("Input 'data' must be an rgbif tibble or rgbif data frame")
  }
  if (!is.logical(remove_no_coords)) {
    stop("'remove_no_coords' must be logical")
  }
  if (!is.numeric(coord_uncertainty_thr)) {
    stop("'coord_uncertainty_thr' must be an integer")
  }
  if (!is.logical(remove_zero_indiv_count)) {
    stop("'remove_zero_indiv_count' must be logical")
  }
  if (!is.numeric(max_indiv_count)) {
    stop("'max_indiv_count' must be an integer")
  }
  if (!is.logical(remove_absent)) {
    stop("'remove_absent' must be logical")
  }
  if (!is.numeric(min_year)) {
    stop("'min_year' must be an integer")
  }
  if (!is.character(issues)) {
    stop("'issues' must be a character vector")
  }

  result <- data
  if (remove_no_coords) {
    if ((!"decimalLatitude" %in% names(result))
        || (!"decimalLongitude" %in% names(result))) {
      stop("Column 'decimalLatitude' and/or 'decimalLongitude' missing")
    }
    result <- subset(result,
                     !is.na(decimalLatitude) & !is.na(decimalLongitude))
  }
  if (!is.null(coord_uncertainty_thr)) {
    if (!"coordinateUncertaintyInMeters" %in% names(result)) {
      stop("Column 'coordinateUncertaintyInMeters' missing")
    }
    result <- subset(result,
                     is.na(coordinateUncertaintyInMeters)
                     | coordinateUncertaintyInMeters <= coord_uncertainty_thr)
  }
  if (!is.null(max_indiv_count)) {
    if (!"individualCount" %in% names(result)) {
      stop("Column 'individualCount' missing")
    }
    result <- subset(result,
                     is.na(individualCount)
                     | individualCount <= max_indiv_count)
  }
  if (remove_zero_indiv_count) {
    if (!"individualCount" %in% names(result)) {
      stop("Column 'individualCount' missing")
    }
    result <- subset(result,
                     is.na(individualCount)
                     | individualCount > 0)
  }
  if (remove_absent) {
    if (!"occurrenceStatus" %in% names(result)) {
      stop("Column 'occurrenceStatus' missing")
    }
    result <- subset(result,
                     occurrenceStatus == "PRESENT")
  }
  if (!is.null(min_year)) {
    if (!"year" %in% names(result)) {
      stop("Column 'year' missing")
    }
    result <- subset(result,
                     is.na(year)
                     | year >= min_year)
  }
  if (!is.null(issues)) {
    if (!"issue" %in% names(result)) {
      stop("Column 'issue' missing")
    }
    result <- subset(result,
                     !grepl(pattern = paste(issues, collapse = "|"), issue))
  }
  result
}


#' Add a species with a country or region centroid to a geographic dataset
#'
#' This function adds a single species record to an existing data frame of
#' geographic occurrences. It uses the `countryref` dataset from the
#' `CoordinateCleaner` package to obtain the centroid coordinates of a specified
#' country or region.
#'
#' @param df A data frame containing a species column, a latitude column and a
#' longitude column.
#' @param species_name Character. The species name to add.
#' @param iso3 Character. The ISO3 code of the country in which the species is
#'   present.
#' @param region Character or NULL (default NULL). Optional region within the
#'   country (e.g., a province/state). If provided, the function will return
#'   the centroid of the region instead of the whole country.
#' @param col_species Character. Name of the column in `df` for species names.
#' @param col_lat Character. Name of the column in `df` for latitude values.
#' @param col_lon Character. Name of the column in `df` for longitude values.
#'
#' @return A data frame with the new species row appended to the input `df`.
#'
#' @details This function was designed for when species are missing from a GBIF
#' dataset and you want to provide a placeholder occurrence at a country/region
#' centroid.
#'
#' @examples
#' # Add Quercus sagrana in Cuba
#' df <- add_species_centroid(gbif_df,
#'                            "Quercus sagrana",
#'                            "CUB")
#' @import CoordinateCleaner
#' @export
add_species_centroid <- function(df,
                                 species_name,
                                 iso3,
                                 region = NULL,
                                 col_species = "species",
                                 col_lon = "lon",
                                 col_lat = "lat") {
  if (!is.data.frame(df)) {
    stop("'df' must be a dataframe")
  }
  if (!is.character(species_name)) {
    stop("'species_name' must be a character")
  }
  if (!is.character(iso3)) {
    stop("'iso3' must be a character indicating iso3 country code")
  }
  if (!is.character(region)) {
    stop("'region' must be a character")
  }
  if (!is.character(col_species)) {
    stop("'col_species' must be a character")
  }
if (!is.character(col_lon) || !is.character(col_lat)) {
    stop("'col_lon' and 'col_lat' must be characters")
  }

  # Subset country/region
  if (is.null(region)) {
    coords <- subset(CoordinateCleaner::countryref,
                     iso3 == iso3
                     & type == "country")[1, c("centroid.lat", "centroid.lon")]
  } else {
    coords <- subset(CoordinateCleaner::countryref,
                     iso3 == iso3
                     & name == region)[1, c("centroid.lat", "centroid.lon")]
  }

  # Combine into a data frame
  new_row <- data.frame(species_name,
                        coords$centroid.lon,
                        coords$centroid.lat)
  colnames(new_row) <- c(col_species, col_lon, col_lat)

  # Append to existing dataset
  rbind(df, new_row)
}


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
#' - `geo_min_dist_mean`: Mean of the minimum distances (km), included if `"mean"`
#' is selected.
#' - `geo_min_dist_median`: Median of the minimum distances (km), included if
#' `"median"` is selected.
#' - `geo_min_dist_mean_norm`, `geo_min_dist_median_norm`: Normalised versions,
#' included if `"norm"` (and `mean` or `median`, respectively) is selected.
#' - `geo_min_dist_mean_log`, `geo_min_dist_median_log`: Log-transformed versions,
#' included if `"log"` (and `mean` or `median`, respectively) is selected.
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
#' results <- generate_dist_metrics_df(expanded_interaction_df = interaction_df,
#'                                     coords_df = coords_df,
#'                                     known_interactions_df = hosts_df,
#'                                     host_taxon_col = "plant_sp",
#'                                     hosted_taxon_col = "agrilus_sp")
#' head(results)
#'
#' @import dplyr
#' @export
generate_dist_metrics_df <- function(expanded_interaction_df,
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
#' metrics (mean and/or median) produced by [generate_dist_metrics_df()]. It can
#' generate normalised and/or log-transformed versions of distance columns,
#' depending on the options specified in `metric_type`. Normalisation parameters
#' can be adjusted using `dnorm_args` and `dnorm_from` to control the shape and
#' bounds of the weighting curve.
#'
#' @param dist_df Data frame containing at least one of the columns
#' `geo_min_dist_mean` and/or `geo_min_dist_median`, as generated by
#' [generate_dist_metrics_df()].
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


#' Plot Geographic Distance Metrics
#'
#' Generates one or more scatter plots of distance metrics from a data frame,
#' automatically sorting and transforming them based on their column names.
#' The function supports any combination of raw, normalised, or log-transformed
#' distance columns.
#' Each distance metric is plotted as a scatter plot, and all plots are
#' arranged in a grid layout using \pkg{gridExtra}.
#'
#' @param dist_data A data frame containing the distance metrics to be plotted.
#' @param dist_cols A character vector specifying the column names within
#' `dist_data` to plot. Columns should correspond to distance metrics such as
#'`geo_min_dist_mean`, `geo_min_dist_mean_norm`, or `geo_min_dist_mean_log`.
#'
#' @return
#' A composite grid of \pkg{ggplot2} scatter plots.
#'
#' @examples
#' \dontrun{
#' plot_geo_dists(
#'   interaction_data,
#'   dist_cols = c("geo_min_dist_mean", "geo_min_dist_mean_norm", "geo_min_dist_mean_log")
#' )
#' }
#'
#' @importFrom gridExtra grid.arrange
#' @export
plot_geo_dists <- function(dist_data, dist_cols) {
  # Checks
  if (!is.data.frame(dist_data)) {
    stop("'dist_data' must be a data frame")
  }
  if (!is.character(dist_cols)) {
    stop("'dist_cols' must be a character vector")
  }
  if (!all(dist_cols %in% colnames(dist_data))) {
    stop("One or more specified columns do not exist in 'dist_data'")
  }

  # Plot colours
  cols <- hcl.colors(length(dist_cols), "Harmonic")

  # Plot list
  dist_plots <- vector("list", length(dist_cols))

  for (i in seq_along(dist_cols)) {
    colname <- dist_cols[i]

    # Plot text
    title_text <- paste0("Distance (", colname, ")")
    if (colname %in% c("geo_min_dist_mean", "geo_min_dist_median")) {
      y_text <- paste0("Distance (km)")
    } else if (colname %in% c("geo_min_dist_mean_norm", "geo_min_dist_median_norm")) {
      y_text <- paste0("Normalised distance")
    } else {
      y_text <- paste0("Log distance")
    }

    # Sort and plot distances
    if (colname %in% c("geo_min_dist_mean", "geo_min_dist_median",
                       "geo_min_dist_mean_log", "geo_min_dist_median_log")) {
      dist_sorted <- sort(dist_data[[colname]], decreasing = FALSE)
      dist_plots[[i]] <- make_dist_plot(dist_sorted, cols[i],
                                        title_text, y_text)
    } else if (colname %in% c("geo_min_dist_mean_norm", "geo_min_dist_median_norm")) {
      dist_sorted <- sort(dist_data[[colname]], decreasing = TRUE)
      dist_plots[[i]] <- make_dist_plot(-dist_sorted, cols[i],
                                        title_text, y_text)
    }
  }

  # Arrange plots
  ncol <- ceiling(length(dist_plots) / 2)
  do.call(gridExtra::grid.arrange, c(dist_plots, ncol = ncol))
}


#' Create a Single Distance Plot
#'
#' Internal helper to generate a scatter plot of a single distance metric.
#'
#' @param values Numeric vector of distances to plot.
#' @param colour Character string specifying point colour.
#' @param title Character string for the plot title.
#' @param ylab_text Character string for the y-axis label.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @import ggplot2
#' @keywords internal
make_dist_plot <- function(values,
                           colour,
                           title,
                           ylab_text) {
  # Checks
  if (!is.numeric(values)) {
    stop("'values' must be a numeric vector")
  }
  if (!is.character(colour)) {
    stop("'colour' must be a character")
  }
  if (!is.character(title)) {
    stop("'title' must be a character")
  }
  if (!is.character(ylab_text)) {
    stop("'ylab_text' must be a character")
  }

  # Plot
  ggplot2::ggplot(data = data.frame(x = seq_along(values), y = values),
                  ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(colour = colour) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    ggplot2::xlab("Index") +
    ggplot2::ylab(ylab_text)
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
#' vars <- c("geo.dist",
#'           "phylo.dist",
#'           "(1 | gr(plant_sp, cov = phylo_cov))")
#'
#' models <- generate_formulas(
#'   variables = vars,
#'   response = "interaction",
#'   include_null = TRUE
#' )
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
#' formulas <- c("y ~ x + (1 | group)")
#' pattern <- "(1 | group)"
#' replacements <- c("(x | group)", "(x + z | group)")
#' add_formula_variations(formulas, pattern, replacements)
#' # [1] "y ~ x + (1 | group)"
#' # [2] "y ~ x + (x | group)"
#' # [3] "y ~ x + (x + z | group)"
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
    formula_rd <- Reduce(function(x, pattern) gsub(pattern, "random", x, perl = TRUE),
                         random_regex,
                         init = formula_strings)
    n_vars <- sapply(formula_rd, function(f) {
      length(strsplit(f,
                      " \\+ ",
                      perl = TRUE)[[1]])
    })
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
