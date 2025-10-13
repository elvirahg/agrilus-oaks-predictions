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
#'                           old_labels = c"Q. maxima",
#'                           new_labels = "Q. rubra"))
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
    stop("Input tree must be of class 'phylo' and contain a 'tip.label' field.")
  }
  if (!is.logical(clean_labels) || length(clean_labels) != 1) {
    stop("'clean_labels' must be a single boolean.")
  }
  if (!is.character(pattern) || length(pattern) != 1) {
    stop("'pattern' must be a single string (regex).")
  }
  if (!is.character(replacement) || length(replacement) != 1) {
    stop("'replacement' must be a single string.")
  }
  if (!is.logical(remove_duplicates) || length(remove_duplicates) != 1) {
    stop("'remove_duplicates' must be a single boolean.")
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
#' @param tree A phylogeny of class \code{phylo} containing \code{tip.label}.
#' @param old_labels A character vector of tip labels to be replaced.
#' @param new_labels A character vector of new tip labels to replace the old ones.
#'
#' @return A \code{phylo} tree with specified tip labels renamed.
#'
#' @keywords internal
rename_labels <- function(tree, old_labels, new_labels) {

  if (!is.character(old_labels)) {
    stop("'old_labels' must be a character vector.")
  }
  if (!is.character(new_labels)) {
    stop("'new_labels' must be a character vector.")
  }
  if (length(old_labels) != length(new_labels)) {
    stop("Length of 'old_labels' and 'new_labels' must be equal.")
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


#' Generate host presence/absence dataframe for `caper::comparative.data()`
#'
#' This function generates a presence/absence dataframe for host species.
#' It takes a character vector of observed hosts and a full list of all
#' species, and returns a dataframe marking species as hosts (1) or non-hosts
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
 #'     \item{host.status}{Integer (0/1) indicating whether the species is a known host.}
#'   }
#'
#' @examples
#' oak_hosts <- c("Quercus_robur", "Quercus_petraea")
#' all_oaks <- c("Quercus_rubra", "Quercus_robur", "Quercus_alba", "Quercus_petraea", ...)
#' generate_pres_abs_df(oak_hosts, all_oaks)
#' #   species host.status
#' # 2 Quercus_rubra 0
#' # 3 Quercus_robur 1
#' # 4 Quercus_petraea 1
#' # 5 Quercus_alba 0
#' # ...
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

  # Generate host status presence/absence dataframe
  host_status <- data.frame(species = species_list,
                            host.status = 0)
  host_status$host.status[host_status[, 1] %in% observed_hosts] <- 1

  host_status
}


#' Retrieve GBIF taxon keys for a list of species
#'
#' This function is a lightweight wrapper around `rgbif::name_lookup()`. It
#' takes a character vector of species names and queries GBIF to retrieve their
#' GBIF taxon keys, and returns a dataframe indicating the GBIF key for each
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
