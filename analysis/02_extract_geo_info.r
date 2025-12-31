#### SET ENVIRONMENT ####
# Custom functions
source("R/phylo_functions.r")
source("R/geo_functions.r")

# Packages
library(dplyr)

# Seed (not needed for actual analyses, only if exact reproducibility of
# results is desired)
set.seed(24601)


#### EXTRACT GBIF OAK DATA ####
# Read in tree
# See Hipp et al., 2020; https://github.com/andrew-hipp/global-oaks-2019
oak_phylo <- ape::read.tree("data/input/tr.singletons.correlated.1.taxaGrepCrown_accepted_names.tre")

# Clean tree
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

# Get GBIF keys
oak_gbif_keys <- get_gbif_keys(oak_phylo$tip.label, higher_taxon = "220")

# Add keys manually for spp. with missing keys
# Species names: oak_gbif_keys$species[which(is.na(oak_gbif_keys[, 2]))]
# Note that:
# * Q. frainetto: synonym according to GBIF (but accepted by WCVP)
# * Q. faginea, Q. petraea, Q. pyrenaica, Q. imbricaria, Q. rubra, Q. palustris,
# and Q. acuta: function above recognises another spp. first instead
# * Q. ×crenata: name does not match
# * Q. litoralis: Atuna excelsa subsp. excelsa in GBIF
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
  is.na(oak_gbif_keys$key) & oak_gbif_keys$species %in% wcvp_v10_names(manual_oak_keys),
  manual_oak_keys[oak_gbif_keys$species],
  oak_gbif_keys$key
)

# Request GBIF data
# rgbif::occ_download(rgbif::pred("taxonKey", 2877951),
#                     format = "SIMPLE_CSV",
#                     user = "<user>",
#                     email = "<email>",
#                     pwd = "<password>")
# rgbif::occ_download_wait('0361079-210914110416597')
oak_geo_initial <- rgbif::occ_download_get("0361079-210914110416597",
                                           path = "data/tmp") |>
  rgbif::occ_download_import(path = "data/tmp") |>
  dplyr::filter(speciesKey %in% oak_gbif_keys$key & !is.na(speciesKey)) |>
  dplyr::mutate(eventDate = as.character(eventDate))

# Find info for missing species
# Note that:
# * Q. corrugata: in DS as part of some Q. lancifolia entries
# * Q. confertifolia: in DS as part of some Q. crassipes entries
# * Q. sagrana: no occurrences in GBIF
# * Q. calophylla: synonym in GBIF (need to use taxonKey)
# * Q. sartorii:synonym in GBIF (need to use taxonKey)
# * Q. litoralis: Atuna excelsa subsp. excelsa in GBIF
missing_keys <- setdiff(unique(oak_gbif_keys$key),
                        unique(oak_geo_initial$speciesKey))
manual_keys <- c(
  "Quercus corrugata"        = missing_keys[2],
  "Quercus sagrana"          = missing_keys[3],
  "Quercus vacciniifolia"    = missing_keys[4],
  "Quercus calophylla"       = missing_keys[5],
  "Quercus confertifolia"    = missing_keys[6],
  "Quercus sartorii"         = missing_keys[7],
  "Quercus acerifolia"       = missing_keys[8],
  "Quercus palustris"        = missing_keys[9]
)

# Request GBIF data for oak species in phylogeny
# selected_keys <- c(4, 5, 7, 8, 9)
# oak_requests <- lapply(
#   manual_keys[selected_keys],
#   function(tk) {
#     rgbif::occ_download(rgbif::pred("taxonKey", tk),
#                         format = "SIMPLE_CSV",
#                         user = "<user>",
#                         email = "<email>",
#                         pwd = "<password>")
#   }
# )

manual_oak_downloads <- list(
  "Quercus vacciniifolia" = "0062534-250920141307145",
  "Quercus calophylla"    = "0361458-210914110416597",
  "Quercus sartorii"      = "0361555-210914110416597",
  "Quercus acerifolia"    = "0062550-250920141307145",
  "Quercus litoralis"     = "0361560-210914110416597"
)

# Retrieve and combine all manually downloaded species
manual_oak_geo <- lapply(manual_oak_downloads, function(download_key) {
  rgbif::occ_download_get(download_key, path = "data/tmp") |>
    rgbif::occ_download_import(path = "data/tmp")
})

# Combine with existing data
oak_geo <- do.call(rbind, c(list(oak_geo_initial), manual_oak_geo))


#### EXTRACT GBIF NON-OAK HOST DATA ####
# Non-oak host species for Agrilus spp. in DS
plant_hosts <- read.table("data/input/non_quercus_hosts.tsv",
                          header = TRUE, sep = "\t")

plant_gbif_keys <- get_gbif_keys(unique(plant_hosts$plant_sp),
                                 higher_taxon = "220")

# Add keys manually for spp. with missing keys (Betula pendula, function above
# recognises another spp. first instead)
# Species names: plant_gbif_keys$species[which(is.na(plant_gbif_keys[, 2]))]
plant_gbif_keys$key[plant_gbif_keys$species == "Betula pendula"] <- 5331916

# Request GBIF data for non-oak hosts
# plant_requests <- lapply(
#   plant_gbif_keys$key,
#   function(tk) {
#     rgbif::occ_download(rgbif::pred("taxonKey", tk),
#                         format = "SIMPLE_CSV",
#                         user = "<user>",
#                         email = "<email>",
#                         pwd = "<password>")
#   }
# )

manual_plant_downloads <- c(
  "Notholithocarpus densiflorus" = "0361580-210914110416597",
  "Betula pendula"               = "0361582-210914110416597",
  "Castanea sativa"              = "0361584-210914110416597",
  "Corylus avellana"             = "0361600-210914110416597",
  "Fagus sylvatica"              = "0361610-210914110416597",
  "Ostrya carpinifolia"          = "0361624-210914110416597",
  "Castanea dentata"             = "0361625-210914110416597",
  "Fagus grandifolia"            = "0361626-210914110416597",
  "Salix nigra"                  = "0361628-210914110416597",
  "Diospyros virginiana"         = "0361642-210914110416597",
  "Gleditsia triacanthos"        = "0361643-210914110416597",
  "Euonymus europaeus"           = "0361646-210914110416597",
  "Carpinus betulus"             = "0361649-210914110416597",
  "Ficus carica"                 = "0361653-210914110416597",
  "Aesculus glabra"              = "0361656-210914110416597",
  "Aesculus pavia"               = "0361661-210914110416597",
  "Betula lenta"                 = "0361664-210914110416597",
  "Carpinus caroliniana"         = "0361668-210914110416597",
  "Ostrya virginiana"            = "0361677-210914110416597"
)

# Retrieve and combine all manually downloaded species
manual_plant_geo <- lapply(manual_plant_downloads, function(download_key) {
  rgbif::occ_download_get(download_key, path = "data/tmp") |>
    rgbif::occ_download_import(path = "data/tmp")
})

# Combine with existing data
plant_geo <- do.call(rbind, c(list(oak_geo), manual_plant_geo))

# Write table
# Note, file not included in data/ due to size (4 GB)
# write.table(plant_geo, "data/results/gbif_geo_plants.tsv",
#             quote = FALSE,
#             row.names = FALSE,
#             col.names = TRUE,
#             sep = "\t")


#### CLEAN GBIF GEO DATASET ####
# See https://cran.r-project.org/web/packages/CoordinateCleaner/vignettes/Cleaning_GBIF_data_with_CoordinateCleaner.html#:~:text=We%20might%20also%20want%20to,it%20with%20the%20following%20code

# Note, file not included in data/ due to size (4 GB)
# plant_geo <- read.table("data/results/gbif_geo_plants.tsv",
#                         header = TRUE,
#                         sep = "\t",
#                         quote = "",
#                         comment.char = "")

plant_geo_clean <- as.data.frame(plant_geo)

# Check missing species
plant_spp <- c(oak_gbif_keys$species, plant_gbif_keys$species)
plant_spp[!(plant_spp %in% unique(plant_geo_clean$species))]

# Re-name species that are reported as scientificName under a different species
# * Q. sagrana has no entries in GBIF
# * Q. corrugata is reported as some Q. lancifolia entries (scientificName)
# * Q. confertifolia is reported as some Q. crassipes entries (scientificName)
plant_geo_clean$species[plant_geo_clean$species == "Quercus lancifolia"
                        & plant_geo_clean$scientificName ==
                          "Quercus corrugata Hook."] <- "Quercus corrugata"
plant_geo_clean$species[plant_geo_clean$species == "Quercus crassipes"
                        & plant_geo_clean$scientificName ==
                          "Quercus confertifolia Bonpl."] <- "Quercus confertifolia"

# Rename species reported as a different species
# * Q. calophylla is reported as Q. candicans
# * Q. sartorii is reported as Q. xalapensis
# * Q. litoralis (Atuna excelsa) is reported as Atuna excelsa
replacements <- c("Quercus candicans" = "Quercus calophylla",
                  "Quercus xalapensis" = "Quercus sartorii",
                  "Atuna excelsa" = "Quercus litoralis (Atuna excelsa)",
                  "Quercus crenata" = "Quercus ×crenata")
for (original_name in wcvp_v10_names(replacements)) {
  plant_geo_clean$species <- gsub(original_name,
                                  replacements[original_name],
                                  plant_geo_clean$species)
}

# Drop for which the "scientificName" is not a synonym of the accepted species
# according to WCVP. This is a conservative yet imperfect approach, but errors
# are unlikely to have a significant impact
gbif_invalid <- read.table("data/input/gbif_synonyms_not_accepted_wcvp.tsv",
                           header = TRUE, sep = "\t")
gbif_invalid$scientificName <- gsub(" [(A-Z].*", "", gbif_invalid$scientificName)

plant_geo_clean$verbatimClean <- clean_taxon_name(plant_geo$verbatimScientificName)

plant_geo_clean <- plant_geo_clean |>
  anti_join(gbif_invalid,
            by = c("species" = "species",
                   "verbatimClean" = "scientificName"))

# Missing spp: Q. sagrana (no records)
plant_spp[!(plant_spp %in% unique(plant_geo_clean$species))]

# Filter data based on coords, counts, year, presence, and common issues
plant_geo_clean <- filter_gbif_data(plant_geo_clean,
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
                                               "FOSSIL_SPECIMEN"))

# Missing spp: Q. sagrana (no records) and Q. yiwuensis (no coordinate info)
plant_spp[!(plant_spp %in% unique(plant_geo_clean$species))]

# Prepare data for coordinate cleaner by converting country code to ISO3c
plant_geo_clean$countryCode <- countrycode::countrycode(plant_geo_clean$countryCode,
                                                        origin =  "iso2c",
                                                        destination = "iso3c")

# Use CoordinateCleaner to flag problematic entries. Issues being flagged:
# * CAPITALS: tests a radius around adm-0 capitals
# * CENTROIDS: tests a radius around country & province centroids
# * COUNTRIES: tests if coords are from the country in the country column;
# SWITCHED OFF as it flags up many records close to borders
# * DUPLICATES: tests for duplicate records (identical coordinates)
# * EQUAL: tests for equal absolute longitude and latitude.
# * GBIF: tests a one-degree radius around the GBIF headquarters
# * INSTITUTIONS: tests a radius around known biodiversity institutions
# * OUTLIERS: tests each species for outlier records; outliers_size = min.
#  number of records in a dataset to run this test
# * SEAS: tests if coordinates fall into the ocean
# * ZEROS: tests for plain zeros, equal lat and lon, and a radius around 0/0
rownames(plant_geo_clean) <- seq_len(nrow(plant_geo_clean))

# Apply cleaning per species (to avoid segmentation fault)
flag_species <- split(plant_geo_clean, plant_geo_clean$species)

flags_list <- lapply(flag_species, function(df) {
  CoordinateCleaner::clean_coordinates(x = df,
                                       lon = "decimalLongitude",
                                       lat = "decimalLatitude",
                                       species = "species",
                                       capitals_rad = 1000,
                                       centroids_rad = 1000,
                                       outliers_size = 50,
                                       inst_rad = 1000,
                                       zeros_rad = 1,
                                       tests = c("capitals",
                                                 "centroids",
                                                 "duplicates",
                                                 "equal",
                                                 "outliers",
                                                 "gbif",
                                                 "institutions",
                                                 "zeros",
                                                 "seas"))
})

# Combine results back into a single data frame
flags <- do.call("rbind", flags_list)

# Remove any flagged entries
# ~ half of data flagged (most are duplicates): summary(flags)
plant_geo_clean <- flags[flags$.summary, ]

# Missing spp: Q. sagrana (no records) and Q. yiwuensis (no coordinate info)
plant_spp[!(plant_spp %in% unique(plant_geo_clean$species))]

# Remove non-relevant fields, only retaining species, latitude, and longitude colums
plant_geo_clean <- plant_geo_clean[, c("species",
                                       "decimalLongitude",
                                       "decimalLatitude")]
colnames(plant_geo_clean) <- c("plant_sp", "lon", "lat")

# Subsample all spp. to =< 100k entries (so that downstream analyses can run)
plant_geo_clean <- plant_geo_clean |>
  group_by(plant_sp) |>
  slice_sample(n = 100000) |>
  ungroup()

# Add a manual entry (country centroid) for missing species
# * Q. sagrana (Cuba: powo.science.kew.org/taxon/urn:lsid:ipni.org:names:216368-2)
plant_geo_clean <- add_species_centroid(df = plant_geo_clean,
                                        species_name = "Quercus sagrana",
                                        iso3 = "CUB",
                                        col_species = "plant_sp")

# * Q. yiwuensis (SC China: powo.science.kew.org/taxon/urn:lsid:ipni.org:names:360253-1)
plant_geo_clean <- add_species_centroid(df = plant_geo_clean,
                                        species_name = "Quercus yiwuensis",
                                        iso3 = "CHN",
                                        region = "Hunan",
                                        col_species = "plant_sp")

# Check missing species
plant_spp[!(plant_spp %in% unique(plant_geo_clean$plant_sp))]

# Write table
# write.table(plant_geo_clean,
#             "data/results/gbif_plants_clean.tsv",
#             quote = FALSE,
#             sep = "\t",
#             col.names = TRUE,
#             row.names = FALSE)
