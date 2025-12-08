##########################################################################################
#####                           AGRILUS OAK HOSTS ANALYSES                          ######
#####                           ELVIRA HERNANDEZ GUTIERREZ                          ######
##########################################################################################

#####       SET ENVIRONMENT                                                          #####
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
#####       RETRIEVE COORDINATE INFO FROM GBIF                                       #####
# A) DOWNLOAD OAK SPP. INFO FOR ALL OAK SPP. PRESENT IN THE PHYLOGENY
# Download information for all oak spp. in GBIF (all downloads take < 1 min.):
# NB, TaxonKey for Quercus in GBIF == 2877951
# occ_download(pred("taxonKey", 2877951), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361079-210914110416597')

geo_info <- occ_download_get('0361079-210914110416597',
                             path = "tmp/gbif") %>%
            occ_download_import(quote = "", path = "tmp/")

# Read the phylogenetic tree of oaks and extract spp. names
# This tree also comes from from Hipp et al. (2019)
oak_tree <- read.tree("analysis/input/quercus_hipp19_singleton_crown_sp_level.tre")
plants <- gsub(oak_tree$tip.label, pattern = "\\|.*", replacement = "")
plants <- gsub(plants, pattern = "_", replacement = " ", fixed = TRUE)

# Check no. of oak spp. in the phylogeny (238) and see if they are all in the GBIF DF:
length(plants)

# Extract their GBIF info (takes ca. 1 min to run)
plant_keys_info <- lapply(plants, name_lookup, rank = "SPECIES",
                          higherTaxonKey = "220", status = "Accepted", limit = 1)
# plant_keys_info[[1]]

# Create a dataframe to store plant species name and their GBIF key. If there is no key
# info for the current sp., save its index in manual_keys. Takes < 1 min to run
plant_keys <- as.data.frame(matrix(ncol = 2, nrow = length(plant_keys_info)))
colnames(plant_keys) <- c("plant", "key")
manual_keys <- c()
for (i in 1:length(plant_keys_info)) {
  if (!is.null(plant_keys_info[[i]][4]$hierarchies)) {
    sp <- plant_keys_info[[i]][2]$data$species
    # So that it doesn't include ""Quercus faginea-mirbeckii" (sp. in "plants" is
    # "Quercus faginea")
    if (sp %in% plants) {
      plant_keys[i,1] <- sp
      plant_keys[i,2] <- plant_keys_info[[i]][2]$data$key
    } else {
      manual_keys <- cbind(manual_keys, i)
    }
  } else {
    manual_keys <- cbind(manual_keys, i)
  }
}

# For all spp. without a key, look them up on GBIF and add them manually (7 oak spp.).
# * = name is a syn. according to GBIF but accepted by WCVP
# ** == Atuna excelsa subsp. excelsa
# *** = loop above recognises another spp. first instead
plants[manual_keys]
no <- grep(plants[manual_keys][1], x = plants)  # "Quercus new species"
plant_keys[no,1] <- plants[manual_keys][1]      # "Quercus new species"
plant_keys[no,2] <- NA                          # "Quercus new species"
no <- grep(plants[manual_keys][2], x = plants)  # "Quercus corrugata"
plant_keys[no,1] <- plants[manual_keys][2]      # "Quercus corrugata"
plant_keys[no,2] <- 2880406                     # "Quercus corrugata"*
no <- grep(plants[manual_keys][3], x = plants)  # "Quercus faginea"
plant_keys[no,1] <- plants[manual_keys][3]      # "Quercus faginea"
plant_keys[no,2] <- 2881480                     # "Quercus faginea"***
no <- grep(plants[manual_keys][4], x = plants)  # "Quercus calophylla"
plant_keys[no,1] <- plants[manual_keys][4]      # "Quercus calophylla"
plant_keys[no,2] <- 2880691                     # "Quercus calophylla"*
no <- grep(plants[manual_keys][5], x = plants)  # "Quercus confertifolia"
plant_keys[no,1] <- plants[manual_keys][5]      # "Quercus confertifolia"
plant_keys[no,2] <- 2878551                     # "Quercus confertifolia"*
no <- grep(plants[manual_keys][6], x = plants)  # "Quercus sartorii"
plant_keys[no,1] <- plants[manual_keys][6]      # "Quercus sartorii"
plant_keys[no,2] <- 2880720                     # "Quercus sartorii"*
no <- grep(plants[manual_keys][7], x = plants)  # "Quercus palustris"
plant_keys[no,1] <- plants[manual_keys][7]      # "Quercus palustris"
plant_keys[no,2] <- 8313153                     # "Quercus palustris"
no <- grep(plants[manual_keys][8], x = plants)  # "Quercus sp. nov."
plant_keys[no,1] <- plants[manual_keys][8]      # "Quercus sp. nov."
plant_keys[no,2] <- NA                          # "Quercus sp. nov."
no <- grep(plants[manual_keys][9], x = plants)  # "Quercus litoralis"
plant_keys[no,1] <- plants[manual_keys][9]      # "Quercus litoralis"
plant_keys[no,2] <- 10839442                    # "Quercus litoralis"**
no <- grep(plants[manual_keys][10], x = plants) # "Quercus acuta"
plant_keys[no,1] <- plants[manual_keys][10]     # "Quercus acuta"
plant_keys[no,2] <- 2876863                     # "Quercus acuta"***

nrow(plant_keys) == length(plants)

# Subset geo_info so as to only keep those in the phylogeny
geo_info <- subset(geo_info, speciesKey %in% plant_keys$key & !is.na(speciesKey))
# unique(geo_info[is.na(geo_info$speciesKey),]$species)
geo_info <- subset(geo_info, speciesKey %in% plant_keys$key)
setdiff(geo_info$speciesKey, plant_keys$key)
setdiff(geo_info$species, plant_keys$plant)
missing_keys <- setdiff(unique(plant_keys$key), unique(geo_info$speciesKey))
missing_sp <- setdiff(unique(plant_keys$plant), unique(geo_info$species))

# For those species not in geo_info, download their info manually:
# Quercus new sp.
missing_sp[1]
missing_keys[1]

# Quercus corrugata: - Don't run! It's already in DS as part of some Q. lancifolia entries
# missing_sp[2]
# # "taxonKey" (instead of "speciesKey") because it is a syn. according to GBIF
# # occ_download(pred("taxonKey", missing_keys[2]), format = "SIMPLE_CSV",
# #              user = "", email = "", pwd = "")
# # occ_download_wait('0361453-210914110416597')
# add_plant <- occ_download_get('0361453-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# geo_info <- rbind(geo_info, add_plant)

# Quercus sagrana: no occurrences in GBIF
missing_sp[3]
missing_keys[3]

# Quercus calophylla
missing_sp[4]
# "taxonKey" because it is a syn. according to GBIF
# occ_download(pred("taxonKey", missing_keys[4]), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361458-210914110416597')
add_plant <- occ_download_get('0361458-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Quercus confertifolia: Don't run! It's already in DS as part of some Q. crassipes entries
# missing_sp[5]
# # "taxonKey" because it is a syn. according to GBIF
# # occ_download(pred("taxonKey", missing_keys[5]), format = "SIMPLE_CSV",
# #              user = "", email = "", pwd = "")
# # occ_download_wait('0361447-210914110416597')
# add_plant <- occ_download_get('0361547-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# geo_info <- rbind(geo_info, add_plant)

# Quercus sartorii
# "taxonKey" because it is a syn. according to GBIF
missing_sp[6]
# occ_download(pred("taxonKey", missing_keys[6]), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361555-210914110416597')
add_plant <- occ_download_get('0361555-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Quercus sp.nov.
missing_sp[7]
missing_keys[1]

# Quercus litoralis (NB, key is actually == 10839442, as this is the key for
# Atuna excelsa subsp. excelsa, which is the accepted name in GBIF)
# "taxonKey" because it is a subsp. according to GBIF
missing_sp[8]
# occ_download(pred("taxonKey", 10839442), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361560-210914110416597')
add_plant <- occ_download_get('0361560-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)



# B) DOWNLOAD PLANT SPP. INFO FOR ALL NON-OAK SPP. THAT HOST AGRILUS THAT USE OAKS
# Read in file with info on other plants that plant Agrilus spp. used by oaks
# NB, I have only retained info for those plants for which there is info to sp. level
plants <- read.table("analysis/input/non_oak_hosts.tsv", header = T, sep = "\t")

# Check no. of non-oak plant spp. used by Agrilus spp. that exploit oaks (19).
length(unique(plants$hosts))

# Extract GBIF info (takes ca. 1 min to run).
plant_keys_info <- lapply(unique(plants$hosts), name_lookup, rank = "SPECIES",
                          higherTaxonKey = "220", status = "Accepted", limit = 1)
# plant_keys_info[[1]]

# Create a dataframe to store plant species name and their GBIF key. If there is no key
# info for the current sp., save its index in manual_keys. Takes < 1 min to run
plant_keys <- as.data.frame(matrix(ncol = 2, nrow = length(plant_keys_info)))
colnames(plant_keys) <- c("plant", "key")
manual_keys <- c()
for (i in 1:length(plant_keys_info)) {
  if (!is.null(plant_keys_info[[i]][4]$hierarchies)) {
    plant_keys[i,1] <- plant_keys_info[[i]][2]$data$species
    plant_keys[i,2] <- plant_keys_info[[i]][2]$data$key
  } else {
    manual_keys <- cbind(manual_keys, i)
  }
}
manual_keys

# Download the info for these extra plant spp. and add them to geo_info
# Notholithocarpus densiflorus
plant_keys[1,]
# occ_download(pred("taxonKey", plant_keys[1,2]), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361580-210914110416597')
add_plant <- occ_download_get('0361580-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Betula pendula
plant_keys[2,]
# occ_download(pred("taxonKey", plant_keys[2,2]), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361582-210914110416597')
add_plant <- occ_download_get('0361582-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Castanea sativa
plant_keys[3,]
# occ_download(pred("taxonKey", plant_keys[3,2]), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361584-210914110416597')
add_plant <- occ_download_get('0361584-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Corylus avellana
plant_keys[4,]
# occ_download(pred("taxonKey", c(plant_keys[4,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361600-210914110416597')
add_plant <- occ_download_get('0361600-210914110416597', path = "\tmp") %>%
             occ_download_import(quote = "", path = "\tmp")
geo_info <- rbind(geo_info, add_plant)

# Fagus sylvatica
plant_keys[5,]
# occ_download(pred("taxonKey", c(plant_keys[5,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361610-210914110416597')
add_plant <- occ_download_get('0361610-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Ostrya carpinifolia
plant_keys[6,]
# occ_download(pred("taxonKey", c(plant_keys[6,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361624-210914110416597')
add_plant <- occ_download_get('0361624-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Castanea dentata
plant_keys[7,]
# occ_download(pred("taxonKey", c(plant_keys[7,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361625-210914110416597')
add_plant <- occ_download_get('0361625-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Fagus grandifolia
plant_keys[8,]
# occ_download(pred("taxonKey", c(plant_keys[8,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361626-210914110416597')
add_plant <- occ_download_get('0361626-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Salix nigra
plant_keys[9,]
# occ_download(pred("taxonKey", c(plant_keys[9,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361628-210914110416597')
add_plant <- occ_download_get('0361628-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Diospyros virginiana
plant_keys[10,]
# occ_download(pred("taxonKey", c(plant_keys[10,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361642-210914110416597')
add_plant <- occ_download_get('0361642-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Gleditsia triacanthos
plant_keys[11,]
# occ_download(pred("taxonKey", c(plant_keys[11,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361643-210914110416597')
add_plant <- occ_download_get('0361643-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Euonymus europaeus
plant_keys[12,]
# occ_download(pred("taxonKey", c(plant_keys[12,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361646-210914110416597')
add_plant <- occ_download_get('0361646-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Carpinus betulus
plant_keys[13,]
# occ_download(pred("taxonKey", c(plant_keys[13,2])), format = "SIMPLE_CSV",
#              user = "", email = "",
#              pwd = "")
# occ_download_wait('0361649-210914110416597')
add_plant <- occ_download_get('0361649-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Ficus carica
plant_keys[14,]
# occ_download(pred("taxonKey", c(plant_keys[14,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361653-210914110416597')
add_plant <- occ_download_get('0361653-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)


# Aesculus glabra
plant_keys[15,]
# occ_download(pred("taxonKey", c(plant_keys[15,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361656-210914110416597')
add_plant <- occ_download_get('0361656-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Aesculus pavia
plant_keys[16,]
# occ_download(pred("taxonKey", c(plant_keys[16,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361661-210914110416597')
add_plant <- occ_download_get('0361661-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Betula lenta
plant_keys[17,]
# occ_download(pred("taxonKey", c(plant_keys[17,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361664-210914110416597')
add_plant <- occ_download_get('0361664-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Carpinus caroliniana
plant_keys[18,]
# occ_download(pred("taxonKey", c(plant_keys[18,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361668-210914110416597')
add_plant <- occ_download_get('0361668-210914110416597', path = "tmp/gbif") %>%
             occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Ostrya virginiana
plant_keys[19,]
# occ_download(pred("taxonKey", c(plant_keys[19,2])), format = "SIMPLE_CSV",
#              user = "", email = "", pwd = "")
# occ_download_wait('0361677-210914110416597')
add_plant <- occ_download_get('0361677-210914110416597', path = "tmp/gbif") %>%
  occ_download_import(quote = "", path = "tmp/")
geo_info <- rbind(geo_info, add_plant)

# Write table
# write.table(geo_info, "analysis/results/gbif_plants.csv", quote = F, row.names = F,
#             col.names = T, sep = "\t")

# NB, check for absolute duplicates on commandline with:
# sort analysis/results/gbif_plants.csv | uniq -c | sort -nr | grep '^\t* *2' | cut -f2



#
#####       CLEAN PLANT GEO DATA FROM GBIF                                           #####
# See https://cran.r-project.org/web/packages/CoordinateCleaner/vignettes/
# Cleaning_GBIF_data_with_CoordinateCleaner.html

# A) EXTRACT PLANT GEO DATA AND REMOVE ENTRIES IN DISCORDANCE WITH WCVP
# Extract data (created above)
geo_info <- read.table("analysis/results/gbif_plants.csv", header = T, sep = "\t",
                       quote = "", comment.char = '')

# Remove any entries that are in gbif_notwcvp_entries.txt, i.e., those for which the
# name under which the sp. was originally reported ("scientificName") is regarded as a
# syn. for the current sp. ("species") by GBIF but not WCVP. Takes ca. 30 min to run.

# NB, conservative approach, some of these species might actually be synonyms... (but
# then they all have very few occurrence records so not likely to make much of an impact)
gbif_notwcvp <- read.table("analysis/input/gbif_notwcvp_auth.txt", header = F, sep = "\t")
head(gbif_notwcvp)
for (i in 1:nrow(geo_info)) {
  if (i%%250000 == 0) {print(i)}
  sp <- geo_info$verbatimScientificName[i]
  sp <- gsub(pattern = "??|-|subsp\\. | var\\.|f\\.| de | ex ", "",  sp)
  # sp <- gsub(pattern = "([A-Z][a-z]* [a-z][a-z]* *[a-z]*[a-z]) *.*", "\\1",  sp)
  sp <- gsub(pattern = "  *", " ",  sp)
  if (sp %in% gbif_notwcvp) {
    # print(sp)
    geo_info$scientificName[i] <- NA
  }
}
tail(geo_info[is.na(geo_info$scientificName), c(10,13,14)], 2)
geo_info <- geo_info[!is.na(geo_info$scientificName),]


# B) CLEAN GEO_INFO
# i) Rename spp. reported under a synonym (i.e., the sp. is a syn. of sth
# else acc. to GBIF, but it's actually an accepted name acc. to WCVP)*
# * I removed these spp. from analysis/input/gbif_notwcvp_auth.txt
# intentionally so that they wouldn't get removed in the previous step

# Get name of all plants_all that should be in my dataset
plants_all <- read.tree("analysis/input/quercus_hipp19_singleton_crown_sp_level.tre")
plants_all <- gsub(plants_all$tip.label, pattern = "Q.", replacement = "Quercus",
                   fixed = T)
plants_all <- gsub(plants_all, pattern = "_", replacement = " ",
                   fixed = T)
plants_all <- gsub(plants_all, pattern = " -Atuna excelsa", replacement = "",
                   fixed = T)

plants_all <- gsub(plants_all, pattern = "?? ", replacement = "", fixed = T)
plants_all <- c(plants_all, unique(read.table("analysis/input/non_oak_hosts.tsv",
                                              header = T, sep = "\t")$hosts))
length(plants_all) == (238 + 19)  # 238 Quercus spp. + 19 extra hosts

# Are they all in geo_info?
plants_all[!(plants_all %in% unique(geo_info$species))]

# Re-name species that are reported as something else
# NB, Q. new species, Q. sagrana and Q. sp.nov. have no entries in GBIF.
# Q. corrugata is reported as Q. lancifolia (which is also in my dataset):
# nrow(subset(geo_info, species == "Quercus lancifolia" &
#               scientificName != "Quercus corrugata Hook."))
geo_info$species[geo_info$species == "Quercus lancifolia" & geo_info$scientificName ==
                 "Quercus corrugata Hook."] <- "Quercus corrugata"

# Q. calophylla is reported as Q. candicans:
# unique(subset(geo_info, taxonKey == 2880691)$species)
geo_info$species <- gsub(pattern = "Quercus candicans",
                         replacement = "Quercus calophylla", x = geo_info$species)

# Q. confertifolia is reported as Q. crassipes (which is also in my dataset):
# unique(subset(geo_info, taxonKey == 2878551)$species)
# nrow(subset(geo_info, species == "Quercus crassipes" &
#             scientificName != "Quercus confertifolia Bonpl."))
geo_info$species[geo_info$species == "Quercus crassipes" & geo_info$scientificName ==
                 "Quercus confertifolia Bonpl."] <- "Quercus confertifolia"

# Q. sartorii is reported as Q. xalapensis:
# unique(subset(geo_info, taxonKey == 2880720)$species)
geo_info$species <- gsub(pattern = "Quercus xalapensis",
                         replacement = "Quercus sartorii", x = geo_info$species)

# Q. litoralis is reported as Atuna excelsa:
geo_info$species <- gsub(pattern = "Atuna excelsa",
                         replacement = "Quercus litoralis", x = geo_info$species)

# Are they all now in geo_info?
# Yes, except for Q. sagrana (no records)
plants_all[!(plants_all %in% unique(geo_info$species))]


# ii) Only keep records of interest
# Only keep entries with coordinate info AND that present a coordinate uncertainty <= 5 km
# length(geo_info$coordinateUncertaintyInMeters[
#        is.na(geo_info$coordinateUncertaintyInMeters)])
# hist(geo_info$coordinateUncertaintyInMeters/1000, breaks = 1000,
#      xlim = c(0, 200), ylim = c(0,10000))
geo_info <- subset(geo_info, !is.na(decimalLatitude) & !is.na(decimalLongitude))
geo_info <- subset(geo_info, is.na(coordinateUncertaintyInMeters) |
                     coordinateUncertaintyInMeters <= 5000)
length(unique(geo_info$species))
nrow(geo_info)

# Remove entries with high individualCounts
# See https://www.gbif.org/data-quality-requirements-occurrences#dcCount
# hist(geo_info$individualCount[!is.na(geo_info$individualCount)])
# length(geo_info$individualCount[geo_info$individualCount >= 100 &
#        !is.na(geo_info$individualCount)])
geo_info <- subset(geo_info, is.na(individualCount) | individualCount <= 10)
geo_info <- subset(geo_info, is.na(individualCount) | individualCount > 0)
length(unique(geo_info$species))

# Remove entries with occurrenceStatus == "ABSENT"
geo_info <- subset(geo_info, occurrenceStatus == "PRESENT")
nrow(geo_info)

# Remove any entries that come from fossil records
geo_info <- subset(geo_info, basisOfRecord != "FOSSIL_SPECIMEN")
geo_info <- subset(geo_info, issue != "FOSSIL_SPECIMEN")
nrow(geo_info)

# Remove any entries where the coordinate info is unreliable
# To see all issues in data, run (bash):
# cut -f18 analysis/results/gbif_plants.csv | sort | uniq | sed -'s/;/\n/g' | sort | uniq
# NB, for issues see https://data-blog.gbif.org/post/issues-and-flags/
issues <- c("COORDINATE_INVALID", "COORDINATE_OUT_OF_RANGE",
            "COORDINATE_REPROJECTION_SUSPICIOUS",
            "COORDINATE_UNCERTAINTY_METERS_INVALID",
            "COUNTRY_COORDINATE_MISMATCH", "PRESUMED_SWAPPED_COORDINATE",
            "PRESUMED_NEGATED_LONGITUDE", "PRESUMED_NEGATED_LATITUDE",
            "ZERO_COORDINATE")
geo_info <- subset(geo_info, !grepl(pattern = paste(issues, collapse =  "|"), issue))
nrow(geo_info)

# Are they all still in geo_info?
# Yes, except for Q. sagrana (no records) and Q. yiwuensis (no coordinate info)
plants_all[!(plants_all %in% unique(geo_info$species))]
length(unique(geo_info$species))

# Remove records pre 1950
geo_info <- subset(geo_info, !is.na(year) & year >= 1950)
nrow(geo_info)

# Are they all still in geo_info?
# Yes, except for Q. sagrana (no records) and Q. yiwuensis (no coordinate info)
plants_all[!(plants_all %in% unique(geo_info$species))]
length(unique(geo_info$species))


# iii) Use CoordinateCleaner
# Check https://mran.microsoft.com/snapshot/2019-02-10/web/packages/CoordinateCleaner/
# vignettes/Tutorial_Cleaning_GBIF_data_with_CoordinateCleaner.html
# First, convert country code from ISO2c to ISO3c
geo_info$countryCode <-  countrycode::countrycode(geo_info$countryCode,
                                                  origin =  'iso2c',
                                                  destination = 'iso3c')

# And plot data to get an overview
wm <- borders("world", colour="gray50", fill="gray50")
ggplot() + coord_fixed()+ wm +
           geom_point(data = geo_info, aes(x = decimalLongitude, y = decimalLatitude),
                      colour = "darkred", size = 0.5)+
           theme_bw()

# Then, flag any problematic entries.
# Issues I am flagging:
# *CAPITALS: tests a radius around adm-0 capitals (radius = capitals_rad, default =
# 10,000 m, here = 1,000 m).
# *CENTROIDS: tests a radius around country centroids (both country and provinces, check
# centroids_detail) (radius = centroids_rad, default = 1,000 m, here = 1,000 m).
# *COUNTRIES: tests if coordinates are from the country in the country column. SWITCHED
# OFF as it flags up many records that are very close to country borders.
# *DUPLICATES: tests for duplicate records (identical coordinates).
# *EQUAL: tests for equal absolute longitude and latitude.
# *GBIF: tests a one-degree radius around the GBIF headquarters in Copenhagen, Denmark.
# *INSTITUTIONS: tests a radius around known biodiversity institutions from institutions
# (radius = inst_rad, default = 100 m, here = 1,000 m).
# *OUTLIERS: tests each species for outlier records. Depending on outliers_method, it
# either flags records that are a min. dist away from all other records (by default,
# outliers_td = 1000 km), or records outside a multiple of the IQR of the min. distances
# to the next record (default; by default, outliers_mtp = 5). The minimum number of
# records in a dataset to run this test is given by outliers_size (default = 7, here
# = 50).
# *SEAS: tests if coordinates fall into the ocean.
# *ZEROS: tests for plain zeros, equal latitude and longitude and a radius in degrees
# around the point 0/0. The radius is zeros_rad (default = 0.5, here = 1).

# Without rownames(geo_info) <- 1:nrow(geo_info) it throws and error message
rownames(geo_info) <- 1:nrow(geo_info)
flags <- clean_coordinates(x = geo_info, lon = "decimalLongitude",
                          lat = "decimalLatitude",
                          # countries = "countryCode",
                          species = "species",
                          capitals_rad = 1000,
                          centroids_rad = 1000,
                          outliers_size = 50,
                          inst_rad = 1000,
                          zeros_rad = 1,
                          tests = c("capitals", "centroids", # "countries",
                                   "duplicates", "equal","outliers", "gbif",
                                    "institutions", "zeros", "seas"))

# Check summary: 52% of data is flagged (48% of all entries are duplicates!!)
summary(flags)
# head(flags[!flags$.con,], 2)
# subset(geo_info, decimalLatitude == 34.73667)[, c("species", "decimalLongitude",
#                                                   "decimalLatitude", "year",
#                                                   "coordinatePrecision", "recordedBy")]

plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

# Remove any flagged entries
geo_info_clean <- geo_info[flags$.summary,]

# # Are they all still in geo_info?
# # Yes, except for Q. sagrana (no records)
# plants_all[!(plants_all %in% unique(geo_info_clean$species))]


# iv) Only keep relevant fields
geo_info_clean <- geo_info_clean[, c("species", "decimalLatitude", "decimalLongitude")]
colnames(geo_info_clean) <- c("plant.sp", "lat", "lon")

length(unique(geo_info_clean$plant.sp))
nrow(geo_info_clean)


# v) Subsample so that no sp. has > 100k entries (or analyses won't run)
sort(table(geo_info_clean$plant.sp)[table(geo_info_clean$plant.sp) >= 100000])
spp <- geo_info_clean$plant.sp
spp <- dimnames(sort(table(spp)[table(spp) >= 100000]))[[1]]

add_all <- subset(geo_info_clean, plant.sp == spp[1])
add <- add_all[sample(nrow(add_all), 100000), ]
add_all <- subset(geo_info_clean, plant.sp == spp[2])
add <- rbind(add, add_all[sample(nrow(add_all), 100000), ])
add_all <- subset(geo_info_clean, plant.sp == spp[3])
add <- rbind(add, add_all[sample(nrow(add_all), 100000), ])
add_all <- subset(geo_info_clean, plant.sp == spp[4])
add <- rbind(add, add_all[sample(nrow(add_all), 100000), ])
add_all <- subset(geo_info_clean, plant.sp == spp[5])
add <- rbind(add, add_all[sample(nrow(add_all), 100000), ])
add_all <- subset(geo_info_clean, plant.sp == spp[6])
add <- rbind(add, add_all[sample(nrow(add_all), 100000), ])

geo_info_clean <- subset(geo_info_clean, !(plant.sp %in% spp))
geo_info_clean <- rbind(geo_info_clean, add)

length(unique(geo_info_clean$plant.sp))
nrow(geo_info_clean)

# vi) Add a "fake" manual entry for the spp. missing (country
# centroid according to CoordinateCleaner)
# Q. sagrana
# Present in Cuba: powo.science.kew.org/taxon/urn:lsid:ipni.org:names:216368-2
add <- data.frame(cbind("Quercus sagrana",
                  subset(countryref, iso3 == "CUB" &
                         type == "country")[1, c("centroid.lat", "centroid.lon")]))
colnames(add) <- colnames(geo_info_clean)
geo_info_clean <- rbind(geo_info_clean, add)

# Q. yiwuensis
# Present in in SC China: powo.science.kew.org/taxon/urn:lsid:ipni.org:names:360253-1
add <- data.frame(cbind("Quercus yiwuensis",
                  subset(countryref, iso3 == "CHN" &
                         name == "Hunan")[1, c("centroid.lat", "centroid.lon")]))
colnames(add) <- colnames(geo_info_clean)
geo_info_clean <- rbind(geo_info_clean, add)
plants_all[!(plants_all %in% unique(geo_info_clean$plant.sp))]
nrow(geo_info_clean)

# Write table
# write.table(geo_info_clean, "analysis/results/gbif_plants_clean.tsv", quote = F, sep = "\t",
#             col.names = T, row.names = F)



#