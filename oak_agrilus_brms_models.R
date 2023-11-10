##########################################################################################
#####                           AGRILUS OAK HOSTS ANALYSES                          ######
#####                           ELVIRA HERNANDEZ GUTIERREZ                          ######
##########################################################################################

#####       SET ENVIRONMENT                                                          #####
## Set directory and clean environment
setwd("~/oak_analyses/01_original_models")   ## Whatever your working directory is
pw <- "/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/"
# rm(list = ls())
# rm(list=setdiff(ls(), c("geo_info", "interaction_data")))

# invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach,
#                  character.only=TRUE, unload=TRUE))

## Load packages
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
##### I)    COMPUTE FRITZ AND PURVIS' D METRIC (PHYLOGENETIC SIGNAL)                 #####
## A) READ IN AND PREPARE PHYLOGENETIC TREE OF OAKS
## Source: Hipp et al. (2019); doi = 10.1111/nph.16162
## See https://github.com/andrew-hipp/global-oaks-2019
## "For dating, samples were pruned to one sample per named species, favoring samples
## with the most loci, except for species in which variable position of samples from
## different populations was deemed to represent cryptic diversity, in which case more
## than one exemplar was retained. The resulting ???singletons tree??? was estimated in
## RAXML using a phylogenetic constraint."
## Read tree (NB, with tmp/tr.*Stem_accepted*.tre results are very similar)
oak_tree <- read.tree(paste0(pw,
                             "tmp/tr.singletons.correlated.1.",
                             "taxaGrepCrown_accepted_names.tre"))

## Check it's rooted (yes)
is.rooted(oak_tree)

## Plot with support values
# plot(oak_tree)
# drawSupportOnEdges(oak_tree$node.label)

## Check no. of tips
# length(oak_tree$tip.label)

## Remove duplicated tips that form a sister clade
## Dulplicated tips (print tree to see, also included subsp. & var.) are:
## "Q. magnoliifolia"   - merge             "Q. crassifolia"     - merge
## "Q. agrifolia"       - merge             "Q. coccifera"       - merge
## "Q. arizonica"       - keep              "Q. laeta"           - keep
## "Q. conzattii"       - keep              "Q. benthamii"       - keep
## "Q. rehderiana"      - keep              "Q. ithanurensis"    - keep
## "Q. parvula"         - keep              "Q. robur"           - keep
## "Q. infectoria"      - keep              "Q. magnoliifolia"   - keep
## "Q. petraea"         - keep
# oak_tree$tip.label[duplicated(oak_tree$tip.label)]
oak_tree <- drop.tip(oak_tree,
                     c("Quercus_magnoliifolia|Mexico|Durango|2017124|QUE002013",
                       "Quercus_crassifolia|Mexico|Durango|2017108|QUE001997",
                       "Quercus_agrifolia|USA|CA|CA-DAV-MH55|QUE000328",
                       "Quercus_coccifera|Italy|Apulia|TUS13-013|QUE001988"))

## Rename tips (only keep Q. spp info)
for (i in 1:length(oak_tree$tip.label)) {
  ## Extract tip label
  tip <- oak_tree$tip.label[i]
  ## Only keep species name info, and change "Quercus" for "Q."
  oak_spp <- gsub(oak_tree$tip.label[i], pattern = "Quercus_([^|]*).*",
                  replacement = "Q._\\1")
  ## Substitute any ".", "._" or "_" for ". ", unless "." is at the end of the string
  oak_spp <- gsub("_", " ", oak_spp)
  ## Replace old label with new one
  oak_tree$tip.label[i] <- oak_spp
}
## Relabel Q. litoralis (accepted name acc. to WCVP v3 is Atuna excelsa subsp. excelsa)
oak_tree$tip.label[220] <- "Q. litoralis (Atuna excelsa)"

## Remove one tip if name is duplicated (even if tips from the same sp. are not sister
## clades). NB, if I remove others instead, results are very similar.
oak_tree$tip.label[duplicated(oak_tree$tip.label)]
grep("Q. arizonica", oak_tree$tip.label)   ## 6, 14
oak_tree$tip.label[6] <- "Q. arizonica2"
grep("Q. laeta", oak_tree$tip.label)       ## 16, 18
oak_tree$tip.label[16] <- "Q. laeta2"
grep("Q. conzattii", oak_tree$tip.label)   ## 114, 116
oak_tree$tip.label[114] <- "Q. conzattii2"
grep("Q. benthamii", oak_tree$tip.label)   ## 131, 132
oak_tree$tip.label[131] <- "Q. benthamii2"
grep("Q. rehderiana", oak_tree$tip.label)  ## 179, 182
oak_tree$tip.label[179] <- "Q. rehderiana2"

oak_tree <- drop.tip(oak_tree,
                     c("Q. arizonica2", "Q. laeta2", "Q. conzattii2",
                       "Q. benthamii2", "Q. rehderiana2",
                       "Q. dentata subsp. yunnanensis",
                       "Q. infectoria subsp. veneris",
                       "Q. ithaburensis subsp. macrolepis",
                       "Q. parvula var. shrevei",
                       "Q. petraea subsp. pinnatiloba",
                       "Q. robur subsp. imeretina"))
oak_tree$tip.label[166] <- "Q. parvula"

## Plot tree and make sure it looks OK
# plot(oak_tree)



## B) READ IN HOST DATA
## Read in host data for oak trees
oak_hosts <- read.table(paste0(pw, "tmp/no_agrilus_quercus_hosts.txt"),
                        sep = "\t")
colnames(oak_hosts) <- c("quercus.sp", "no.agrilus.sp")
oak_hosts <- as.data.frame(lapply(oak_hosts, gsub, pattern = "Quercus",
                                  replacement = "Q.", fixed = TRUE))
str(oak_hosts)

## Create a 0/1 host status data df
oak_host_status <- as.data.frame(matrix(nrow = length(oak_tree$tip.label), ncol = 2))
colnames(oak_host_status) <- c("quercus.sp", "host.status")
oak_host_status$quercus.sp <- oak_tree$tip.label
oak_host_status$host.status <- 0

for (i in 1:nrow(oak_host_status)) {
  if (oak_host_status$quercus.sp[i] %in% oak_hosts$quercus.sp) {
    oak_host_status$host.status[i] <- 1
  }
}


## C) COMPUTE SIGNAL
## Create a comparative data object to match the tree spp. names in the phylo tree to
## those in  the DS.
oak_tree$node.label <- NULL
oak_comparative <- comparative.data(phy = oak_tree,
                                    data = oak_host_status,
                                    ## Name of the column with spp. names
                                    names.col = quercus.sp,
                                    ## Store a variance covariance matrix of your tree
                                    vcv = TRUE,
                                    ## Don't remove spp. w/o data for certain variables
                                    na.omit = FALSE,
                                    ## Shows dropped spp.
                                    warn.dropped = TRUE)

## Now compute strength of the phylogenetic signal (Fritz & Purvis' D metric)
d_oaks <- phylo.d(data = oak_comparative,
                  binvar = host.status, permut = 1000)
d_oaks
## D = 0.55496, which suggests clumping is moderately strong (according to F&P 2010,
## 0.56 is "moderately strong").
## Phylogenetic pattern differs significantly from the Brownian expectation (prob. of
## estimating D under Brownian evolution = 0.001)
## Also differs significantly from 1 (prob of estimating D with random phylogenetic
## structure = 0).

## See this graphically
# svg("results/d_oaks.svg")
plot(d_oaks, lwd=5, cex.axis = 2.5, cex.lab = 2)
axis(side = 2, lwd = 3, lab=F)
axis(side = 1, lwd = 3, lab=F)
box(lwd=3)
# dev.off()




#
##### II)   PLOT HOT CLADES                                                          #####
## A) EXTRACT NO. ENTRIES IN GBIF FOR OAK SPP. IN PHYLOGENY
# ## i) V1: use no. of entries according to GBIF
# ## Read the phylogenetic tree of oaks and extract spp. names.
# oaks <- read.tree("input/quercus_hipp19_singleton_crown_sp_level.tre")
# oaks <- gsub(oaks$tip.label, pattern = "Q.", replacement = "Quercus", fixed = T)
# oaks <- gsub(oaks, pattern = "_", replacement = " ", fixed = T)
# oaks <- gsub(oaks, pattern = " -Atuna excelsa", replacement = "", fixed = T)
# oaks <- gsub(oaks, pattern = "?? ", replacement = "", fixed = T)
#
# ## Extract their GBIF info (takes ca. 1 min to run).
# oak_keys_info <- lapply(oaks, name_lookup, rank = "SPECIES",
#                         higherTaxonKey = "220", status = "Accepted", limit = 1)
#
# ## Create a dataframe to store oak species name and their GBIF key. If there is no key
# ## info for the current sp., save its index in manual_keys. Takes < 1 min to run.
# oak_keys <- as.data.frame(matrix(ncol = 2, nrow = length(oak_keys_info)))
# colnames(oak_keys) <- c("oak", "key")
# manual_keys <- c()
# for (i in 1:length(oak_keys_info)) {
#   if (!is.null(oak_keys_info[[i]][4]$hierarchies)) {
#     sp <- oak_keys_info[[i]][2]$data$species
#     ## So that it doesn't include ""Quercus faginea-mirbeckii" (sp. in "oaks" is
#     ## "Quercus faginea")
#     if (sp %in% oaks) {
#       oak_keys[i,1] <- sp
#       oak_keys[i,2] <- oak_keys_info[[i]][2]$data$key
#     } else {
#       manual_keys <- cbind(manual_keys, i)
#     }
#   } else {
#     manual_keys <- cbind(manual_keys, i)
#   }
# }
# rm(oak_keys_info, i, sp)
#
# ## For all spp. without a key, look them up on GBIF and add them manually (7 oak spp.).
# ## * = name is a syn. according to GBIF but accepted by WCVP
# ## ** == Atuna excelsa subsp. excelsa
# ## *** = loop above recognises another spp. first instead
# oaks[manual_keys]
# no <- grep(oaks[manual_keys][1], x = oaks)    ## "Quercus new species"
# oak_keys[no,1] <- oaks[manual_keys][1]        ## "Quercus new species"
# oak_keys[no,2] <- NA                          ## "Quercus new species"
# no <- grep(oaks[manual_keys][2], x = oaks)    ## "Quercus corrugata"
# oak_keys[no,1] <- oaks[manual_keys][2]        ## "Quercus corrugata"
# oak_keys[no,2] <- 2880406                     ## "Quercus corrugata"*
# no <- grep(oaks[manual_keys][3], x = oaks)    ## "Quercus faginea"
# oak_keys[no,1] <- oaks[manual_keys][3]        ## "Quercus faginea"
# oak_keys[no,2] <- 2881480                     ## "Quercus faginea"***
# no <- grep(oaks[manual_keys][4], x = oaks)    ## "Quercus calophylla"
# oak_keys[no,1] <- oaks[manual_keys][4]        ## "Quercus calophylla"
# oak_keys[no,2] <- 2880691                     ## "Quercus calophylla"*
# no <- grep(oaks[manual_keys][5], x = oaks)    ## "Quercus confertifolia"
# oak_keys[no,1] <- oaks[manual_keys][5]        ## "Quercus confertifolia"
# oak_keys[no,2] <- 2878551                     ## "Quercus confertifolia"*
# no <- grep(oaks[manual_keys][6], x = oaks)    ## "Quercus sartorii"
# oak_keys[no,1] <- oaks[manual_keys][6]        ## "Quercus sartorii"
# oak_keys[no,2] <- 2880720                     ## "Quercus sartorii"*
# no <- grep(oaks[manual_keys][7], x = oaks)    ## "Quercus palustris"
# oak_keys[no,1] <- oaks[manual_keys][7]        ## "Quercus palustris"
# oak_keys[no,2] <- 8313153                     ## "Quercus palustris"
# no <- grep(oaks[manual_keys][8], x = oaks)    ## "Quercus sp. nov."
# oak_keys[no,1] <- oaks[manual_keys][8]        ## "Quercus sp. nov."
# oak_keys[no,2] <- NA                          ## "Quercus sp. nov."
# no <- grep(oaks[manual_keys][9], x = oaks)    ## "Quercus litoralis"
# oak_keys[no,1] <- oaks[manual_keys][9]        ## "Quercus litoralis"
# oak_keys[no,2] <- 10839442                    ## "Quercus litoralis"**
# no <- grep(oaks[manual_keys][10], x = oaks)   ## "Quercus acuta"
# oak_keys[no,1] <- oaks[manual_keys][10]       ## "Quercus acuta"
# oak_keys[no,2] <- 2876863                     ## "Quercus acuta"***
#
# nrow(oak_keys) == length(oaks)
# rm(manual_keys, no)
#
# entries <- c()
# for (i in 1:nrow(oak_keys)) {
#   entries <- c(entries, occ_count(taxonKey = oak_keys$key[i],
#                                   georeferenced = TRUE, from = 1950, type = "count"))
# }
# entries <- data.frame(cbind(oak_keys$oak, entries))
# colnames(entries) <- c("oak.sp", "no.entries")
# rm(oaks, i, oak_keys, oak_keys_info)
# entries$no.entries <- as.numeric(entries$no.entries)
# entries$oak.sp <- gsub("Quercus", "Q.", entries$oak.sp, fixed = T)
# entries$oak.sp <- gsub(pattern = "-", replacement = "",
#                        entries$oak.sp, fixed = T)
#
# entries[entries$oak.sp == "Q. sp.nov.", ]$no.entries <- 0
# entries[entries$oak.sp == "Q. new species", ]$no.entries <- 0
#
#
# ## ii) : use no. of GBIF entries in my cleaned data used for the model analyses
# ## Read table created in V)
# entries <- read.table("tmp/interaction_data_oak_hosts.txt",
#                       header = T, sep = "\t")
# entries <- unique(entries[, c(1,ncol(entries))])
# colnames(entries) <- c("oak.sp", "no.entries")
# entries$oak.sp <- gsub("Quercus", "Q.", entries$oak.sp, fixed = T)
# entries <- rbind(entries, c("Q. sp.nov.", "0"))
# entries <- rbind(entries, c("Q. new species", "0"))
# entries$no.entries <- as.numeric(entries$no.entries)
# str(entries)
#
#
#
# ## B) LOAD  REST OF DATA
# ## Read nexus tree created with phylocom nodesig (first tree: host status info)
# oak_nodesig <- read.nexus("input/quercus_nodesig_result.nex")[[1]]
#
# ## Repace "Quercus" for "Q." (tip labels) and "_" for " "
# oak_nodesig$tip.label <- gsub(pattern = "Quercus_", replacement = "Q. ",
#                               oak_nodesig$tip.label)
# oak_nodesig$tip.label <- gsub(pattern = "_", replacement = " ",
#                               oak_nodesig$tip.label)
#
# ## Read file with info on number of Agrilus species hosted per Quercus species
# hosts <- read.table("tmp/agrilus_quercus_hosts.txt", sep = "\t", header = F)
# hosts <- data.frame(cbind(hosts$, hosts$V1))
# colnames(hosts) <- c("quercus.sp", "no.agrilus")
# hosts$quercus.sp <- as.character(hosts$quercus.sp)
# hosts <- as.data.frame(table(hosts$quercus.sp))
# colnames(hosts) <- c("quercus.sp", "no.agrilus")
# str(hosts)
#
# ## Replace "Quercus" with "Q."
# hosts$quercus.sp <- gsub(pattern = "Quercus", replacement = "Q.", fixed = T,
#                          x = hosts$quercus.sp)
#
#
#
# ## C) PLOT TREE WITH HOT CLADES
# ## Check with Mesquite (in Windows) which clades are hot/cold, and colour them
# ## accordingly with R
# oak_nodesig$node.label <- gsub(pattern = "'", replacement = "",
#                                oak_nodesig$node.label)
#
# oak_nodesig <- groupClade(.data = oak_nodesig, .node = c(
#                           "N246", "N247",
#                           "N293", "N294", "N295", "N296", "N297", "N298",
#                           "N316", "N317",
#                           "N342", "N343", "N344",
#                           "N381",
#                           "N404",
#                           "N419",
#                           "N431",
#                           "N443", "N444", "N445", "N446"))
#
# cols <- c("grey10",
#           brewer.pal(n = 9, name = "Blues")[4:5],
#           brewer.pal(n = 9, name = "YlOrRd")[4:9],
#           brewer.pal(n = 9, name = "YlOrRd")[5:6],
#           brewer.pal(n = 9, name = "Blues")[4:6],
#           brewer.pal(n = 9, name = "YlOrRd")[4],
#           brewer.pal(n = 9, name = "YlOrRd")[4],
#           brewer.pal(n = 9, name = "YlOrRd")[4],
#           brewer.pal(n = 9, name = "YlOrRd")[4],
#           brewer.pal(n = 9, name = "Blues")[4:7])
#
# oak_tree <- ggtree(oak_nodesig,
#                    aes(color = group), size = 1.3, layout = "circular") +
#             scale_color_manual(values = cols) +
#             geom_tiplab(color = "black", offset = 15, size = 1.7, fontface = "italic") +
#             # theme(legend.position = "none") +
#             geom_point2(aes(subset = (label == "N246")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N247")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N293")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N294")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N295")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N296")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N297")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N298")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N316")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N317")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N342")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N343")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N344")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N381")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N404")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N419")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N431")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N443")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N444")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N445")), shape = 20, size = 5) +
#             geom_point2(aes(subset = (label == "N446")), shape = 20, size = 5)
# # plot(oak_tree)
#
# ## NB, re stat = "identity": For geom_bar(), the default (stat = "count") is to count
# ## the rows for each x value. With stat = "identity", you're telling ggplot2 to skip
# ## the aggregation and that you'll provide the y values. From: stackoverflow.com/
# ## questions/59008974/why-is-stat-identity-necessary-in-geom-bar-in-ggplot
# oak_tree2 <- oak_tree + ggnewscale::new_scale_fill() +
#              geom_fruit(data = hosts,
#                         geom = geom_bar,
#                         mapping = aes(y = quercus.sp, x =  no.agrilus, size = 100),
#                         offset = 0.15,
#                         pwidth = 0.1,
#                         stat = "identity",                   ## see above
#                         orientation = "y"                    ## axis orientation
#                         # axis.params = list(axis = "x",     ## add axis text
#                         #               text.angle = -45,    ## text size of axis
#                         #               hjust = 0            ## adjust text horiz. pos.
#                         #               )
#                         # grid.params = list()               ## add grid line
#                         )
#
# ## Save tree
# # svg("results/figures/hot_clades_oaks.svg")
# # plot(oak_tree2)
# # dev.off()
#
#
#
# ## D) ADD INFO ON NO. GBIF ENTRIES
# # oak_tree2$data$label[1:238][!(oak_tree2$data$label[1:238] %in% entries$oak.sp)]
# oak_tree3  <- oak_tree + ggnewscale::new_scale_fill()
# # entries[entries$no.entries >= 100000,]$no.entries <- 100000
#
# oak_tree3 <- oak_tree2 %<+% entries +
#   ggnewscale::new_scale_fill() +
#   geom_tippoint(aes(x = x+4, size = no.entries), # , fill = no.entries
#                  shape = 21, fill = alpha("#86929C", .5), colour = "grey") +
#   scale_size_continuous(range = c(0.1, 10),
#                         limits = c(1, 500000),
#                         breaks = c(1, 1000, 10000,
#                                    100000, 500000)) +
#   guides(color = "none") +
#   theme(legend.position = c(0.9, 0.95))
#   # scale_fill_gradient(low = "lightgoldenrod1",
#   #                     high = "red3",
#   #                     na.value = NA,
#   #                     limits = c(0.1, 8))
# # svg("results/figures/hot_clades_oaks.svg")
# plot(oak_tree3)
# # dev.off()
#
#
#
#
# ## E) ADD INDIVIDUAL AGRILUS INFO
# hosts <- read.table("tmp/agrilus_quercus_hosts.txt", sep = "\t", header = F)
# hosts <- data.frame(cbind(hosts$, hosts$V1))
# colnames(hosts) <- c("quercus.sp", "agrilus.sp")
# hosts$quercus.sp <- gsub(pattern = "Quercus", replacement = "Q.", fixed = T,
#                          x = hosts$quercus.sp)
# hosts$agrilus.sp <- gsub(pattern = "Agrilus", replacement = "A.", fixed = T,
#                          x = hosts$agrilus.sp)
# hosts <- with(hosts, table(quercus.sp, agrilus.sp))
# hosts <- as.data.frame.matrix(hosts)
# rest <- oak_tree3$data$label[1:238][!(oak_tree3$data$label[1:238] %in% rownames(hosts))]
#
# hostnames <- c(rownames(hosts), rest)
# for (i in 1:length(rest)) {
#   hosts <- rbind(hosts, rep(0, ncol(hosts)))
# }
# rownames(hosts) <- hostnames
#
# for (i in 1:ncol(hosts)) {
#   for (j in 1:nrow(hosts)) {
#     if (hosts[j,i] == 1) {
#       hosts[j,i] <- as.character(i)
#     } else {
#       hosts[j,i] <- NA
#     }
#   }
# }
#
# oak_tree4 <- oak_tree3 + ggnewscale::new_scale_fill() +
#              gheatmap(oak_tree4, hosts, colnames_angle = 90,
#                       color = "white", width = 0.5,
#                       offset = 15, font.size = 0.5) +
#              scale_fill_viridis_d(option="D", na.value = "white") +
#              theme(legend.position = 'none')
# plot(oak_tree4)
#
#
#
# #
##### III)  RETRIEVE COORDINATE INFO FROM GBIF                                       #####
# ## A) DOWNLOAD OAK SPP. INFO FOR ALL OAK SPP. PRESENT IN THE PHYLOGENY
# ## Download information for all oak spp. in GBIF (all downloads take < 1 min.):
# ## NB, TaxonKey for Quercus in GBIF == 2877951
# # occ_download(pred("taxonKey", 2877951), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361079-210914110416597')
#
# geo_info <- occ_download_get('0361079-210914110416597', path = "tmp/gbif") %>%
#             occ_download_import(quote = "", path = "tmp/")
#
# # geo_info <- geo_info[, c("gbifID", "species", "scientificName",
# #                          "verbatimScientificName", "countryCode", "occurrenceStatus",
# #                          "individualCount", "decimalLatitude", "decimalLongitude",
# #                          "coordinateUncertaintyInMeters", "coordinatePrecision",
# #                          "year", "taxonKey", "speciesKey", "basisOfRecord",
# #                          "recordNumber", "recordedBy", "issue")]
#
# ## Read the phylogenetic tree of oaks and extract spp. names.
# oak_tree <- read.tree("input/quercus_hipp19_singleton_crown_sp_level.tre")
# plants <- gsub(oak_tree$tip.label, pattern = "Q.", replacement = "Quercus", fixed = T)
# plants <- gsub(plants, pattern = "_", replacement = " ", fixed = T)
# plants <- gsub(plants, pattern = " -Atuna excelsa", replacement = "", fixed = T)
# plants <- gsub(plants, pattern = "?? ", replacement = "", fixed = T)
#
# ## Check no. of oak spp. in the phylogeny (238, which is right) and see if they are
# ## all in the GBIF df:
# length(plants)
# rm(oak_tree)
#
# ## Extract their GBIF info (takes ca. 1 min to run).
# plant_keys_info <- lapply(plants, name_lookup, rank = "SPECIES",
#                           higherTaxonKey = "220", status = "Accepted", limit = 1)
# # plant_keys_info[[1]]
#
# ## Create a dataframe to store plant species name and their GBIF key. If there is no key
# ## info for the current sp., save its index in manual_keys. Takes < 1 min to run.
# plant_keys <- as.data.frame(matrix(ncol = 2, nrow = length(plant_keys_info)))
# colnames(plant_keys) <- c("plant", "key")
# manual_keys <- c()
# for (i in 1:length(plant_keys_info)) {
#   if (!is.null(plant_keys_info[[i]][4]$hierarchies)) {
#     sp <- plant_keys_info[[i]][2]$data$species
#     ## So that it doesn't include ""Quercus faginea-mirbeckii" (sp. in "plants" is
#     ## "Quercus faginea")
#     if (sp %in% plants) {
#       plant_keys[i,1] <- sp
#       plant_keys[i,2] <- plant_keys_info[[i]][2]$data$key
#     } else {
#       manual_keys <- cbind(manual_keys, i)
#     }
#   } else {
#     manual_keys <- cbind(manual_keys, i)
#   }
# }
# rm(plant_keys_info, i, sp)
#
# ## For all spp. without a key, look them up on GBIF and add them manually (7 oak spp.).
# ## * = name is a syn. according to GBIF but accepted by WCVP
# ## ** == Atuna excelsa subsp. excelsa
# ## *** = loop above recognises another spp. first instead
# plants[manual_keys]
# no <- grep(plants[manual_keys][1], x = plants)  ## "Quercus new species"
# plant_keys[no,1] <- plants[manual_keys][1]      ## "Quercus new species"
# plant_keys[no,2] <- NA                          ## "Quercus new species"
# no <- grep(plants[manual_keys][2], x = plants)  ## "Quercus corrugata"
# plant_keys[no,1] <- plants[manual_keys][2]      ## "Quercus corrugata"
# plant_keys[no,2] <- 2880406                     ## "Quercus corrugata"*
# no <- grep(plants[manual_keys][3], x = plants)  ## "Quercus faginea"
# plant_keys[no,1] <- plants[manual_keys][3]      ## "Quercus faginea"
# plant_keys[no,2] <- 2881480                     ## "Quercus faginea"***
# no <- grep(plants[manual_keys][4], x = plants)  ## "Quercus calophylla"
# plant_keys[no,1] <- plants[manual_keys][4]      ## "Quercus calophylla"
# plant_keys[no,2] <- 2880691                     ## "Quercus calophylla"*
# no <- grep(plants[manual_keys][5], x = plants)  ## "Quercus confertifolia"
# plant_keys[no,1] <- plants[manual_keys][5]      ## "Quercus confertifolia"
# plant_keys[no,2] <- 2878551                     ## "Quercus confertifolia"*
# no <- grep(plants[manual_keys][6], x = plants)  ## "Quercus sartorii"
# plant_keys[no,1] <- plants[manual_keys][6]      ## "Quercus sartorii"
# plant_keys[no,2] <- 2880720                     ## "Quercus sartorii"*
# no <- grep(plants[manual_keys][7], x = plants)  ## "Quercus palustris"
# plant_keys[no,1] <- plants[manual_keys][7]      ## "Quercus palustris"
# plant_keys[no,2] <- 8313153                     ## "Quercus palustris"
# no <- grep(plants[manual_keys][8], x = plants)  ## "Quercus sp. nov."
# plant_keys[no,1] <- plants[manual_keys][8]      ## "Quercus sp. nov."
# plant_keys[no,2] <- NA                          ## "Quercus sp. nov."
# no <- grep(plants[manual_keys][9], x = plants)  ## "Quercus litoralis"
# plant_keys[no,1] <- plants[manual_keys][9]      ## "Quercus litoralis"
# plant_keys[no,2] <- 10839442                    ## "Quercus litoralis"**
# no <- grep(plants[manual_keys][10], x = plants) ## "Quercus acuta"
# plant_keys[no,1] <- plants[manual_keys][10]     ## "Quercus acuta"
# plant_keys[no,2] <- 2876863                     ## "Quercus acuta"***
#
# nrow(plant_keys) == length(plants)
# rm(manual_keys, no)
#
# ## Subset geo_info so as to only keep those in the phylogeny
# geo_info <- subset(geo_info, speciesKey %in% plant_keys$key & !is.na(speciesKey))
# # unique(geo_info[is.na(geo_info$speciesKey),]$species)
# geo_info <- subset(geo_info, speciesKey %in% plant_keys$key)
# setdiff(geo_info$speciesKey, plant_keys$key)
# setdiff(geo_info$species, plant_keys$plant)
# missing_keys <- setdiff(unique(plant_keys$key), unique(geo_info$speciesKey))
# missing_sp <- setdiff(unique(plant_keys$plant), unique(geo_info$species))
#
# ## For those species not in geo_info, download their info manually:
# ## Quercus new sp.
# missing_sp[1]
# missing_keys[1]
#
# ## Quercus corrugata: - Don't run! It's already in dataset as part of some Q. lancifolia
# ## entries
# # missing_sp[2]
# # ## "taxonKey" (instead of "speciesKey") because it is a syn. according to GBIF
# # # occ_download(pred("taxonKey", missing_keys[2]), format = "SIMPLE_CSV",
# # #              user = "elvhergut", email = "elvira.herg@gmail.com",
# # #              pwd = "agrilusoakgbif")
# # # occ_download_wait('0361453-210914110416597')
# # add_plant <- occ_download_get('0361453-210914110416597', path = "tmp/gbif") %>%
# #              occ_download_import(quote = "", path = "tmp/")
# # # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# # #                            "verbatimScientificName", "countryCode",
# # #                            "occurrenceStatus", "individualCount", "decimalLatitude",
# # #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# # #                            "coordinatePrecision", "year", "taxonKey", "speciesKey",
# # #                            "basisOfRecord", "recordNumber", "recordedBy", "issue")]
# # geo_info <- rbind(geo_info, add_plant)
# # rm(add_plant)
#
# ## Quercus sagrana: no occurrences in GBIF
# missing_sp[3]
# missing_keys[3]
#
# ## Quercus calophylla
# missing_sp[4]
# ## "taxonKey" because it is a syn. according to GBIF
# # occ_download(pred("taxonKey", missing_keys[4]), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361458-210914110416597')
# add_plant <- occ_download_get('0361458-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName", "countryCode",
# #                            "occurrenceStatus", "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey", "speciesKey",
# #                            "basisOfRecord", "recordNumber", "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Quercus confertifolia: Don't run! It's already in dataset as part of some Q.
# ## crassipes entries
# # missing_sp[5]
# # ## "taxonKey" because it is a syn. according to GBIF
# # # occ_download(pred("taxonKey", missing_keys[5]), format = "SIMPLE_CSV",
# # #              user = "elvhergut", email = "elvira.herg@gmail.com",
# # #              pwd = "agrilusoakgbif")
# # # occ_download_wait('0361447-210914110416597')
# # add_plant <- occ_download_get('0361547-210914110416597', path = "tmp/gbif") %>%
# #              occ_download_import(quote = "", path = "tmp/")
# # # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# # #                            "verbatimScientificName",
# # #                            "countryCode", "occurrenceStatus",
# # #                            "individualCount", "decimalLatitude",
# # #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# # #                            "coordinatePrecision", "year", "taxonKey",
# # #                            "speciesKey", "basisOfRecord", "recordNumber",
# # #                            "recordedBy", "issue")]
# # geo_info <- rbind(geo_info, add_plant)
# # rm(add_plant)
#
# ## Quercus sartorii
# ## "taxonKey" because it is a syn. according to GBIF
# missing_sp[6]
# # occ_download(pred("taxonKey", missing_keys[6]), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361555-210914110416597')
# add_plant <- occ_download_get('0361555-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Quercus sp.nov.
# missing_sp[7]
# missing_keys[1]
#
# ## Quercus litoralis (NB, key is actually == 10839442, as this is the key for
# ## Atuna excelsa subsp. excelsa, which is the accepted name in GBIF)
# ## "taxonKey" because it is a subsp. according to GBIF
# missing_sp[8]
# # occ_download(pred("taxonKey", 10839442), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361560-210914110416597')
# add_plant <- occ_download_get('0361560-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName", "countryCode",
# #                            "occurrenceStatus", "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey", "speciesKey",
# #                            "basisOfRecord", "recordNumber", "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
#
#
# ## B) DOWNLOAD PLANT SPP. INFO FOR ALL NON-OAK SPP. THAT HOST AGRILUS THAT USE OAKS
# ## Read in file with info on other plants that plant Agrilus spp. used by oaks
# ## NB, I have only retained info for those plants for which there is info to sp. level
# ## (for some there is only info to genus level).
# plants <- read.table("tmp/non_oak_hosts.tsv", header = T, sep = "\t")
#
# ## Check no. of non-oak plant spp. used by Agrilus spp. that exploit oaks (19).
# length(unique(plants$hosts))
#
# ## Extract GBIF info (takes ca. 1 min to run).
# plant_keys_info <- lapply(unique(plants$hosts), name_lookup, rank = "SPECIES",
#                           higherTaxonKey = "220", status = "Accepted", limit = 1)
# # plant_keys_info[[1]]
#
# ## Create a dataframe to store plant species name and their GBIF key. If there is no key
# ## info for the current sp., save its index in manual_keys. Takes < 1 min to run.
# plant_keys <- as.data.frame(matrix(ncol = 2, nrow = length(plant_keys_info)))
# colnames(plant_keys) <- c("plant", "key")
# manual_keys <- c()
# for (i in 1:length(plant_keys_info)) {
#   if (!is.null(plant_keys_info[[i]][4]$hierarchies)) {
#     plant_keys[i,1] <- plant_keys_info[[i]][2]$data$species
#     plant_keys[i,2] <- plant_keys_info[[i]][2]$data$key
#   } else {
#     manual_keys <- cbind(manual_keys, i)
#   }
# }
# manual_keys
# rm(plant_keys_info, i)
#
# ## Download the info for these extra plant spp. and add them to geo_info
# ## Notholithocarpus densiflorus
# plant_keys[1,]
# # occ_download(pred("taxonKey", plant_keys[1,2]), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361580-210914110416597')
# add_plant <- occ_download_get('0361580-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Betula pendula
# plant_keys[2,]
# # occ_download(pred("taxonKey", plant_keys[2,2]), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361582-210914110416597')
# add_plant <- occ_download_get('0361582-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Castanea sativa
# plant_keys[3,]
# # occ_download(pred("taxonKey", plant_keys[3,2]), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361584-210914110416597')
# add_plant <- occ_download_get('0361584-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Corylus avellana
# plant_keys[4,]
# # occ_download(pred("taxonKey", c(plant_keys[4,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361600-210914110416597')
# add_plant <- occ_download_get('0361600-210914110416597', path = "\tmp") %>%
#              occ_download_import(quote = "", path = "\tmp")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Fagus sylvatica
# plant_keys[5,]
# # occ_download(pred("taxonKey", c(plant_keys[5,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361610-210914110416597')
# add_plant <- occ_download_get('0361610-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Ostrya carpinifolia
# plant_keys[6,]
# # occ_download(pred("taxonKey", c(plant_keys[6,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361624-210914110416597')
# add_plant <- occ_download_get('0361624-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Castanea dentata
# plant_keys[7,]
# # occ_download(pred("taxonKey", c(plant_keys[7,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361625-210914110416597')
# add_plant <- occ_download_get('0361625-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Fagus grandifolia
# plant_keys[8,]
# # occ_download(pred("taxonKey", c(plant_keys[8,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361626-210914110416597')
# add_plant <- occ_download_get('0361626-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Salix nigra
# plant_keys[9,]
# # occ_download(pred("taxonKey", c(plant_keys[9,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361628-210914110416597')
# add_plant <- occ_download_get('0361628-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Diospyros virginiana
# plant_keys[10,]
# # occ_download(pred("taxonKey", c(plant_keys[10,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361642-210914110416597')
# add_plant <- occ_download_get('0361642-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Gleditsia triacanthos
# plant_keys[11,]
# # occ_download(pred("taxonKey", c(plant_keys[11,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361643-210914110416597')
# add_plant <- occ_download_get('0361643-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Euonymus europaeus
# plant_keys[12,]
# # occ_download(pred("taxonKey", c(plant_keys[12,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361646-210914110416597')
# add_plant <- occ_download_get('0361646-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Carpinus betulus
# plant_keys[13,]
# # occ_download(pred("taxonKey", c(plant_keys[13,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361649-210914110416597')
# add_plant <- occ_download_get('0361649-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Ficus carica
# plant_keys[14,]
# # occ_download(pred("taxonKey", c(plant_keys[14,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361653-210914110416597')
# add_plant <- occ_download_get('0361653-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
#
# ## Aesculus glabra
# plant_keys[15,]
# # occ_download(pred("taxonKey", c(plant_keys[15,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361656-210914110416597')
# add_plant <- occ_download_get('0361656-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Aesculus pavia
# plant_keys[16,]
# # occ_download(pred("taxonKey", c(plant_keys[16,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361661-210914110416597')
# add_plant <- occ_download_get('0361661-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Betula lenta
# plant_keys[17,]
# # occ_download(pred("taxonKey", c(plant_keys[17,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361664-210914110416597')
# add_plant <- occ_download_get('0361664-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Carpinus caroliniana
# plant_keys[18,]
# # occ_download(pred("taxonKey", c(plant_keys[18,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361668-210914110416597')
# add_plant <- occ_download_get('0361668-210914110416597', path = "tmp/gbif") %>%
#              occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# ## Ostrya virginiana
# plant_keys[19,]
# # occ_download(pred("taxonKey", c(plant_keys[19,2])), format = "SIMPLE_CSV",
# #              user = "elvhergut", email = "elvira.herg@gmail.com",
# #              pwd = "agrilusoakgbif")
# # occ_download_wait('0361677-210914110416597')
# add_plant <- occ_download_get('0361677-210914110416597', path = "tmp/gbif") %>%
#   occ_download_import(quote = "", path = "tmp/")
# # add_plant <- add_plant[, c("gbifID", "species", "scientificName",
# #                            "verbatimScientificName",
# #                            "countryCode", "occurrenceStatus",
# #                            "individualCount", "decimalLatitude",
# #                            "decimalLongitude", "coordinateUncertaintyInMeters",
# #                            "coordinatePrecision", "year", "taxonKey",
# #                            "speciesKey", "basisOfRecord", "recordNumber",
# #                            "recordedBy", "issue")]
# geo_info <- rbind(geo_info, add_plant)
# rm(add_plant)
#
# # write.table(geo_info, "tmp/gbif_plants.csv", quote = F, row.names = F,
# #             col.names = T, sep = "\t")
#
# ## NB, check for absolute duplicates with:
# # sort tmp/gbif_plants.csv | uniq -c | sort -nr | grep '^\t* *2' | cut -f2
#
#
#
# ## C) CLEAN PLANT GEO DATA FROM GBIF
# ## See https://cran.r-project.org/web/packages/CoordinateCleaner/vignettes/
# ## Cleaning_GBIF_data_with_CoordinateCleaner.html
# ## i) Extract plant geographical data and remove any entries in discordance with WCVP
# ## Extract data (created in A) and B))
# geo_info <- read.table("tmp/gbif_plants.csv", header = T, sep = "\t",
#                        quote = "", comment.char = '')
#
# ## Remove any entries that are in gbif_notwcvp_entries.txt, i.e., those for which the
# ## name under which the sp. was originally reported ("scientificName") is regarded as a
# ## syn. for the current sp. ("species") by GBIF but not WCVP. Takes ca. 30 min to run.
# ## NB, conservative approach, some of these species might actually be synonyms... (but
# ## then I think they all have very few occurrence records so not likely to make much
# ## of an impact anyways).
# gbif_notwcvp <- read.table("tmp/gbif_notwcvp_auth.txt", header = F, sep = "\t")
# head(gbif_notwcvp)
# for (i in 1:nrow(geo_info)) {
#   if (i%%250000 == 0) {print(i)}
#   sp <- geo_info$verbatimScientificName[i]
#   sp <- gsub(pattern = "??|-|subsp\\. | var\\.|f\\.| de | ex ", "",  sp)
#   # sp <- gsub(pattern = "([A-Z][a-z]* [a-z][a-z]* *[a-z]*[a-z]) *.*", "\\1",  sp)
#   sp <- gsub(pattern = "  *", " ",  sp)
#   if (sp %in% gbif_notwcvp) {
#     # print(sp)
#     geo_info$scientificName[i] <- NA
#   }
# }
# tail(geo_info[is.na(geo_info$scientificName), c(10,13,14)], 2)
# geo_info <- geo_info[!is.na(geo_info$scientificName),]
#
#
# ## ii) Clean geo_info: rename spp. that I am interested in but that are reported as
# ## sth else (i.e., the sp. I am interested in is a syn. of sth else acc. to GBIF,
# ## but it's actually an accepted name acc. to WCVP)*.
# ## *I removed these spp. from tmp/gbif_notwcvp_auth.txt intentionally so that they
# ## wouldn't get removed in the previous step.
# ## Get name of all plants_all that should be in my dataset
# plants_all <- read.tree("input/quercus_hipp19_singleton_crown_sp_level.tre")
# plants_all <- gsub(plants_all$tip.label, pattern = "Q.", replacement = "Quercus",
#                    fixed = T)
# plants_all <- gsub(plants_all, pattern = "_", replacement = " ",
#                    fixed = T)
# plants_all <- gsub(plants_all, pattern = " -Atuna excelsa", replacement = "",
#                    fixed = T)
#
# plants_all <- gsub(plants_all, pattern = "?? ", replacement = "", fixed = T)
# plants_all <- c(plants_all, unique(read.table("tmp/non_oak_hosts.tsv",
#                                               header = T, sep = "\t")$hosts))
# length(plants_all) == (238 + 19)  ## 238 Quercus spp. + 19 extra hosts
#
# ## Are they all in geo_info?
# plants_all[!(plants_all %in% unique(geo_info$species))]
#
# ## Re-name species that are reported as something else
# ## NB, Q. new species, Q. sagrana and Q. sp.nov. have no entries in GBIF.
# ## Q. corrugata is reported as Q. lancifolia (which is also in my dataset):
# # nrow(subset(geo_info, species == "Quercus lancifolia" &
# #               scientificName != "Quercus corrugata Hook."))
# geo_info$species[geo_info$species == "Quercus lancifolia" & geo_info$scientificName ==
#                  "Quercus corrugata Hook."] <- "Quercus corrugata"
#
# ## Q. calophylla is reported as Q. candicans:
# # unique(subset(geo_info, taxonKey == 2880691)$species)
# geo_info$species <- gsub(pattern = "Quercus candicans",
#                          replacement = "Quercus calophylla", x = geo_info$species)
#
# ## Q. confertifolia is reported as Q. crassipes (which is also in my dataset):
# # unique(subset(geo_info, taxonKey == 2878551)$species)
# # nrow(subset(geo_info, species == "Quercus crassipes" &
# #             scientificName != "Quercus confertifolia Bonpl."))
# geo_info$species[geo_info$species == "Quercus crassipes" & geo_info$scientificName ==
#                  "Quercus confertifolia Bonpl."] <- "Quercus confertifolia"
#
# ## Q. sartorii is reported as Q. xalapensis:
# # unique(subset(geo_info, taxonKey == 2880720)$species)
# geo_info$species <- gsub(pattern = "Quercus xalapensis",
#                          replacement = "Quercus sartorii", x = geo_info$species)
#
# ## Q. litoralis is reported as Atuna excelsa:
# geo_info$species <- gsub(pattern = "Atuna excelsa",
#                          replacement = "Quercus litoralis", x = geo_info$species)
#
# ## Are they all now in geo_info?
# ## Yes, except for Q. sagrana (no records)
# plants_all[!(plants_all %in% unique(geo_info$species))]
#
# # ## Write file:
# # write.table(geo_info, "tmp/gbif_plants_nameschecked.tsv", sep = "\t",
# #             col.names = T, row.names = F, quote = F)
#
#
# ## iii) Clean geo_info: only keep records of interest
# geo_info_o <- read.table("tmp/gbif_plants_nameschecked.tsv", header = T, sep = "\t",
#                          quote = "", comment.char = '')
# length(unique(geo_info_o$species))
# nrow(geo_info_o)
#
# ## Only keep entries with coordinate info AND that present a coordinate uncertainty
# ## of <= 5 km
# # length(geo_info$coordinateUncertaintyInMeters[
# #        is.na(geo_info$coordinateUncertaintyInMeters)])
# # hist(geo_info$coordinateUncertaintyInMeters/1000, breaks = 1000,
# #      xlim = c(0, 200), ylim = c(0,10000))
# geo_info <- subset(geo_info_o, !is.na(decimalLatitude) & !is.na(decimalLongitude))
# geo_info <- subset(geo_info, is.na(coordinateUncertaintyInMeters) |
#                      coordinateUncertaintyInMeters <= 5000)
# length(unique(geo_info$species))
# nrow(geo_info)
#
# ## Remove entries with high individualCounts
# ## See https://www.gbif.org/data-quality-requirements-occurrences#dcCount
# # hist(geo_info$individualCount[!is.na(geo_info$individualCount)])
# # length(geo_info$individualCount[geo_info$individualCount >= 100 &
# #        !is.na(geo_info$individualCount)])
# geo_info <- subset(geo_info, is.na(individualCount) | individualCount <= 10)
# geo_info <- subset(geo_info, is.na(individualCount) | individualCount > 0)
# length(unique(geo_info$species))
#
# ## Remove entries with occurrenceStatus == "ABSENT"
# geo_info <- subset(geo_info, occurrenceStatus == "PRESENT")
# nrow(geo_info)
#
# ## Remove any entries that come from fossil records
# geo_info <- subset(geo_info, basisOfRecord != "FOSSIL_SPECIMEN")
# geo_info <- subset(geo_info, issue != "FOSSIL_SPECIMEN")
# nrow(geo_info)
#
# ## Remove any entries where the coordinate info is unreliable
# ## To see all issues in data, run (bash):
# ## cut -f18 tmp/gbif_plants.csv | sort | uniq | sed -'s/;/\n/g' | sort | uniq
# ## NB, for issues see https://data-blog.gbif.org/post/issues-and-flags/
# issues <- c("COORDINATE_INVALID", "COORDINATE_OUT_OF_RANGE",
#             "COORDINATE_REPROJECTION_SUSPICIOUS",
#             "COORDINATE_UNCERTAINTY_METERS_INVALID",
#             "COUNTRY_COORDINATE_MISMATCH", "PRESUMED_SWAPPED_COORDINATE",
#             "PRESUMED_NEGATED_LONGITUDE", "PRESUMED_NEGATED_LATITUDE",
#             "ZERO_COORDINATE")
# geo_info <- subset(geo_info, !grepl(pattern = paste(issues, collapse =  "|"), issue))
# rm(issues)
# nrow(geo_info)
#
# ## Are they all still in geo_info?
# ## Yes, except for Q. sagrana (no records) and Q. yiwuensis (no coordinate info)
# plants_all[!(plants_all %in% unique(geo_info$species))]
# length(unique(geo_info$species))
#
# ## Remove records pre 1950
# geo_info <- subset(geo_info, !is.na(year) & year >= 1950)
# nrow(geo_info)
#
# ## Are they all still in geo_info?
# ## Yes, except for Q. sagrana (no records) and Q. yiwuensis (no coordinate info)
# plants_all[!(plants_all %in% unique(geo_info$species))]
# length(unique(geo_info$species))
#
# ## Save data
# rm(geo_info_o)
# # save.image("tmp/geo_info_clean.RData")
#
#
# ## iv) Clean geo_info: use CoordinateCleaner
# ## Check https://mran.microsoft.com/snapshot/2019-02-10/web/packages/CoordinateCleaner/
# ## vignettes/Tutorial_Cleaning_GBIF_data_with_CoordinateCleaner.html
# load("tmp/geo_info_clean.RData")
#
# ## First, convert country code from ISO2c to ISO3c
# geo_info$countryCode <-  countrycode::countrycode(geo_info$countryCode,
#                                                   origin =  'iso2c',
#                                                   destination = 'iso3c')
#
# ## And plot data to get an overview
# wm <- borders("world", colour="gray50", fill="gray50")
# ggplot() + coord_fixed()+ wm +
#            geom_point(data = geo_info, aes(x = decimalLongitude, y = decimalLatitude),
#                       colour = "darkred", size = 0.5)+
#            theme_bw()
#
# ## Then, flag any problematic entries.
# ## Issues I am flagging:
# ## *CAPITALS: tests a radius around adm-0 capitals (radius = capitals_rad, default =
# ## 10,000 m, here = 1,000 m).
# ## *CENTROIDS: tests a radius around country centroids (both country and provinces, check
# ## centroids_detail) (radius = centroids_rad, default = 1,000 m, here = 1,000 m).
# ## *COUNTRIES: tests if coordinates are from the country in the country column. SWITCHED
# ## OFF as it flags up many records that are very close to country borders.
# ## *DUPLICATES: tests for duplicate records (identical coordinates).
# ## *EQUAL: tests for equal absolute longitude and latitude.
# ## *GBIF: tests a one-degree radius around the GBIF headquarters in Copenhagen, Denmark.
# ## *INSTITUTIONS: tests a radius around known biodiversity institutions from institutions
# ## (radius = inst_rad, default = 100 m, here = 1,000 m).
# ## *OUTLIERS: tests each species for outlier records. Depending on outliers_method, it
# ## either flags records that are a min. dist away from all other records (by default,
# ## outliers_td = 1000 km), or records outside a multiple of the IQR of the min. distances
# ## to the next record (default; by default, outliers_mtp = 5). The minimum number of
# ## records in a dataset to run this test is given by outliers_size (default = 7, here
# ## = 50).
# ## *SEAS: tests if coordinates fall into the ocean.
# ## *ZEROS: tests for plain zeros, equal latitude and longitude and a radius in degrees
# ## around the point 0/0. The radius is zeros_rad (default = 0.5, here = 1).
# source("input/coordinate_cleaner_edited.R")
#
# ## Without rownames(geo_info) <- 1:nrow(geo_info) it throws and error message
# rownames(geo_info) <- 1:nrow(geo_info)
# flags <- clean_coordinates2(x = geo_info, lon = "decimalLongitude",
#                             lat = "decimalLatitude",
#                             # countries = "countryCode",
#                             species = "species",
#                             capitals_rad = 1000,
#                             centroids_rad = 1000,
#                             outliers_size = 50,
#                             inst_rad = 1000,
#                             zeros_rad = 1,
#                             tests = c("capitals", "centroids", # "countries",
#                                      "duplicates", "equal","outliers", "gbif",
#                                       "institutions", "zeros", "seas"))
#
# ## Check summary: 52% of data is flagged (48% of all entries are duplicates!!)
# summary(flags)
# # head(flags[!flags$.con,], 2)
# # subset(geo_info, decimalLatitude == 34.73667)[, c("species", "decimalLongitude",
# #                                                   "decimalLatitude", "year",
# #                                                   "coordinatePrecision", "recordedBy")]
#
# plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")
#
# ## Remove any flagged entries
# geo_info_clean <- geo_info[flags$.summary,]
# # # ## Are they all still in geo_info?
# # ## Yes, except for Q. sagrana (no records)
# # plants_all[!(plants_all %in% unique(geo_info_clean$species))]
#
# ## Save data
# # save.image("tmp/geo_info_clean.RData")
#
#
# ## v) Only keep relevant fields
# load("tmp/geo_info_clean.RData")
# geo_info <- geo_info_clean
# rm(flags, geo_info_clean)
# geo_info <- geo_info[, c("species", "decimalLatitude", "decimalLongitude")]
# colnames(geo_info) <- c("plant.sp", "lat", "lon")
# length(unique(geo_info$plant.sp))
# nrow(geo_info)
#
#
# ## vi) Clean geo_info: subsample so that no sp. has > 100k entries (otherwise it
# ## doesn't run)
# sort(table(geo_info$plant.sp)[table(geo_info$plant.sp) >= 100000])
# spp <- dimnames(sort(table(geo_info$plant.sp)[table(geo_info$plant.sp) >= 100000]))[[1]]
# add_all <- subset(geo_info, plant.sp == spp[1])
# add <- add_all[sample(nrow(add_all), 100000), ]
# add_all <- subset(geo_info, plant.sp == spp[2])
# add <- rbind(add, add_all[sample(nrow(add_all), 100000), ])
# add_all <- subset(geo_info, plant.sp == spp[3])
# add <- rbind(add, add_all[sample(nrow(add_all), 100000), ])
# add_all <- subset(geo_info, plant.sp == spp[4])
# add <- rbind(add, add_all[sample(nrow(add_all), 100000), ])
# add_all <- subset(geo_info, plant.sp == spp[5])
# add <- rbind(add, add_all[sample(nrow(add_all), 100000), ])
# add_all <- subset(geo_info, plant.sp == spp[6])
# add <- rbind(add, add_all[sample(nrow(add_all), 100000), ])
#
# geo_info <- subset(geo_info, !(plant.sp %in% spp))
# geo_info <- rbind(geo_info, add)
# length(unique(geo_info$plant.sp))
# nrow(geo_info)
# rm(add, add_all, spp)
#
# ## vii) Clean geo_info: add a "fake" manual entry for the spp. missing (country centroid
# ## according to CoordinateCleaner).
# ## Q. sagrana
# ## Present in Cuba: powo.science.kew.org/taxon/urn:lsid:ipni.org:names:216368-2
# add <- data.frame(cbind("Quercus sagrana",
#                   subset(countryref, iso3 == "CUB" &
#                          type == "country")[1, c("centroid.lat", "centroid.lon")]))
# colnames(add) <- colnames(geo_info)
# geo_info <- rbind(geo_info, add)
#
# ## Q. yiwuensis
# ## Present in in SC China: powo.science.kew.org/taxon/urn:lsid:ipni.org:names:360253-1
# add <- data.frame(cbind("Quercus yiwuensis",
#                   subset(countryref, iso3 == "CHN" &
#                          name == "Hunan")[1, c("centroid.lat", "centroid.lon")]))
# colnames(add) <- colnames(geo_info)
# geo_info <- rbind(geo_info, add)
# plants_all[!(plants_all %in% unique(geo_info$plant.sp))]
# nrow(geo_info)
# rm(add)
#
# ## Write table
# # write.table(geo_info, "tmp/gbif_plants_clean.tsv", quote = F, sep = "\t",
# #             col.names = T, row.names = F)
#
#
#
# #
##### IV)   CREATE A DF WITH INFO ON OAK - INTERACTION DATA                          #####
# ## A) READ PLANT GEO DATA FROM GBIF
# ## Read in file
# geo_info <- read.table(paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/",
#                               "01_original_models/tmp/gbif_plants_clean.tsv"),
#                        header = T, sep = "\t")
# geo_info <- geo_info[, c("plant.sp", "lon", "lat")]
# str(geo_info)
#
# ## Check that all plant species that should be there are present
# plants <- read.tree(paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/",
#                                "01_original_models/input/",
#                                "quercus_hipp19_singleton_crown_sp_level.tre"))
# plants <- gsub(plants$tip.label, pattern = "Q.", replacement = "Quercus",
#                    fixed = T)
# plants <- gsub(plants, pattern = "_", replacement = " ",
#                    fixed = T)
# plants <- gsub(plants, pattern = " -Atuna excelsa", replacement = "",
#                    fixed = T)
#
# plants <- gsub(plants, pattern = "?? ", replacement = "", fixed = T)
# plants <- c(plants,
#                 unique(read.table(paste0("/data/SBCS-NicholsLab/elvirahg/",
#                                          "oak_analyses/01_original_models/",
#                                          "tmp/non_oak_hosts.tsv"),
#                                   header = T, sep = "\t")$hosts))
# length(plants) == (238 + 19)  ## 238 Quercus spp. + 19 extra hosts
#
#
#
# ## B) EXTRACT INFO ON AGRILUS SPP. THAT EXPLOIT OAKS AND THEIR HOSTS
# agrilus_hosts <- read.table(paste0("/data/SBCS-NicholsLab/elvirahg/",
#                                    "oak_analyses/01_original_models/",
#                                    "input/agrilus_quercus_hosts.txt"),
#                             header = F, sep = "\t")
# colnames(agrilus_hosts) <- c("agrilus", "hosts")
# agrilus_hosts <- rbind(agrilus_hosts,
#                        read.table(paste0("/data/SBCS-NicholsLab/elvirahg/",
#                                          "oak_analyses/01_original_models/",
#                                          "tmp/non_oak_hosts.tsv"),
#                                   header = T, sep = "\t"))
# agrilus_hosts <- agrilus_hosts[,c(2,1)]
# colnames(agrilus_hosts) <- c("plant.sp", "agrilus.sp")
# agrilus_hosts$agrilus.sp <- gsub("Agrilus", "A.", agrilus_hosts$agrilus.sp)
#
# ## Remove any hosts that are not in the phylogeny (only Q. gambelii) from agrilus_hosts
# agrilus_hosts <- agrilus_hosts[!grepl("gambelii", agrilus_hosts$plant.sp),]
#
#
# ## IMPORTANT (!): I have also decided to remove any reported hosts that are non-native
# ## for those Agrilus species that have been introduced into new areas.
# ## A. auroguttatus: introduced from AZ (USA) to CA (USA). In CA, reported novel larval
# ## hosts are Q. kelloggii, Q. chrysolepis & Q. agrifolia.
# ## A. bilineatus: introduced from N America to Turkey (no new reported hosts so far).
# ## A. sulcicollis: introduced from Europe to N America, where Q. macrocarpa (in my DS)
# ## is a new reported larval host.
# head(agrilus_hosts[agrilus_hosts$plant.sp == "Quercus macrocarpa" &
#                    agrilus_hosts$agrilus.sp == "A. sulcicollis",])
#
# agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Quercus agrifolia"
#                                & agrilus_hosts$agrilus.sp == "A. auroguttatus")), ]
# agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Quercus macrocarpa"
#                                & agrilus_hosts$agrilus.sp == "A. sulcicollis")), ]
#
# ## IMPORTANT (!): Let's also remove any larval hosts according to J&P(2014) (i.e.,
# ## there's at least one reference with '!') but with a  confidence index < 3.
# ## A. graminis: Corylus avellana, Ostrya carpinifolia, Euonymus europaeus and Castanea
# ## sativa.
# ## A. obscuricollis: Ficus carica
# ## A. relegatus (reported as A. dualis)	Quercus robur
# agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Ficus carica"
#                                  & agrilus_hosts$agrilus.sp == "A. obscuricollis")), ]
#
# agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Quercus robur"
#                                  & agrilus_hosts$agrilus.sp == "A. relegatus")), ]
#
# agrilus_hosts <- agrilus_hosts[-(which((agrilus_hosts$plant.sp == "Corylus avellana" |
#                                  agrilus_hosts$plant.sp == "Ostrya carpinifolia" |
#                                  agrilus_hosts$plant.sp == "Euonymus europaeus" |
#                                  agrilus_hosts$plant.sp == "Castanea sativa") &
#                                  agrilus_hosts$agrilus.sp == "A. graminis")), ]
#
#
#
# ## C) CREATE INTERACTION DF
# ## Extract all plant and Agrilus spp.
# oaks <- grep("Quercus", unique(geo_info$plant.sp), value = TRUE)
# agrilus <- unique(agrilus_hosts$agrilus.sp)
#
# ## Check which oak spp from the phylogeny & which non-oak hosts are missing in "oaks"
# ## (they won't be present in the final dataset).
# setdiff(plants, oaks)
# rm(plants)
#
# ## Create a DF with all Agrilus sp. - oak sp. combinations (interaction_data)
# interaction_data <- expand.grid(oaks, agrilus)
# colnames(interaction_data) <- c("plant.sp", "agrilus.sp")
# nrow(interaction_data)
#
# ## Append info on Agrilus hosts (i.e., oak spp. and the beetle sp. they host) at the
# ## end of interaction_data
# interaction_data <- rbind(interaction_data,
#                           agrilus_hosts[grep("Quercus", agrilus_hosts$plant.sp),])
#
# ## Make it so that hosts = 1, non-hosts = 0 (i.e., count number of entries for each
# ## plant - Agrilus pair and if == 1 => 0, 2 => 1)
# interaction_data <- ddply(interaction_data, .(plant.sp, agrilus.sp),  nrow)
# colnames(interaction_data) <- c("quercus.sp", "agrilus.sp", "interaction")
# interaction_data$interaction <- interaction_data$interaction - 1
# nrow(interaction_data) == length(agrilus)*length(oaks)
#
# rm(agrilus, oaks)
#
#
#
# #
##### V)    CREATE A DF WITH INFO ON OAK - ADD GEO OVERLAP DATA                      #####
# ## A) ADD DISTANCE INFORMATION
# ## For each plant sp. - Agrilus sp. pair, find the mean min. distance to the nearest
# ## host (from a different plant sp.)  of the current Agrilus sp.
#
# ## NB: "geodist() accepts only one or two primary arguments, which must be rectagular
# ## objects with unambiguously labelled longitude and latitude columns (i.e., some
# ## variant of lon/lat, or x/y)." Input(s) "can be in arbitrary rectangular format."
# ## See cran.r-project.org/web/packages/geodist/vignettes/geodist.html
#
# ## NB: Geodesic distance = "a curve representing the shortest path (arc) between two
# ## points in a surface", tse3.mm.bing.net/th?id=OIP.nAhtZINY0zbx_wX-xqzdJwHaEr&pid=Api
# ## measure = "haversine" ("geodesic" is too computationally demanding)
#
# ## Add future columns to interaction_data
# interaction_data$min.dist.mean <- rep(NA, nrow(interaction_data))
# interaction_data$min.dist.mean.norm <- rep(NA, nrow(interaction_data))
# interaction_data$min.dist.mean.log <- rep(NA, nrow(interaction_data))
# interaction_data$min.dist.median <- rep(NA, nrow(interaction_data))
# interaction_data$min.dist.median.norm <- rep(NA, nrow(interaction_data))
# interaction_data$min.dist.median.log <- rep(NA, nrow(interaction_data))
# str(interaction_data)
#
# ## Loop is not the most efficient: takes ca. 24 h to run (ask for 5 GB/iteration).
# ## Most iterations take < 1 h but some (30 or so??) take longer.
# ## Need to run this inside a jobscript with -t 1-236.
# oak <- as.numeric(commandArgs(trailingOnly=TRUE))
# oak <- as.character(unique(interaction_data$quercus.sp)[oak])
# first <- grep(oak, interaction_data$quercus.sp)[1]
# print(oak)
# print(nrow(interaction_data))
#
# ## For each oak - Agrilus pair (one oak at a time)
# options(warn = 2)
# for(i in first:(first+31)) {  ## 1:nrow(interaction_data)
#   ## Get current plant sp
#   oak <- interaction_data$quercus.sp[i]
#
#   ## Get current agrilus sp.
#   agrilus <- interaction_data$agrilus.sp[i]
#
#   ## Find all host spp. (that are not the current oak sp., if it happens to be a host)
#   hosts <- subset(agrilus_hosts, agrilus.sp == agrilus &
#                     plant.sp != oak)$plant.sp
#   print(paste0(i, ": ", oak, " - ", agrilus), )
#
#   ## If there are no hosts for the current Agrilus sp. other than the current oak sp.,
#   ## let's set min.dist.mean to 30k km (Earth's circumference = ca. 40k km)
#   if (length(hosts) == 0) {
#     min.dist.mean <- 3e4
#     min.dist.median <- 3e4
#   }
#   ## Otherwise
#   else {
#     ## Extract geo info
#     ## i) For all individuals from the current oak sp.
#     oak_geo <- subset(geo_info, plant.sp == oak)
#
#     ## ii) For all [other] plant host spp. of the current Agrilus sp.
#     hosts_geo <- subset(geo_info, plant.sp %in% hosts)
#
#     ## Now, for all individuals from this oak sp. with coord. info from GBIF, find the
#     ## closest host in hosts_geo
#     min.dists <- c()
#
#     for (j in 1:nrow(oak_geo)){
#       ## Compute min dist (in km) to find the closest host in hosts_geo
#       ## Suppressed message when using "cheap" measure: "Maximum distance is > 100km.
#       ## The 'cheap' measure is inaccurate over such large distances, you'd likely be
#       ## better using a different 'measure'."
#       suppressMessages(
#         dists <- geodist(x = oak_geo[j,], y = hosts_geo, measure = "cheap")/1e3
#       )
#       min.dist <- min(dists)
#       min.dists <- c(min.dists, min.dist)
#     }
#
#     ## For this oak sp., find the mean min.dist and median min.dist
#     min.dist.mean <- mean(min.dists)
#     min.dist.median <- median(min.dists)
#   }
#
#   ## Append the info for the current oak.sp to interaction.data
#   interaction_data[i,]$min.dist.mean <- min.dist.mean
#   interaction_data[i,]$min.dist.median <- min.dist.median
# }
# options(warn=1)
#
# ## Now transform distances
# ## First, multiply values * normal distr so that things that are close have a stronger
# ## weight. Use dnorm(), with mean = 0 and SD = 50 so that pairs that are > 100 km
# ## km away have very low values. The distr is rescaled so that scores look prettier.
# ## NB, dnorm(0, 0, 50) [max val] and dnorm(1e5, 0, 50) [min val] are just to make
# ## sure that all pairs are rescaled equally
# ## Mean
# l <- length(interaction_data$min.dist.mean)
# min.dist.mean.norm <- scales::rescale(-c(dnorm(0, 0, 50),
#                                          dnorm(1e5, 0, 50),
#                                          dnorm(interaction_data$min.dist.mean,
#                                                0, 50)))[3:(l+2)]
# interaction_data$min.dist.mean.norm <- min.dist.mean.norm; rm(min.dist.mean.norm, l)
#
# ## Median
# l <- length(interaction_data$min.dist.median)
# min.dist.median.norm <- scales::rescale(-c(dnorm(0, 0, 50),
#                                            dnorm(1e5, 0, 50),
#                                            dnorm(interaction_data$min.dist.median,
#                                                  0, 50)))[3:(l+2)]
# interaction_data$min.dist.median.norm <- min.dist.median.norm
# rm(min.dist.median.norm, l)
#
# ## Now, log distances to check if a log distribution fits the data better
# ## Mean
# interaction_data$min.dist.mean.log <- scales::rescale(log(
#   interaction_data$min.dist.mean + 0.1))
# ## Median
# interaction_data$min.dist.median.log <- scales::rescale(log(
#   interaction_data$min.dist.median + 0.1))
#
# ## Finally, change "Quercus" for "Q.":
# interaction_data$quercus.sp <-  gsub("Quercus", "Q.",
#                                      interaction_data$quercus.sp,
#                                      fixed = TRUE)
#
#
# ## Write table by appending (only add header if file doesn't already exist)
# interaction_data <- interaction_data[first:(first+31), ]
# nrow(interaction_data)
# write.table(x = interaction_data,
#             file = "tmp/interaction_data_oak_hosts.txt",
#             row.names = FALSE,
#             col.names = !file.exists("tmp/interaction_data_oak_hosts.txt"),
#             sep = "\t", quote = FALSE, append = TRUE)
# print("done")
#
# # ## NB: Note that the distance calculated is, inf fact, very different for large dists
# # ## when using the 'cheap' metric.
# # ## E.g., in our DF, the biggest imputed mean min. dist is ca. 23k km for Q. litoralis
# # ## and A. chiricahuae (only host = Q. hypoleucoides). However, if we use "harvestine"
# # ## or "geodesic" as our measure, we get ca. 14k km. Still, the transformed distances
# # ## performed better in our models (see analyses below), so this is not likely to be
# # ## something we'd need to worry about...
# # ggplot() +
# #   borders("world", colour = "gray50", fill = "gray50") +
# #   # coord_cartesian(xlim = c(-10, 50), ylim = c(37, 70)) +
# #   geom_point(aes(lon, lat),
# #              data = geo_info[geo_info$plant.sp == "Quercus litoralis",]) +
# #   geom_point(aes(lon, lat),
# #              data = geo_info[geo_info$plant.sp == "Quercus hypoleucoides",],
# #              col = "red")
# #
# # dists <- geodist(x = geo_info[geo_info$plant.sp == "Quercus litoralis",],
# #                  y = geo_info[geo_info$plant.sp == "Quercus hypoleucoides",],
# #                  measure = "haversine")/1e3
# # m <- c()
# # for (i in 1:nrow(dists)) {
# #   m <- c(m, min(dists[i, ]))
# # }
# # mean(m); rm(dists, m)
#
# # ## NB2: Also note that, when using measure = "cheap", the distance values depend on
# # ## the number of locations that the distances are being computed for, e.g.: [note that
# # ## this doesn't happen if measure is set to "haversine" or "geodesic"]
# # geodist(x = geo_info[geo_info$plant.sp == "Quercus uxoris",][1,],
# #         y = geo_info[geo_info$plant.sp  == "Quercus imbricaria", ][1,],
# #         measure = "cheap")/1e3
# #
# # geodist(x = geo_info[geo_info$plant.sp == "Quercus uxoris",][1,],
# #         y = geo_info[geo_info$plant.sp  == "Quercus imbricaria", ][1:5,],
# #         measure = "cheap")/1e3
#
#
# ## B) EXPLORE THE DATA
# interaction_data <- read.table("tmp/interaction_data_oak_hosts.txt",
#                                header = T, sep = "\t")
#
# ## i) Have a look at min.dist.mean and min.dist.mean
# ## min.dist.mean
# min(interaction_data$min.dist.mean)
# subset(interaction_data, min.dist.mean == min(interaction_data$min.dist.mean))
# hist(interaction_data$min.dist.mean[interaction_data$interaction == 1],
#      seq(0, 3e04, length.out = 30), main = "mean(min.dist) hosts - another host sp.
#      of the current Agrilus sp.", xlab = "mean(dist.min)")
# hist(interaction_data$min.dist.mean[interaction_data$interaction == 0],
#      seq(0, 3e04, length.out = 30), main = "mean(min.dist) non-hosts - host sp.
#      of the current Agrilus sp.", xlab = "mean(dist.min)")
#
# ## min.dist.median
# min(interaction_data$min.dist.median)
# subset(interaction_data, min.dist.median == min(interaction_data$min.dist.median))
# hist(interaction_data$min.dist.median[interaction_data$interaction == 1],
#      seq(0, 3e04, length.out = 30), main = "median(min.dist) hosts - another host sp.
#      of the current Agrilus sp.", xlab = "median(dist.min)")
# hist(interaction_data$min.dist.median[interaction_data$interaction == 0],
#      seq(0, 3e04, length.out = 30), main = "median(min.dist) non-hosts - host sp.
#      of the current Agrilus sp.", xlab = "median(dist.min)")
#
#
# ## ii) Explore max. difference between 'mean' and 'median' distances
# ## Plot histogram for min.dist.mean - min.dist.median
# hist(interaction_data$min.dist.mean - interaction_data$min.dist.median,
#      breaks = 100)
#
# ## Plot histogram for $min.dist.mean.norm - min.dist.median.norm
# hist(interaction_data$min.dist.mean.norm - interaction_data$min.dist.median.norm,
#      breaks = 100)
#
# ## Check max. difference between mean and median min. distances [could also do min.]
# ## Turns out it's for Q. acutissima [mostly present in NE America and SE Asia] - A.
# ## albocomus (hosts = Q. emoryi & Q. grisea [both NW America]) [in this case, preds
# ## using median perform better than w/ mean]
# interaction_data[which.max(interaction_data$min.dist.mean -
#                              interaction_data$min.dist.median), ]
#
# ## Plot distribution of hosts
# ggplot(geo_info[geo_info$plant.sp == "Quercus emoryi" |
#                   geo_info$plant.sp == "Quercus grisea",],
#        aes(x = lon, y = lat, colour = plant.sp)) +
#   coord_fixed() + borders("world", colour = "lightgrey", fill = "lightgrey") +
#   geom_point() +
#   scale_color_manual(values = c("Quercus emoryi" = scales::alpha('darkblue', 1),
#                                 'Quercus grisea' = scales::alpha('coral', 0.3))) +
#   # geom_point(data = poe_geo, aes(x = lon, y = lat),
#   # colour = scales::alpha('darkblue', 1)) +
#   labs(colour = "Quercus spp.") +
#   theme_bw() +
#   theme(legend.text = element_text(size = 15),
#         legend.title = element_text(size = 20))
#
# ## Have a look at the min. distances
# dists <- geodist(x = subset(geo_info, plant.sp == "Quercus acutissima"),
#                  y = subset(geo_info, plant.sp == "Quercus emoryi" |
#                               plant.sp == "Quercus grisea"),
#                  measure = "haversine")/1e3
# min_d <- c()
# for (i in 1:nrow(dists)) {
#   min_d <- c(min_d, min(dists[i,]))
# }
#
# ## Plot a histogram to see what's going on
# hist(min_d, breaks = 50) ## makes sense
# abline(v = interaction_data[which.max(interaction_data$min.dist.mean -
#                                         interaction_data$min.dist.median),
#                             ]$min.dist.mean,
#        col = "red", lwd = 2, lty = "dashed")
# abline(v = interaction_data[which.max(interaction_data$min.dist.mean -
#                                         interaction_data$min.dist.median),
#                             ]$min.dist.median,
#        col = "blue", lwd = 2, lty = "dashed")
#
#
# ## iii) Explore max. difference between norm 'mean' and norm 'median' distances
# ## Check max. diff. between norm mean and median min. dists [could also do min.]
# ## Turns out it's for Q. acutissima [mostly present in NE America and SE Asia] - A.
# ## osburni (hosts = Q. virginiana [NE America]  & Q. phellos [W Europe & NE America])
# interaction_data[which.max(interaction_data$min.dist.mean.norm -
#                              interaction_data$min.dist.median.norm), ]
#
# ## Plot geo location of Q. acutissima (black) and hosts (blue)
# plot(geo_info[geo_info$plant.sp == "Quercus acutissima", 2:3],
#      xlim = c(-180, 180), ylim = c(-90, 90), pch = 16)
# points(geo_info[geo_info$plant.sp == "Quercus virginiana", 2:3],
#        col = scales::alpha("steelblue", 0.8), pch = 8)
# points(geo_info[geo_info$plant.sp == "Quercus phellos", 2:3],
#        col = scales::alpha("steelblue", 0.8), pch = 8)
#
# ## Have a look at the min. distances
# dists <- geodist(x = subset(geo_info, plant.sp == "Quercus acutissima"),
#                  y = subset(geo_info, plant.sp == "Quercus virginiana" |
#                               plant.sp == "Quercus phellos"),
#                  measure = "haversine")/1e3
# min_d <- c()
# for (i in 1:nrow(dists)) {
#   min_d <- c(min_d, min(dists[i,]))
# }
#
# ## Plot a histogram to see what's going on
# hist(min_d, breaks = 50) ## makes sense
# abline(v = interaction_data[which.max(interaction_data$min.dist.mean.norm -
#                                         interaction_data$min.dist.median.norm),
#                             ]$min.dist.mean,
#        col = "red", lwd = 2, lty = "dashed")
# abline(v = interaction_data[which.max(interaction_data$min.dist.mean.norm -
#                                         interaction_data$min.dist.median.norm),
#                             ]$min.dist.median,
#        col = "blue", lwd = 2, lty = "dashed")
#
#
# ## iv) Have a look at A. ussuricola (prediction is better using the mean instead of
# ## the median)
# model_predictions[model_predictions$agrilus.sp == "A. ussuricola" &
#                     model_predictions$host.status == 1,]
#
# min.dists <- c()
# for (j in 1:nrow(geo_info[geo_info$plant.sp == "Quercus acutissima",])) {
#   suppressMessages(
#     dists <- geodist(x = geo_info[geo_info$plant.sp == "Quercus acutissima",][j,],
#                      y = geo_info[geo_info$plant.sp == "Quercus serrata",],
#                      measure = "cheap")/1e3
#   )
#   min.dist <- min(dists)
#   min.dists <- c(min.dists, min.dist)
# }
# hist(min.dists, breaks = 50)     ## bimodal
# abline(v = mean(min.dists), col = "red")
# abline(v = median(min.dists), col = "blue")
#
#
#
## C) PLOT THE 3 GEO DIST METRCIS TO COMPARE THEM
geo_dists <- data.frame(matrix(c(
                        sort(interaction_data$min.dist.mean, decreasing = F),
                        sort(interaction_data$min.dist.mean.norm, decreasing = T),
                        sort(interaction_data$min.dist.mean.log, decreasing = F)),
                        ncol = 3))
colnames(geo_dists) <- c("raw.dist", "norm.dist", "log.dist")


## i) min.dist.mean
# svg("results/figures/distance_comparisons.svg")
raw_plot <- ggplot(data = geo_dists,
                   aes(y = raw.dist, x = seq(1, length(raw.dist)))) +
  geom_point() +
  geom_point(color='steelblue') +
  ggtitle("Distance (km)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("index") + ylab("distance (km)") +
  theme(panel.border = element_blank(),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())


## ii) min.dist.mean.norm
norm_plot <- ggplot(data = geo_dists,
                    aes(y = -norm.dist, x = seq(1, length(norm.dist)))) +
  geom_point() +
  geom_point(color='cornsilk3') +
  ggtitle("Distance (normal-weighted)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("index") + ylab("distance (normal-weighted)") +
  theme(panel.border = element_blank(),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())



## iii) min.dist.mean.log
log_plot <- ggplot(data = geo_dists,
                   aes(y = log.dist, x = seq(1, length(log.dist)))) +
  geom_point() +
  geom_point(color='coral') +
  ggtitle("Distance (log-transformed)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("index") + ylab("distance (log-transformed)") +
  theme(panel.border = element_blank(),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

gridExtra::grid.arrange(raw_plot, norm_plot, log_plot, nrow = 2)
# dev.off()



#
##### VI)   RUN MODELS: BRMS LOG-ODDS                                                #####
## A) PREPARE DATA
## i) Read in interaction data
interaction_data <- read.table(paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/",
                                      "01_original_models/",
                                      "tmp/interaction_data_oak_hosts.txt"),
                               header = T, sep = "\t")
nrow(interaction_data) == 236*32  ## 236 oak spp.*32 Agrilus spp.
interaction_data$quercus.sp <- gsub("-", "", interaction_data$quercus.sp, fixed = T)


## ii) Add info on total no. of total GBIF entries for each oak species
## First, read (cleaned) GBIF geo info DF and only retain info for oaks
geo_info <- read.table(paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/",
                              "01_original_models/",
                              "tmp/gbif_plants_clean.tsv"),
                       header = T, sep = "\t")
geo_info <- geo_info[, c("plant.sp", "lon", "lat")]
geo_info <- geo_info[grepl("Quercus", geo_info$plant.sp),]

## Extract no. occurrences
oak_occurrences <- data.frame(table(geo_info$plant.sp))
colnames(oak_occurrences) <- c("quercus.sp", "gbif.entries")
oak_occurrences$quercus.sp <- as.character(oak_occurrences$quercus.sp)
oak_occurrences$quercus.sp <- gsub("-", "", oak_occurrences$quercus.sp, fixed = T)
oak_occurrences$quercus.sp <-  gsub("Quercus", "Q.", oak_occurrences$quercus.sp,
                                    fixed = TRUE)

## Append this info to interaction_data
interaction_data$ID <- 1:nrow(interaction_data)
interaction_data <- merge(interaction_data, oak_occurrences, by = "quercus.sp",
                          all = T)
# interaction_data[is.na(interaction_data$gbif.entries),]$gbif.entries <- 0
interaction_data <- interaction_data[order(interaction_data$ID), ]
interaction_data$ID <- NULL

rm(geo_info, oak_occurrences)

## iii) Retrieve info on Agrilus hosts
agrilus_hosts <- read.table(paste0("/data/SBCS-NicholsLab/elvirahg/",
                                   "oak_analyses/01_original_models/",
                                   "input/agrilus_quercus_hosts.txt"),
                            header = F, sep = "\t")
colnames(agrilus_hosts) <- c("agrilus", "hosts")
agrilus_hosts <- rbind(agrilus_hosts,
                       read.table(paste0("/data/SBCS-NicholsLab/elvirahg/",
                                         "oak_analyses/01_original_models/",
                                         "tmp/non_oak_hosts.tsv"),
                                  header = T, sep = "\t"))
agrilus_hosts <- agrilus_hosts[,c(2,1)]
colnames(agrilus_hosts) <- c("plant.sp", "agrilus.sp")
agrilus_hosts$agrilus.sp <- gsub("Agrilus", "A.", agrilus_hosts$agrilus.sp)

## Remove any hosts that are not in the phylogeny (only Q. gambelii) from agrilus_hosts
agrilus_hosts <- agrilus_hosts[!grepl("gambelii", agrilus_hosts$plant.sp),]

## IMPORTANT (!): I have also decided to remove any reported hosts that are non-native
## for those Agrilus species that have been introduced into new areas.
## A. auroguttatus: introduced from AZ (USA) to CA (USA). In CA, reported novel larval
## hosts are Q. kelloggii (not in my DS), Q. chrysolepis (not in my DS) & Q. agrifolia
## (in my DS).
## A. bilineatus: introduced from N America to Turkey (no new reported hosts so far).
## A. sulcicollis: introduced from Europe to N America, where Q. macrocarpa (in my DS)
## is a new reported larval host.
# head(agrilus_hosts[agrilus_hosts$plant.sp == "Quercus macrocarpa" &
#                    agrilus_hosts$agrilus.sp == "A. sulcicollis",])

agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Quercus agrifolia"
                                 & agrilus_hosts$agrilus.sp == "A. auroguttatus")), ]
agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Quercus macrocarpa"
                                 & agrilus_hosts$agrilus.sp == "A. sulcicollis")), ]

## IMPORTANT (!): Let's also remove an larval hosts according to J&P(2014) (i.e.,
## there's at least one reference with '!') but with a  confidence index < 3.
## A. graminis: Corylus avellana, Ostrya carpinifolia, Euonymus europaeus and Castanea
## sativa.
## A. obscuricollis: Ficus carica
## A. relegatus (reported as A. dualis)	Quercus robur
agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Ficus carica"
                                 & agrilus_hosts$agrilus.sp == "A. obscuricollis")), ]

agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Quercus robur"
                                 & agrilus_hosts$agrilus.sp == "A. relegatus")), ]

agrilus_hosts <- agrilus_hosts[-(which((agrilus_hosts$plant.sp == "Corylus avellana" |
                                        agrilus_hosts$plant.sp == "Ostrya carpinifolia" |
                                        agrilus_hosts$plant.sp == "Euonymus europaeus" |
                                        agrilus_hosts$plant.sp == "Castanea sativa") &
                                        agrilus_hosts$agrilus.sp == "A. graminis")), ]


## iv) Create phylogenetic cov matrix
## Read nexus tree created with phylocom nodesig (first tree: host status info)
oak_phylo <- paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/",
                    "01_original_models/input/quercus_nodesig_result.nex")
oak_phylo <- read.nexus(oak_phylo)[[1]]
# plot(oak_phylo)

## Edit node label information
oak_phylo$node.label <- gsub("'", "", oak_phylo$node.label)
oak_phylo$node.label <- gsub("N", "", oak_phylo$node.label)

## Replace "_" for " "
oak_phylo$tip.label <- gsub(pattern = "_", replacement = " ", oak_phylo$tip.label)
oak_phylo$tip.label <- gsub(pattern = "Quercus", replacement = "Q.",
                            oak_phylo$tip.label, fixed = TRUE)

## Create a phylo object with info on the relationship between species. This can
## be used to construct a cov matrix of species (Hadfield & Nakagawa, 2010)
oak_phylo_cov <- vcv.phylo(oak_phylo)
rm(oak_phylo)


## v) Add phylo distance info to interaction_data
## Create a new column for phylo distance info
interaction_data$phylo.dist.mean <- NA
interaction_data$phylo.dist.min <- NA

## For every possible Quercus - Agrilus interaction, find the **oak** host species of
## the current Agrilus species, and then find the mean phylogenetic distance of the
## focal Quercus species to the hosts, and normalise it. If the focal Quercus species
## happens to be the only host, then set the distance to the max distance.
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

## Now normalise distances
# plot(phylo.dist.mean ~ phylo.dist.min, data = interaction_data)
interaction_data$phylo.dist.mean <- scales::rescale(-c(interaction_data$phylo.dist.mean))
interaction_data$phylo.dist.min <- scales::rescale(-c(interaction_data$phylo.dist.min))

## And compare
# plot(phylo.dist.mean ~ phylo.dist.min, data = interaction_data)


## B) RUN MODELS - PRELIMINARY MODELS
# ## i) Geographic distance models: which geographic distance metric performs best?
# ## Prepare model formulas (all possible variable combinations)
# formulas <- list(
#   oak_modA1 = "interaction ~ min.dist.mean",
#   oak_modB1 = "interaction ~ min.dist.median",
#   oak_modC1 = "interaction ~ min.dist.mean.norm",
#   oak_modD1 = "interaction ~ min.dist.median.norm",
#   oak_modE1 = "interaction ~ min.dist.mean.log",
#   oak_modF1 = "interaction ~ min.dist.median.log")
#
# ## Run models
# for (i in 1:length(formulas)){
#   oak_mod <- brm(formulas[[i]],
#                  data = interaction_data,
#                  family = "bernoulli",
#                  iter = 7000, cores = 4,
#                  control = list(adapt_delta = 0.995,  max_treedepth = 15),
#                  data2 = list(oak_phylo_cov = oak_phylo_cov),
#                  save_pars = save_pars(all = TRUE)
#   )
#   assign(names(formulas[i]), oak_mod)
# }
#
# ## Compare models and choose best-performing geographic distance metric
# loo_geo_dist <- loo(oak_modA1, oak_modB1, oak_modC1,
#                     oak_modD1, oak_modE1, oak_modF1,
#                     save_psis = T, moment_match = T,
#                     reloo = T, compare = T)
# loo_geo_dist$diff
#
# ## Best-fitting models = min.dist.mean.norm & min.dist.median.norm
# #             elpd_diff   se_diff
# # oak_modC1    0.0        0.0
# # oak_modD1   -4.7        6.0
# # oak_modF1   -57.9       9.4
# # oak_modE1   -65.3       9.7
# # oak_modA1   -146.0      17.1
# # oak_modB1   -147.7      17.1
#
# ## Transform object into DF and add model diff information
# loo_geo_dist_df <- data.frame(matrix(nrow = length(loo_geo_dist$diffs[,1]),
#                                      ncol = 3))
# colnames(loo_geo_dist_df) <- c("elpd_diff", "se_diff", "sig")
# rownames(loo_geo_dist_df) <- names(loo_geo_dist$diffs[,1])
# loo_geo_dist_df$elpd_diff <- loo_geo_dist$diffs[,1]
# loo_geo_dist_df$se_diff <- loo_geo_dist$diffs[,2]
# loo_geo_dist_df$sig <- FALSE
#
# ## Now add model-specific information
# loo_geo_dist_df <-  loo_geo_dist_df[order(rownames(loo_geo_dist_df)),]
# loo_geo_dist_df$elpd_loo.estimate <- NA
# loo_geo_dist_df$elpd_loo.se <- NA
# loo_geo_dist_df$p_loo.estimate <- NA
# loo_geo_dist_df$p_loo.se <- NA
# loo_geo_dist_df$looic.estimate <- NA
# loo_geo_dist_df$looic.se <- NA
#
# for (i in 1:nrow(loo_geo_dist_df)) {
#   loo_geo_dist_df$elpd_loo.estimate[i] <- loo_geo_dist$loos[[i]]$estimates[1,1]
#   loo_geo_dist_df$elpd_loo.se[i] <- loo_geo_dist$loos[[i]]$estimates[1,1]
#   loo_geo_dist_df$p_loo.estimate[i] <- loo_geo_dist$loos[[i]]$estimates[2,1]
#   loo_geo_dist_df$p_loo.se[i] <- loo_geo_dist$loos[[i]]$estimates[2,2]
#   loo_geo_dist_df$looic.estimate[i] <- loo_geo_dist$loos[[i]]$estimates[3,1]
#   loo_geo_dist_df$looic.se[i] <- loo_geo_dist$loos[[i]]$estimates[3,2]
# }
# loo_geo_dist_df <- loo_geo_dist_df[order(loo_geo_dist_df$elpd_diff,
#                                          decreasing = TRUE),]
# head(loo_geo_dist_df)
#
# # write.table(loo_geo_dist_df, "results/looic_model_values_geo.tsv",
# #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
#
#
# ## ii) Fixed, Agrilus-species specific oak phylogenetic distance models: which metric
# ## performs best?
# ## This is just exploratory, out of curiosity
# ## Prepare model formulas (all possible variable combinations)
# formulas <- list(
#   oak_modA2 = "interaction ~ phylo.dist.mean",
#   oak_modB2 = "interaction ~ phylo.dist.min",
#   oak_modC2 = "interaction ~ phylo.dist.mean + phylo.dist.min",
#   oak_modD2 = "interaction ~ phylo.dist.mean + (phylo.dist.mean | agrilus.sp)",
#   oak_modE2 = "interaction ~ phylo.dist.min + (phylo.dist.mean | agrilus.sp)",
#   oak_modF2 = paste("interaction ~ phylo.dist.mean + phylo.dist.min",
#                     "+ (phylo.dist.mean + phylo.dist.min | agrilus.sp)"),
#   oak_modG2 = "interaction ~ phylo.dist.mean + (1 | agrilus.sp)",
#   oak_modH2 = "interaction ~ phylo.dist.min + (1 | agrilus.sp)",
#   oak_modI2 = paste("interaction ~ phylo.dist.mean + phylo.dist.min",
#                     "+ (1 | agrilus.sp)")
# )
#
# ## Run models
# for (i in 7:length(formulas)){
#   oak_mod <- brm(formulas[[i]],
#                  data = interaction_data,
#                  family = "bernoulli",
#                  iter = 7000, cores = 4,
#                  control = list(adapt_delta = 0.995,  max_treedepth = 15),
#                  data2 = list(oak_phylo_cov = oak_phylo_cov),
#                  save_pars = save_pars(all = TRUE)
#   )
#   assign(names(formulas[i]), oak_mod)
# }
#
# ## Compare models and choose best-performing geographic distance metric
# loo_phy_dist <- loo(oak_modA2, oak_modB2, oak_modC2,
#                     oak_modD2, oak_modE2, oak_modF2,
#                     oak_modG2, oak_modH2, oak_modI2,
#                     save_psis = T, moment_match = T,
#                     reloo = T, compare = T)
# loo_phy_dist$diff
#
# ## Best-fitting models = mean + min + (mean + min | agrilus) // min + (min | agrilus)
# #           elpd_diff     se_diff
# # oak_modF2    0.0        0.0
# # oak_modE2   -9.7        4.5
# # oak_modH2   -23.2       6.0
# # oak_modI2   -23.5       6.0
# # oak_modC2   -24.1       6.6
# # oak_modB2   -24.6       6.8
# # oak_modD2   -24.9       8.6
# # oak_modG2   -45.3       10.0
# # oak_modA2   -64.6       13.6
#
# ## Transform object into DF and add model diff information
# loo_phy_dist_df <- data.frame(matrix(nrow = length(loo_phy_dist$diffs[,1]),
#                                      ncol = 3))
# colnames(loo_phy_dist_df) <- c("elpd_diff", "se_diff", "sig")
# rownames(loo_phy_dist_df) <- names(loo_phy_dist$diffs[,1])
# loo_phy_dist_df$elpd_diff <- loo_phy_dist$diffs[,1]
# loo_phy_dist_df$se_diff <- loo_phy_dist$diffs[,2]
# loo_phy_dist_df$sig <- FALSE
#
# ## Now add model-specific information
# loo_phy_dist_df <-  loo_phy_dist_df[order(rownames(loo_phy_dist_df)),]
# loo_phy_dist_df$elpd_loo.estimate <- NA
# loo_phy_dist_df$elpd_loo.se <- NA
# loo_phy_dist_df$p_loo.estimate <- NA
# loo_phy_dist_df$p_loo.se <- NA
# loo_phy_dist_df$looic.estimate <- NA
# loo_phy_dist_df$looic.se <- NA
#
# for (i in 1:nrow(loo_phy_dist_df)) {
#   loo_phy_dist_df$elpd_loo.estimate[i] <- loo_phy_dist$loos[[i]]$estimates[1,1]
#   loo_phy_dist_df$elpd_loo.se[i] <- loo_phy_dist$loos[[i]]$estimates[1,1]
#   loo_phy_dist_df$p_loo.estimate[i] <- loo_phy_dist$loos[[i]]$estimates[2,1]
#   loo_phy_dist_df$p_loo.se[i] <- loo_phy_dist$loos[[i]]$estimates[2,2]
#   loo_phy_dist_df$looic.estimate[i] <- loo_phy_dist$loos[[i]]$estimates[3,1]
#   loo_phy_dist_df$looic.se[i] <- loo_phy_dist$loos[[i]]$estimates[3,2]
# }
# loo_phy_dist_df <- loo_phy_dist_df[order(loo_phy_dist_df$elpd_diff,
#                                          decreasing = TRUE),]
# head(loo_phy_dist_df)

# write.table(loo_phy_dist_df, "results/looic_model_values_phy.tsv",
#             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
# rm(a, c, d.mean, d.min, h, i, q, r, oak_mod, formulas)
# # save.image(paste("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/",
# #                  "results/oak_models_preliminary.RData", sep = ""))
#
# load(paste("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/",
#            "results/oak_models_preliminary.RData", sep = ""))


## C) RUN MODELS
## i) Make formulas for all models
## Vector with all variables that we're going to use (minus ""variations"")
variables <- c("min.dist.mean.norm", "phylo.dist.min", "gbif.entries",
               "(1 | gr(quercus.sp, cov = oak_phylo_cov))", "(1 | agrilus.sp)")

## All variable combinations
## NB, got this bit from https://stackoverflow.com/questions/40049313/generate-all-
## combinations-of-all-lengths-in-r-from-a-vector
formulas <- do.call("c",
                    lapply(seq_along(variables),
                           function(i) combn(variables, i, FUN = list)))
formulas <- lapply(formulas, paste, collapse = " + ")

## We established above that "min.dist.mean.norm" and "min.dist.median.norm" are the
## best fitting geographic variables - so let's only use there two for the subsequent
## models. Here, every time there is "min.dist.mean.norm" in our model formula, make
## another formula that's the exact same, but with "min.dist.median.norm" instead
## Note that with all 6 geographic distance metrics we'd have 308 models to run, as:
## sum(grepl("min.dist.mean.norm", formulas), na.rm = TRUE)*4 + length(formulas)
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

## We want to compare the performance of the different Quercus sp. - Agrilus sp.
## "specific" phylogenetic distances, and check whether we should add a random slope
## effect to our model. Here, every time there is "phylo.dist.min" in our model formula,
## make another formula that's the exact same, but with "phylo.dist.mean" instead, and
## another one that's "phylo.dist.min + phylo.dist.mean"
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

## Here, every time there is "phylo.dist.min|phylo.dist.mean" plus an Agrilus random
## effect in our model formula, make another formula that's the exact same, but with
## a random slope instead
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

## Add response variable to our formulas
formulas <- gsub("^", "interaction ~ ", formulas)

## Add intercept null model as the first one
formulas <- c(list("interaction ~ 1"), formulas)

## Rename our formulas with useful model names
names(formulas) <- paste("oak_mod", 1:length(formulas))
names(formulas) <- gsub(" ([0-9])$", "00\\1", names(formulas))
names(formulas) <- gsub(" ([0-9][0-9])$", "0\\1", names(formulas))
names(formulas) <- gsub(" ", "", names(formulas))


# ## ii) Run models
# ## Now run models inside a jobscript: ask for 133 iterations, 2G/core, 4 cores, 240 h
# ## (the longest-running jobs take ca. 3 days to run)
# i <- as.numeric(commandArgs(trailingOnly=TRUE))
# oak_mod <- brm(formulas[[i]],
#                data = interaction_data,
#                family = "bernoulli",
#                iter = 7000, cores = 4,
#                control = list(adapt_delta = 0.995,  max_treedepth = 15),
#                data2 = list(oak_phylo_cov = oak_phylo_cov),
#                save_pars = save_pars(all = TRUE)
# )
#
# assign(names(formulas[i]), oak_mod)
#
# # save.image(paste("/data/scratch/btx840/results/", names(formulas[i]), ".RData",
# #                  sep = ""))
#
#
# ## iii) Save models in one .RData file
# files <- paste("/data/scratch/btx840/results/",
#                list.files(path = "/data/scratch/btx840/results/"),
#                sep = "")
# lapply(files, load, .GlobalEnv)
# rm(formulas, oak_mod, a, c, d.mean, d.min, files, h, i, q, r, v, variables)
#
# load(paste("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/",
#            "results/oak_models_preliminary.RData", sep = ""))
#
# # save.image(paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/",
# #                   "01_original_models/results/oak_models_all.RData"))


## iv) Save all summary tables
# load(paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/results/",
#             "oak_models_all.RData"))
# models <- mget(ls(pattern = '^oak_mod[0-9]')); rm(list=setdiff(ls(), c("models")))
# sink(paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/results/",
#             "/summary_tables.txt"), append = TRUE)
# for (i in 1:length(models)) {
#   print(paste("MODEL", i))
#   print(summary(models[[i]]))
# }
# sink()
# rm(models, i)
#
## Make the file look prettier with bash:
# cat("sed -i -e 's/\\[1\\] //g' -e 's/\"//g' -e 's/^MODEL/\\n\\nMODEL/g' ",
#     "-e 's/^   *Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS/",
#     "\\tEstimate\\tEst.Error\\tl-95% CI\\tu-95% CI\\tRhat\\tBulk_ESS\\tTail_ESS/g' ",
#     "'-e s/\\([0-9][0-9\\.][0-9\\.]*\\)  */\\1\\t/g' ",
#     "-e 's/Intercept  */Intercept\\t/g' -e 's/sd(Intercept) */sd(Intercept)\\t/g' ",
#     "-e 's/^\\([a-z\\.][a-z\\.]*\\)  *\\([0-9\\-][0-9\\.]\\)/\\1\\t\\2/g' ",
#     "-e '/^Draws/d' -e '/^and/d' -e '/^scale/d' -e 's/^  *//g' ",
#     "results/summary_tables.txt", sep = "")
# cat("sed -i 'N;s/\\ntotal/ total/;P;D' results/summary_tables.txt")



#
##### VII)  CHECK LOOIC                                                              ####
## NB, see https://discourse.mc-stan.org/t/understanding-looic/13409
## And https://mc-stan.org/loo/reference/loo-glossary.html#se-diff
## Lower LOOIC values indicate a better fit
# load(paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/results/",
#             "oak_models_all.RData"))

## NB, if I do WAIC (waic())), I get a message suggesting to use loo instead:
## "p_waic estimates greater than 0.4. We recommend trying loo instead."

# ## A) RUN LOO ANALYSES
# ## i) Make three sets of comparisons
# # cat(ls(pattern = "oak_mod"), sep = ", ")
# i <- as.numeric(commandArgs(trailingOnly=TRUE))
#
# if (i == 1) {
#   ## Takes about 1 day to run, uses ca. 80 GB
#   rm(list=setdiff(ls(),
#                   c("oak_mod001", "oak_mod002", "oak_mod003", "oak_mod004", "oak_mod005",
#                     "oak_mod006", "oak_mod007", "oak_mod008", "oak_mod009", "oak_mod010",
#                     "oak_mod011", "oak_mod012", "oak_mod013", "oak_mod014", "oak_mod015",
#                     "oak_mod016", "oak_mod017", "oak_mod018", "oak_mod019", "oak_mod020",
#                     "oak_mod021", "oak_mod022", "oak_mod023", "oak_mod024", "oak_mod025",
#                     "oak_mod026", "oak_mod027", "oak_mod028", "oak_mod029", "oak_mod030",
#                     "oak_mod031", "oak_mod032", "oak_mod033", "oak_mod034", "oak_mod035",
#                     "oak_mod036", "oak_mod037", "oak_mod038", "oak_mod039", "oak_mod040",
#                     "oak_mod041", "oak_mod042", "oak_mod043", "oak_mod044", "i")))
#
#   print("start")
#   loo_oak_mods1 <- loo(oak_mod001, oak_mod002, oak_mod003, oak_mod004, oak_mod005,
#                        oak_mod006, oak_mod007, oak_mod008, oak_mod009, oak_mod010,
#                        oak_mod011, oak_mod012, oak_mod013, oak_mod014, oak_mod015,
#                        oak_mod016, oak_mod017, oak_mod018, oak_mod019, oak_mod020,
#                        oak_mod021, oak_mod022, oak_mod023, oak_mod024, oak_mod025,
#                        oak_mod026, oak_mod027, oak_mod028, oak_mod029, oak_mod030,
#                        oak_mod031, oak_mod032, oak_mod033, oak_mod034, oak_mod035,
#                        oak_mod036, oak_mod037, oak_mod038, oak_mod039, oak_mod040,
#                        oak_mod041, oak_mod042, oak_mod043, oak_mod044,
#                        save_psis = TRUE, moment_match = TRUE,
#                        reloo = TRUE, compare = TRUE)
# } else if (i == 2) {
#   ## Takes about 2 days to run, uses ca. 80 GB
#   rm(list=setdiff(ls(),
#                   c("oak_mod045", "oak_mod046", "oak_mod047", "oak_mod048", "oak_mod049",
#                     "oak_mod050", "oak_mod051", "oak_mod052", "oak_mod053", "oak_mod054",
#                     "oak_mod055", "oak_mod056", "oak_mod057", "oak_mod058", "oak_mod059",
#                     "oak_mod060", "oak_mod061", "oak_mod062", "oak_mod063", "oak_mod064",
#                     "oak_mod065", "oak_mod066", "oak_mod067", "oak_mod068", "oak_mod069",
#                     "oak_mod070", "oak_mod071", "oak_mod072", "oak_mod073", "oak_mod074",
#                     "oak_mod075", "oak_mod076", "oak_mod077", "oak_mod078", "oak_mod079",
#                     "oak_mod080", "oak_mod081", "oak_mod082", "oak_mod083", "oak_mod084",
#                     "oak_mod085", "oak_mod086", "oak_mod087", "oak_mod088", "i")))
#
#   print("start")
#   loo_oak_mods2 <- loo(oak_mod045, oak_mod046, oak_mod047, oak_mod048, oak_mod049,
#                        oak_mod050, oak_mod051, oak_mod052, oak_mod053, oak_mod054,
#                        oak_mod055, oak_mod056, oak_mod057, oak_mod058, oak_mod059,
#                        oak_mod060, oak_mod061, oak_mod062, oak_mod063, oak_mod064,
#                        oak_mod065, oak_mod066, oak_mod067, oak_mod068, oak_mod069,
#                        oak_mod070, oak_mod071, oak_mod072, oak_mod073, oak_mod074,
#                        oak_mod075, oak_mod076, oak_mod077, oak_mod078, oak_mod079,
#                        oak_mod080, oak_mod081, oak_mod082, oak_mod083, oak_mod084,
#                        oak_mod085, oak_mod086, oak_mod087, oak_mod088,
#                        save_psis = TRUE, moment_match = TRUE,
#                        reloo = TRUE, compare = TRUE)
# } else if (i == 3) {
#   ## Takes about 5 days to run, uses ca. 100 GB
#   rm(list=setdiff(ls(),
#                   c("oak_mod089", "oak_mod090", "oak_mod091", "oak_mod092", "oak_mod093",
#                     "oak_mod094", "oak_mod095", "oak_mod096", "oak_mod097", "oak_mod098",
#                     "oak_mod099", "oak_mod100", "oak_mod101", "oak_mod102", "oak_mod103",
#                     "oak_mod104", "oak_mod105", "oak_mod106", "oak_mod107", "oak_mod108",
#                     "oak_mod109", "oak_mod110", "oak_mod111", "oak_mod112", "oak_mod113",
#                     "oak_mod114", "oak_mod115", "oak_mod116", "oak_mod117", "oak_mod118",
#                     "oak_mod119", "oak_mod120", "oak_mod121", "oak_mod122", "oak_mod123",
#                     "oak_mod124", "oak_mod125", "oak_mod126", "oak_mod127", "oak_mod128",
#                     "oak_mod129", "oak_mod130", "oak_mod131", "oak_mod132", "i")))
#
#   print("start")
#   loo_oak_mods3 <- loo(oak_mod089, oak_mod090, oak_mod091, oak_mod092, oak_mod093,
#                        oak_mod094, oak_mod095, oak_mod096, oak_mod097, oak_mod098,
#                        oak_mod099, oak_mod100, oak_mod101, oak_mod102, oak_mod103,
#                        oak_mod104, oak_mod105, oak_mod106, oak_mod107, oak_mod108,
#                        oak_mod109, oak_mod110, oak_mod111, oak_mod112, oak_mod113,
#                        oak_mod114, oak_mod115, oak_mod116, oak_mod117, oak_mod118,
#                        oak_mod119, oak_mod120, oak_mod121, oak_mod122, oak_mod123,
#                        oak_mod124, oak_mod125, oak_mod126, oak_mod127, oak_mod128,
#                        oak_mod129, oak_mod130, oak_mod131, oak_mod132,
#                        save_psis = TRUE, moment_match = TRUE,
#                        reloo = TRUE, compare = TRUE)
# }
#
# ## Save image
# # save.image(paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/",
# #                   "results/oak_models_looic", i, ".RData"))
# # print("saved")



## B) EXPLORE RESULTS
# ## i) Group 1
# ## Load data
# load("/data/scratch/btx840/results/oak_models_looic1.RData")
#
# ## Explore results
# # loo_oak_mods1
# ##             elpd_diff  se_diff
# ## oak_mod043   0.0       0.0
# ## oak_mod044  -4.4       3.2     *
# ## oak_mod018  -9.9       4.9
# ## oak_mod019  -24.7      8.4
# ## oak_mod037  -35.7      9.3
#
# ## Transform object into DF and add model diff information
# loo_oak_mods1_df <- data.frame(matrix(nrow = length(loo_oak_mods1$diffs[,1]),
#                                      ncol = 2))
# colnames(loo_oak_mods1_df) <- c("elpd_diff", "se_diff")
# rownames(loo_oak_mods1_df) <- names(loo_oak_mods1$diffs[,1])
# loo_oak_mods1_df$elpd_diff <- loo_oak_mods1$diffs[,1]
# loo_oak_mods1_df$se_diff <- loo_oak_mods1$diffs[,2]
# # loo_oak_mods1_df$sig <- FALSE
# #
# # ## Check significance
# # for (i in 1:nrow(loo_oak_mods1_df)) {
# #   sig1 <- loo_oak_mods1_df$elpd_diff[i] + loo_oak_mods1_df$se_diff[i]
# #   sig2 <- loo_oak_mods1_df$elpd_diff[i] - loo_oak_mods1_df$se_diff[i]
# #
# #   if((sign(sig1) == 1 & sign(sig2) == -1) | sig1 == 0 | sig2 == 0) {
# #     loo_oak_mods1_df$sig[i] <- TRUE
# #   } else {
# #     break
# #   }
# # }; rm(i, sig1, sig2)
#
# ## Now add model-specific information
# loo_oak_mods1_df <-  loo_oak_mods1_df[order(rownames(loo_oak_mods1_df)),]
# loo_oak_mods1_df$elpd_loo.estimate <- NA
# loo_oak_mods1_df$elpd_loo.se <- NA
# loo_oak_mods1_df$p_loo.estimate <- NA
# loo_oak_mods1_df$p_loo.se <- NA
# loo_oak_mods1_df$looic.estimate <- NA
# loo_oak_mods1_df$looic.se <- NA
#
# for (i in 1:nrow(loo_oak_mods1_df)) {
#   loo_oak_mods1_df$elpd_loo.estimate[i] <- loo_oak_mods1$loos[[i]]$estimates[1,1]
#   loo_oak_mods1_df$elpd_loo.se[i] <- loo_oak_mods1$loos[[i]]$estimates[1,2]
#   loo_oak_mods1_df$p_loo.estimate[i] <- loo_oak_mods1$loos[[i]]$estimates[2,1]
#   loo_oak_mods1_df$p_loo.se[i] <- loo_oak_mods1$loos[[i]]$estimates[2,2]
#   loo_oak_mods1_df$looic.estimate[i] <- loo_oak_mods1$loos[[i]]$estimates[3,1]
#   loo_oak_mods1_df$looic.se[i] <- loo_oak_mods1$loos[[i]]$estimates[3,2]
# }
# loo_oak_mods1_df <- loo_oak_mods1_df[order(loo_oak_mods1_df$elpd_diff,
#                                          decreasing = TRUE),]
# head(loo_oak_mods1_df)
#
# # write.table(loo_oak_mods1_df, "results/looic_model1_values.tsv",
# #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
# # write.table(loo_oak_mods1$ic_diffs__, "results/looic_model1_comparisons.tsv",
# #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
#
#
# ## ii) Group 2
# ## Load data
# # load("/data/scratch/btx840/results/oak_models_looic2.RData")
#
# ## Explore results
# # loo_oak_mods2
# ##              elpd_diff  se_diff
# ## oak_mod045    0.0       0.0
# ## oak_mod083   -0.3       3.0
# ## oak_mod085   -1.1       3.0
# ## oak_mod084   -4.5       4.3
# ## oak_mod065   -7.9       4.6
# ## oak_mod061   -8.7       5.6     *
# ## oak_mod046   -11.1      5.5
# ## oak_mod086   -11.6      6.6
# ## oak_mod048   -11.9      5.5
#
# ## Transform object into DF and add model diff information
# loo_oak_mods2_df <- data.frame(matrix(nrow = length(loo_oak_mods2$diffs[,1]),
#                                       ncol = 2))
# colnames(loo_oak_mods2_df) <- c("elpd_diff", "se_diff")
# rownames(loo_oak_mods2_df) <- names(loo_oak_mods2$diffs[,1])
# loo_oak_mods2_df$elpd_diff <- loo_oak_mods2$diffs[,1]
# loo_oak_mods2_df$se_diff <- loo_oak_mods2$diffs[,2]
# # loo_oak_mods2_df$sig <- FALSE
# #
# # ## Check significance
# # for (i in 1:nrow(loo_oak_mods2_df)) {
# #   sig1 <- loo_oak_mods2_df$elpd_diff[i] + loo_oak_mods2_df$se_diff[i]
# #   sig2 <- loo_oak_mods2_df$elpd_diff[i] - loo_oak_mods2_df$se_diff[i]
# #
# #   if((sign(sig1) == 1 & sign(sig2) == -1) | sig1 == 0 | sig2 == 0) {
# #     loo_oak_mods2_df$sig[i] <- TRUE
# #   } else {
# #     break
# #   }
# # }; rm(i, sig1, sig2)
#
# ## Now add model-specific information
# loo_oak_mods2_df <-  loo_oak_mods2_df[order(rownames(loo_oak_mods2_df)),]
# loo_oak_mods2_df$elpd_loo.estimate <- NA
# loo_oak_mods2_df$elpd_loo.se <- NA
# loo_oak_mods2_df$p_loo.estimate <- NA
# loo_oak_mods2_df$p_loo.se <- NA
# loo_oak_mods2_df$looic.estimate <- NA
# loo_oak_mods2_df$looic.se <- NA
#
# for (i in 1:nrow(loo_oak_mods2_df)) {
#   loo_oak_mods2_df$elpd_loo.estimate[i] <- loo_oak_mods2$loos[[i]]$estimates[1,1]
#   loo_oak_mods2_df$elpd_loo.se[i] <- loo_oak_mods2$loos[[i]]$estimates[1,2]
#   loo_oak_mods2_df$p_loo.estimate[i] <- loo_oak_mods2$loos[[i]]$estimates[2,1]
#   loo_oak_mods2_df$p_loo.se[i] <- loo_oak_mods2$loos[[i]]$estimates[2,2]
#   loo_oak_mods2_df$looic.estimate[i] <- loo_oak_mods2$loos[[i]]$estimates[3,1]
#   loo_oak_mods2_df$looic.se[i] <- loo_oak_mods2$loos[[i]]$estimates[3,2]
# }
# loo_oak_mods2_df <- loo_oak_mods2_df[order(loo_oak_mods2_df$elpd_diff,
#                                            decreasing = TRUE),]
# head(loo_oak_mods2_df)
#
# # write.table(loo_oak_mods2_df, "results/looic_model2_values.tsv",
# #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
# # write.table(loo_oak_mods2$ic_diffs__, "results/looic_model2_comparisons.tsv",
# #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
#
#
#
# ## iii) Group 3
# ## Load data
# # load("/data/scratch/btx840/results/oak_models_looic3.RData")
#
# ## Explore results
# # loo_oak_mods3
# #             elpd_diff  se_diff
# # oak_mod126   0.0       0.0
# # oak_mod106  -0.2       3.0
# # oak_mod102  -0.3       3.1
# # oak_mod122  -1.2       1.0
# # oak_mod101  -3.0       4.6
# # oak_mod105  -3.9       4.6
# # oak_mod121  -4.5       4.0
# # oak_mod104  -4.5       4.7
# # oak_mod125  -4.9       3.9
# # oak_mod103  -5.1       4.8
# # oak_mod124  -5.4       3.7
# # oak_mod123  -6.1       4.0
# # oak_mod132  -7.0       6.0
# # oak_mod112  -7.2       6.5
# # oak_mod128  -8.6       5.9
# # oak_mod108  -8.9       6.5
# # oak_mod113  -10.6      5.6  *
# # oak_mod127  -15.1      6.7
# # oak_mod107  -15.2      7.0
# # oak_mod131  -15.6      6.7
#
# ## Transform object into DF and add model diff information
# loo_oak_mods3_df <- data.frame(matrix(nrow = length(loo_oak_mods3$diffs[,1]),
#                                       ncol = 2))
# colnames(loo_oak_mods3_df) <- c("elpd_diff", "se_diff")
# rownames(loo_oak_mods3_df) <- names(loo_oak_mods3$diffs[,1])
# loo_oak_mods3_df$elpd_diff <- loo_oak_mods3$diffs[,1]
# loo_oak_mods3_df$se_diff <- loo_oak_mods3$diffs[,2]
# # loo_oak_mods3_df$sig <- FALSE
# #
# # ## Check significance
# # for (i in 1:nrow(loo_oak_mods3_df)) {
# #   sig1 <- loo_oak_mods3_df$elpd_diff[i] + loo_oak_mods3_df$se_diff[i]
# #   sig2 <- loo_oak_mods3_df$elpd_diff[i] - loo_oak_mods3_df$se_diff[i]
# #
# #   if((sign(sig1) == 1 & sign(sig2) == -1) | sig1 == 0 | sig2 == 0) {
# #     loo_oak_mods3_df$sig[i] <- TRUE
# #   } else {
# #     break
# #   }
# # }; rm(i, sig1, sig2)
#
# ## Now add model-specific information
# loo_oak_mods3_df <-  loo_oak_mods3_df[order(rownames(loo_oak_mods3_df)),]
# loo_oak_mods3_df$elpd_loo.estimate <- NA
# loo_oak_mods3_df$elpd_loo.se <- NA
# loo_oak_mods3_df$p_loo.estimate <- NA
# loo_oak_mods3_df$p_loo.se <- NA
# loo_oak_mods3_df$looic.estimate <- NA
# loo_oak_mods3_df$looic.se <- NA
#
# for (i in 1:nrow(loo_oak_mods3_df)) {
#   loo_oak_mods3_df$elpd_loo.estimate[i] <- loo_oak_mods3$loos[[i]]$estimates[1,1]
#   loo_oak_mods3_df$elpd_loo.se[i] <- loo_oak_mods3$loos[[i]]$estimates[1,2]
#   loo_oak_mods3_df$p_loo.estimate[i] <- loo_oak_mods3$loos[[i]]$estimates[2,1]
#   loo_oak_mods3_df$p_loo.se[i] <- loo_oak_mods3$loos[[i]]$estimates[2,2]
#   loo_oak_mods3_df$looic.estimate[i] <- loo_oak_mods3$loos[[i]]$estimates[3,1]
#   loo_oak_mods3_df$looic.se[i] <- loo_oak_mods3$loos[[i]]$estimates[3,2]
# }
# loo_oak_mods3_df <- loo_oak_mods3_df[order(loo_oak_mods3_df$elpd_diff,
#                                            decreasing = TRUE),]
# head(loo_oak_mods3_df)
#
# # write.table(loo_oak_mods3_df, "results/looic_model3_values.tsv",
# #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
# # write.table(loo_oak_mods3$ic_diffs__, "results/looic_model3_comparisons.tsv",
# #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)



## C) COMPARE BEST PERFORMING MODELS FROM EACH SET TO FIND THE SIMPLEST BEST-FIT MODEL
# ## Ask for 100 GB, takes ca 1 day to run
#
# ## i) Using all best performing models from all 3 groups
# load(paste("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/results",
#            "/oak_models_all.RData", sep = ""))
# rm(list=setdiff(ls(), c("oak_mod018", ## out of curiosity
#                         "oak_mod043", "oak_mod044",
#                         "oak_mod045", "oak_mod083", "oak_mod085", "oak_mod084",
#                         "oak_mod065", "oak_mod061",
#                         "oak_mod101", "oak_mod102", "oak_mod103", "oak_mod104",
#                         "oak_mod105", "oak_mod106", "oak_mod108", "oak_mod112",
#                         "oak_mod113", "oak_mod121", "oak_mod122", "oak_mod123",
#                         "oak_mod124", "oak_mod125", "oak_mod126", "oak_mod128",
#                         "oak_mod132")))
#
# ## The "best-performing" models I selected are models for which, in each group:
# ## *Their value overlaps with that of the lowest with the lowest LOOIC in their group,
# ## i.e. 0 not contained in [eldp_diff - se_diff, eldp_diff + se_diff]
# ## *If the value does not overlap, then eldp_diff/se_diff < 2??
# ## See https://discourse.mc-stan.org/t/understanding-looic/13409/3
# ## "'Say [eldp_diff/se_diff] was only twice as big??? would that mean there is a negligible
# ## difference between the models?' 'Not necessarily negligible difference, but certainly
# ## there is non-negligible uncertainty about the difference'": because model 43 might be
# ## performing better than 18, then we keep the more complex model 43 just in case (if
# ## that makes sense, that's my reasoning)
#
# ## Compare models
# print("start")
# loo_oak_mods <- loo(oak_mod018,             ## out of curiosity
#                     oak_mod043, oak_mod044,
#                     oak_mod045, oak_mod061, oak_mod065,
#                     oak_mod083, oak_mod084, oak_mod085,
#                     oak_mod101, oak_mod102, oak_mod103, oak_mod104,
#                     oak_mod105, oak_mod106, oak_mod108, oak_mod112,
#                     oak_mod113, oak_mod121, oak_mod122, oak_mod123,
#                     oak_mod124, oak_mod125, oak_mod126, oak_mod128,
#                     oak_mod132,
#                     save_psis = TRUE, moment_match = TRUE,
#                     reloo = TRUE, compare = TRUE)
# loo_oak_mods$diffs
#
# ## Save image
# # save.image(paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/",
# #                   "results/oak_models_looic.RData"))
# # print("saved")
#
#
# ## ii) Using only the best performing models from all 3 groups that only had variables
# ## with significant effects
# ## Ask for 80 GB, takes ca 1 day to run
# # load(paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/",
# #             "results/oak_models_looic.RData"))
#
# ## The reasoning behind this is that, so, for example - model 45 was one of the best
# ## performing models from group 3, but I am not including it here, because:
# ## Model 43: host.status ~ geo.mean + phylo.min + (1 | phylo)
# ## Model 45: host.status ~ geo.mean + phylo.min + phylo.mean + (1 | phylo)
# ## In model 43, all variables have a significant effect (i.e. they do not overlap with
# ## 0) but in model 45, phylo.mean does not. As it does not perform better than model 43,
# ## and when removing its non-significant variable it is == model 43, I decided to remove
# ## it from this table.
# ## For all models in the final comparison above with non-significant  variables, their
# ## "respective" model with only significant variables is also included above (i.e.,
# ## there is a "model 43" for every "model 45"), so I decided to exclude the models with
# ## non-significant terms.
#
# ## Compare models
# print("start")
# loo_oak_mods_sig <- loo(oak_mod018, ## out of curiosity
#                         oak_mod043, oak_mod044,
#                         oak_mod061, oak_mod065, oak_mod083, oak_mod084,
#                         oak_mod101, oak_mod103, oak_mod108, oak_mod113,
#                         oak_mod121, oak_mod123, oak_mod128,
#                         save_psis = TRUE, moment_match = TRUE,
#                         reloo = TRUE, compare = TRUE)
# loo_oak_mods$diffs
#
# ## Save image
# # save.image(paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/",
# #                   "results/oak_models_looic.RData"))
# # print("saved")



## D) CHECK COMPARISONS
# load(paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/",
#             "results/oak_models_looic.RData"))
#
# ## i) All models
# ## Transform object into DF and add model diff information
# loo_oak_mods_df <- data.frame(matrix(nrow = length(loo_oak_mods$diffs[,1]),
#                                      ncol = 2))
# colnames(loo_oak_mods_df) <- c("elpd_diff", "se_diff")
# rownames(loo_oak_mods_df) <- names(loo_oak_mods$diffs[,1])
# loo_oak_mods_df$elpd_diff <- loo_oak_mods$diffs[,1]
# loo_oak_mods_df$se_diff <- loo_oak_mods$diffs[,2]
# # loo_oak_mods_df$sig <- FALSE
# #
# # ## Check significance
# # for (i in 1:nrow(loo_oak_mods_df)) {
# #   sig1 <- loo_oak_mods_df$elpd_diff[i] + loo_oak_mods_df$se_diff[i]
# #   sig2 <- loo_oak_mods_df$elpd_diff[i] - loo_oak_mods_df$se_diff[i]
# #
# #   if((sign(sig1) == 1 & sign(sig2) == -1) | sig1 == 0 | sig2 == 0) {
# #     loo_oak_mods_df$sig[i] <- TRUE
# #   } else {
# #     break
# #   }
# # }; rm(i, sig1, sig2)
#
# ## Now add model-specific information
# loo_oak_mods_df <-  loo_oak_mods_df[order(rownames(loo_oak_mods_df)),]
# loo_oak_mods_df$elpd_loo.estimate <- NA
# loo_oak_mods_df$elpd_loo.se <- NA
# loo_oak_mods_df$p_loo.estimate <- NA
# loo_oak_mods_df$p_loo.se <- NA
# loo_oak_mods_df$looic.estimate <- NA
# loo_oak_mods_df$looic.se <- NA
#
# for (i in 1:nrow(loo_oak_mods_df)) {
#   loo_oak_mods_df$elpd_loo.estimate[i] <- loo_oak_mods$loos[[i]]$estimates[1,1]
#   loo_oak_mods_df$elpd_loo.se[i] <- loo_oak_mods$loos[[i]]$estimates[1,2]
#   loo_oak_mods_df$p_loo.estimate[i] <- loo_oak_mods$loos[[i]]$estimates[2,1]
#   loo_oak_mods_df$p_loo.se[i] <- loo_oak_mods$loos[[i]]$estimates[2,2]
#   loo_oak_mods_df$looic.estimate[i] <- loo_oak_mods$loos[[i]]$estimates[3,1]
#   loo_oak_mods_df$looic.se[i] <- loo_oak_mods$loos[[i]]$estimates[3,2]
# }
# loo_oak_mods_df <- loo_oak_mods_df[order(loo_oak_mods_df$elpd_diff,
#                                          decreasing = TRUE),]
# head(loo_oak_mods_df)
#
# # write.table(loo_oak_mods_df, "results/looic_model_values.tsv",
# #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
# # write.table(loo_oak_mods$ic_diffs__, "results/looic_model_comparisons.tsv",
# #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
#
#
# ## ii) 'Significant' models
# ## Transform object into DF and add model diff information
# loo_oak_mods_sig_df <- data.frame(matrix(nrow = length(loo_oak_mods_sig$diffs[,1]),
#                                      ncol = 2))
# colnames(loo_oak_mods_sig_df) <- c("elpd_diff", "se_diff")
# rownames(loo_oak_mods_sig_df) <- names(loo_oak_mods_sig$diffs[,1])
# loo_oak_mods_sig_df$elpd_diff <- loo_oak_mods_sig$diffs[,1]
# loo_oak_mods_sig_df$se_diff <- loo_oak_mods_sig$diffs[,2]
# # loo_oak_mods_sig_df$sig <- FALSE
# #
# # ## Check significance
# # for (i in 1:nrow(loo_oak_mods_sig_df)) {
# #   sig1 <- loo_oak_mods_sig_df$elpd_diff[i] + loo_oak_mods_sig_df$se_diff[i]
# #   sig2 <- loo_oak_mods_sig_df$elpd_diff[i] - loo_oak_mods_sig_df$se_diff[i]
# #
# #   if((sign(sig1) == 1 & sign(sig2) == -1) | sig1 == 0 | sig2 == 0) {
# #     loo_oak_mods_sig_df$sig[i] <- TRUE
# #   } else {
# #     break
# #   }
# # }; rm(i, sig1, sig2)
#
# ## Now add model-specific information
# loo_oak_mods_sig_df <-  loo_oak_mods_sig_df[order(rownames(loo_oak_mods_sig_df)),]
# loo_oak_mods_sig_df$elpd_loo.estimate <- NA
# loo_oak_mods_sig_df$elpd_loo.se <- NA
# loo_oak_mods_sig_df$p_loo.estimate <- NA
# loo_oak_mods_sig_df$p_loo.se <- NA
# loo_oak_mods_sig_df$looic.estimate <- NA
# loo_oak_mods_sig_df$looic.se <- NA
#
# for (i in 1:nrow(loo_oak_mods_sig_df)) {
#   loo_oak_mods_sig_df$elpd_loo.estimate[i] <- loo_oak_mods_sig$loos[[i]]$estimates[1,1]
#   loo_oak_mods_sig_df$elpd_loo.se[i] <- loo_oak_mods_sig$loos[[i]]$estimates[1,2]
#   loo_oak_mods_sig_df$p_loo.estimate[i] <- loo_oak_mods_sig$loos[[i]]$estimates[2,1]
#   loo_oak_mods_sig_df$p_loo.se[i] <- loo_oak_mods_sig$loos[[i]]$estimates[2,2]
#   loo_oak_mods_sig_df$looic.estimate[i] <- loo_oak_mods_sig$loos[[i]]$estimates[3,1]
#   loo_oak_mods_sig_df$looic.se[i] <- loo_oak_mods_sig$loos[[i]]$estimates[3,2]
# }
# loo_oak_mods_sig_df <- loo_oak_mods_sig_df[order(loo_oak_mods_sig_df$elpd_diff,
#                                          decreasing = TRUE),]
# head(loo_oak_mods_sig_df)
#
# # write.table(loo_oak_mods_sig_df, "results/looic_model_sig_values.tsv",
# #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
# # write.table(loo_oak_mods_sig$ic_diffs__, "results/looic_model_sig_comparisons.tsv",
# #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


#
##### VIII) EXPLORE THE BEST-FITTING MODELS                                          #####
## Only keep best-performing models + other models we wish to explore
# load(paste("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/",
#            "/results/oak_models_all.RData", sep = ""))
# rm(list=setdiff(ls(), c("oak_mod002", "oak_mod004", "oak_mod008", "oak_mod010",
#                         "oak_mod018", "oak_mod025")))
# load(paste("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/",
#            "/results/oak_models_looic.RData", sep = ""))
# rm(l, loo_oak_mods, loo_oak_mods_sig, i,
#    oak_mod045, oak_mod085, oak_mod102, oak_mod104, oak_mod105, oak_mod106,
#    oak_mod112, oak_mod122, oak_mod124, oak_mod125, oak_mod126, oak_mod132)
# save.image(paste("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/",
#                  "/results/oak_models_final.RData", sep = ""))
load(paste("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/",
           "/results/oak_models_final.RData", sep = ""))

## A) EXPLORE BEST-PERFORMING MODELS
## i) oak_mod043
oak_mod043$formula

## Plot
# svg("results/mcmc_plot_oak_mod043.svg")
mcmc_plot(oak_mod043) +
  geom_vline(xintercept = 0, linetype="dashed", colour = "maroon") +
  ggtitle(expression(bold("Posterior distributions"))) +
  theme(text = element_text(family = "Arial", size = 14),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  theme_bw()
# dev.off()

## Summary
## NB, see https://stats.stackexchange.com/questions/573042/how-is-the-standard-
## deviation-of-random-effects-estimated
## And https://stats.stackexchange.com/questions/238005/what-is-the-intuition-on-fixed
## -and-random-effects-models
## And maybe https://www.researchgate.net/post/How-do-I-report-the-results-of-a-linear-
## mixed-models-analysis
summary(oak_mod043)
head(ranef(oak_mod043)$quercus.sp[,,1])
# write.table(ranef(oak_mod043)$quercus.sp[,,1], "results/mod43_ranef.tsv",
#             se= "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


## ii) oak_mod065
oak_mod065$formula

## Plot
mcmc_plot(oak_mod065) +
  geom_vline(xintercept = 0, linetype="dashed", colour = "maroon") +
  ggtitle(expression(bold("Posterior distributions"))) +
  theme(text = element_text(family = "Arial", size = 14),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  theme_bw()

## Summary
## NB, see https://stats.stackexchange.com/questions/573042/how-is-the-standard-
## deviation-of-random-effects-estimated
## And https://stats.stackexchange.com/questions/238005/what-is-the-intuition-on-fixed
## -and-random-effects-models
## And maybe https://www.researchgate.net/post/How-do-I-report-the-results-of-a-linear-
## mixed-models-analysis
summary(oak_mod065)
# ranef(oak_mod065)$quercus.sp[,,1]
# ranef(oak_mod065)$agrilus.sp[,,1]

hist(ranef(oak_mod065)$quercus.sp[,,1][,1])
# hist(ranef(oak_mod065)$quercus.sp[,,1][,2])

# svg("results/oak_model065_Agrilus_random_effect_estimates.svg")
hist(ranef(oak_mod065)$agrilus.sp[,,1][,1],
     breaks = 6,
     col = "darkslategray4",
     xlim = c(-0.7, 0.5),
     ylim = c(0, 20),
     xlab = substitute(paste("Estimate of the intercept value for the random ",
                             "effect for each ", italic(Agrilus), " species")),
     ylab = "Frequency",
     main = "Model 65")
# dev.off()
# hist(ranef(oak_mod065)$agrilus.sp[,,1][,2])



## B) EXPLORE THE SIMPLEST BEST-FITTING MODEL (SBFM) IN DEPTH (oak_mod047)
## This is a binary logistic regression (link: mu = logit), i.e., the response ==
## log-odds of something being a host. To convert it to probability: exp(y)/(1+exp(y))
oak_mod043$formula

## NB, to chech variable names: variables(oak_mod043)[1:5]

## i) First, have a look at the summary
summary(oak_mod043)

## Intercept: log-odds of the whole population being a host with no predictor variables.
## To convert it to probability:
exp(-3.16)/(1+exp(-3.16))*100

## Decrease in % prob of being a host per unit increase in norm.geo.dist (I think):
(exp(-3.16)/(1+exp(-3.16))*100) - (exp(-3.16-3.61*1)/(1+exp(-3.16-3.61*1))*100)

## Decrease in % prob of being a host per unit increase in phylo.dist.min (I think):
(exp(-3.16)/(1+exp(-3.16))*100) - (exp(-3.16-1.51*1)/(1+exp(-3.16-1.51*1))*100)


## ii) Diagnostics
## Compare the observed outcome variable y to simulated datasets yrep from the
## posterior predictive distribution (looks good).
## yrep: predicted scores from single replication; y: observed outcome
pp_check(oak_mod043)                ## looks OK I think
pp_check(oak_mod043, "stat")        ## mc-stan.org/rstanarm/reference/

## Diagnostic plots (Inspect chains - Trace and Density Plots)
## b: population-level
plot(oak_mod043, N = 2, ask = FALSE)


## iii) Predict interaction status (as log-odds) for each Q-A pair:
## Note:
## *fitted(, scale = ""response") (== predict()) uses the response scale, i.e., results
## after applying the inverse link function
## *fitted(, scale = "linear") return the value pre inverse link transformation
pred_prob <- fitted(oak_mod043, scale = "response")
pred_lodds <- fitted(oak_mod043, scale = "linear")

## If something was next to sth else both in phylogeny and geography
## logit(p) = intercept + random_effect + slope1*fixed_effect1 + slope2*fixed_effect2
r_oak_mod043 <- ranef(oak_mod043)$quercus.sp[,,1][,1]
f_oak_mod043 <- summary(oak_mod043)$fixed$Estimate
f_oak_mod043[1] + max(r_oak_mod043) + f_oak_mod043[2]*0 + f_oak_mod043[3]*0
## And the max result we get is
max(pred_lodds[,1])
which.max(r_oak_mod043)
interaction_data[which(pred_lodds[,1] %in%
                       sort(pred_lodds[,1], decreasing = TRUE)[1:4]), c(1:3, 5, 12)]
interaction_data[interaction_data$quercus.sp == "Q. scytophylla" &
                 interaction_data$agrilus.sp == "A. auroguttatus",]

## If something was super far from sth else both in phylogeny and geography
f_oak_mod043[1] + min(r_oak_mod043) + f_oak_mod043[2]*1  + f_oak_mod043[3]*1
## And the min result we get is (literally the same!!) - and it makes sense
min(pred_lodds[,1])
which.min(r_oak_mod043)
interaction_data[which.min(pred_lodds[,1]),]

## If something was next to sth else in phylogeny but not geography
## (Value still over cutoff[2] from next section)
f_oak_mod043[1] + max(r_oak_mod043) + f_oak_mod043[2]*1 + f_oak_mod043[3]*0

## If something was next to sth else in geography but not phylogeny
## (Value under cutoff[2] from next section)
f_oak_mod043[1] + min(r_oak_mod043) + f_oak_mod043[2]*0 + f_oak_mod043[3]*1
rm(r_oak_mod043, f_oak_mod043)

## Understanding probability vs. log-odds
plot(pred_prob[,1] ~ pred_lodds[,1])
plot((exp(pred_lodds[,1])/(1+exp(pred_lodds[,1]))) ~ pred_prob[,1])
rm(pred_prob, pred_lodds)

## Check the effect of the different variables in our model
## NB, fitted(oak_mod10, s = "l") == fitted(oak_mod043, re_formula = NA, s = "l")
comparisons <- data.frame(
  host.status = as.factor(interaction_data$interaction),
  oak = as.factor(interaction_data$quercus.sp),
  mod043 = fitted(oak_mod043, scale = "linear")[,1],  ## full
  mod002 = fitted(oak_mod002, scale = "linear")[,1],  ## geo
  mod004 = fitted(oak_mod004, scale = "linear")[,1],  ## phy1
  mod008 = fitted(oak_mod008, scale = "linear")[,1],  ## (1 | gr(oak | phy2))
  mod010 = fitted(oak_mod010, scale = "linear")[,1],  ## geo + phy1
  mod018 = fitted(oak_mod018, scale = "linear")[,1],  ## geo + (1 | gr(oak | phy2))
  mod025 = fitted(oak_mod025, scale = "linear")[,1]   ## phy1 + (1 | gr(oak | phy2))
  )

## Full model vs. geographic-distance-only model
summary(lm(comparisons$mod043 ~ comparisons$mod002))

# svg("results/full_vs_geo.svg")
ggplot(data = comparisons,
  aes(y = mod043, x = mod002, col = host.status)) +
  scale_colour_manual(values = c("lightgrey", "steelblue")) +
  geom_point() +
  geom_abline(intercept = -0.70436, slope = 1.24781,
              colour = "black", linewidth = 2, linetype = 'dashed') +
  ylab("Full model (model 43)") +
  xlab("Geographic distance model (model 2)") +
  theme_bw() +
  theme(legend.position = "none")
# dev.off()
# hist(interaction_data$min.dist.mean.norm, breaks = 100, ylim = c(0, 500))

## Full model vs. phylogenetic-distance-(fixed)-only model
summary(lm(comparisons$mod043 ~ comparisons$mod004))

# svg("results/full_vs_phylofx.svg")
ggplot(data = comparisons,
       aes(y = mod043, x = mod004, col = host.status)) +
  scale_colour_manual(values = c("lightgrey", "steelblue")) +
  geom_point() +
  geom_abline(intercept = -0.77562, slope = 1.44191,
              colour = "black", linewidth = 2, linetype = 'dashed') +
  ylab("Full model (model 43)") +
  xlab("Phylogenetic distance (fixed) model (model 4)") +
  theme_bw() +
  theme(legend.position = "none")
# dev.off()
# hist(interaction_data$phylo.dist.min, breaks = 100, ylim = c(0, 500))

## Full model vs. phylogenetic-random-only model
summary(lm(comparisons$mod043 ~ comparisons$mod008))

# svg("results/full_vs_phylord.svg")
ggplot(data = comparisons,
       aes(y = mod043, x = mod008, col = host.status)) +
  scale_colour_manual(values = c("lightgrey", "steelblue")) +
  geom_point() +
  geom_abline(intercept = -1.30016, slope = 1.00775,
              colour = "black", linewidth = 2, linetype = 'dashed') +
  ylab("Full model (model 43)") +
  xlab("Phylogenetic random model (model 8)") +
  theme_bw() +
  theme(legend.position = "none")
# dev.off()

## Full model vs. phylogenetic-only (both fixed + random) model
summary(lm(comparisons$mod043 ~ comparisons$mod025))

# svg("results/full_vs_phylo.svg")
ggplot(data = comparisons,
       aes(y = mod043, x = mod025, col = host.status)) +
  scale_colour_manual(values = c("lightgrey", "steelblue")) +
  geom_point() +
  geom_abline(intercept = -1.093212, slope = 0.924603,
              colour = "black", linewidth = 2, linetype = 'dashed') +
  ylab("Full model (model 43)") +
  xlab(paste("Phylogenetic-only model\n",
             "(model 25, both fixed + random phylogenetic effects)")) +
  theme_bw() +
  theme(legend.position = "none")
# dev.off()

## Full model vs. geographic + phylogenetic (fixed) distance model
summary(lm(comparisons$mod043 ~ comparisons$mod010))

# svg("results/full_vs_geophylofx.svg")
## Clustered by oak species
ggplot(data = comparisons,
       aes(y = mod043, x = mod010, col = host.status)) +
  scale_colour_manual(values = c("lightgrey", "steelblue")) +
  geom_point() +
  geom_abline(intercept = -0.24889, slope = 1.30971,
              colour = "black", linewidth = 2, linetype = 'dashed') +
  ylab("Full model (model 43)") +
  xlab(paste("Geographic and phylogenetic (fixed) distance model\n",
             "(model 10, no phylogenetic random effect)")) +
  theme_bw() +
  theme(legend.position = "none")
# dev.off()

## Full model vs. geographic distance + phylogenetic random effect model
summary(lm(comparisons$mod043 ~ comparisons$mod018))

# svg("results/full_vs_geophylord.svg")
ggplot(data = comparisons,
       aes(y = mod043, x = mod018, col = host.status)) +
  scale_colour_manual(values = c("lightgrey", "steelblue")) +
  geom_point() +
  geom_abline(intercept = -1.093212, slope = 0.924603,
              colour = "black", linewidth = 2, linetype = 'dashed') +
  ylab("Full model (model 43)") +
  xlab(paste("Geographic distance and phylogenetic-random-effect model\n",
             "(model 10, no phylogenetic distance fixed effect)")) +
  theme_bw() +
  theme(legend.position = "none")
# dev.off()
rm(comparisons)

# ## Explore interaction ~ (1 | gr(oak | phylo)) (oak_mod008)
# plot(fitted(oak_mod008, scale = "linear")[,1])
# abline(-6.35, 0, col = "blue", lwd = 2)
# abline(a = (-6.33+sd(fitted(oak_mod008, scale = "linear")[,1])),
#        b = 0, col = "red", lwd = 2)
# abline(a = (-6.33-sd(fitted(oak_mod008, scale = "linear")[,1])),
#        b = 0, col = "red", lwd = 2)

## Save model predictions into a DF
predictions <- data.frame(cbind(interaction_data$quercus.sp,
                                interaction_data$agrilus.sp,
                                interaction_data$interaction,
                                fitted(oak_mod043, scale = "response")[,1],
                                fitted(oak_mod043, scale = "linear")[,1]))
colnames(predictions) <- c("quercus.sp", "agrilus.sp", "host.status",
                           "prediction.prob", "prediction.logodds")
predictions$host.status <- as.factor(predictions$host.status)
predictions$prediction.prob <- as.numeric(predictions$prediction.prob)
predictions$prediction.logodds <- as.numeric(predictions$prediction.logodds)

## Save data
# write.table(predictions,
#             "results/preditions_all_quercus.tsv",
#             sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

## Make a very simple plot
# plot(predictions$prediction.prob ~ predictions$host.status)
# plot(predictions$prediction.logodds ~ predictions$host.status)



#
##### IX)   FIND A PREDICTION THRESHOLD FOR THE SBFM AND PLOT PREDICTIONS            #####
## A) COMPUTE THRESHOLD
## i) Compute range of thresholds to test
## First, find the minimum and maximum predicted values
## Minimum predicted value (rounded to closest lower .0 or 0.5 decimal)
min_v <- min(predictions$prediction.logodds); min_v

## If round == floor, leave as round(); else, leave as round() - 0.5
if ((round(min_v)) == floor(min_v)) {
  min_v <- round(min_v)
} else {
  min_v <- round(min_v) - 0.5
}; min_v

## Maximum predicted value (rounded to closest upper .0 or 0.5 decimal)
max_v <- max(predictions$prediction.logodds); max_v

## If round == floor, leave as round() + 0.5; else, leave as round()
if (round(max_v) == floor(max_v)) {
  max_v <- round(max_v) + 0.5
} else {
  max_v <- round(max_v)
}; max_v

## Next, get a sequence of the threshold values we're going to be looking at (i.e.,
## values from the min. val to the max. val, with 0.5 increases)
thrs <- seq(from = min_v, to = max_v, by = 0.5); rm(min_v, max_v)


## ii) Compute threshold value
## See https://stackoverflow.com/questions/23240182/deciding-threshold-for-glm-logis
## tic-regression-model-in-r
## And https://statinfer.com/203-4-2-calculating-sensitivity-and-specificity-in-r/

## Initialise an empty vector for:
sens <- c()            ## Sensitivity (i.e., True Positive Rate = TP/(TP+FN))
spec <- c()            ## Specificity (i.e., True Negative Rate = TN/(TN+FP))
fpr <- c()             ## FPR (i.e., False Positive Rate = FP/(FP + TN))

## Now, for every threshold value
for (thr in thrs) {
  ## First, compute the Confusion Matrix (i.e., "a summary of prediction results on a
  ## classification problem") under current threshold
  ## Rows = Predicted values classified as 1s or 0s according to the current threshold
  ## Cols = Observed values
  ## Matrix will look like:
  ##            0 (Obs)   1 (Obs)
  ## 0 (Pred)   TN        FN
  ## 1 (pred)   FP        TP
  conf_matrix <- table(ifelse(predictions$prediction.logodds > thr, 1, 0),   ## rows
                       predictions$host.status)                              ## cols
  # print(conf_matrix)

  ## Then, compute sensitivity (TPR), specificity (TNR) and FPR
  ## Let's deal with the exceptions first (i.e., if all predicted values are either 1
  ## or 0 under the current threshold)
  if(nrow(conf_matrix) == 1) {
    ## If all pred vals are == 1, then TPR (sens) and FPR = 1, and TNR (spec) = 0
    if(dimnames(conf_matrix)[[1]] == "1") {
      sens <- c(sens, 1); spec <- c(spec, 0); fpr <- c(fpr, 1)
      ## If all pred vals are == 0, then TPR (sens) and FPR = 0, and TNR (spec) = 1
    } else {
      sens <- c(sens, 0); spec <- c(spec, 1); fpr <- c(fpr, 0)
    }

    ## Now, if the current threshold divides predicted values into "1"s and "0"s
  } else {
    ## Sensitivity (i.e., True Positive Rate) = TP/(TP+FN)
    sens <- c(sens, conf_matrix[2,2]/(conf_matrix[2,2] + conf_matrix[1,2]))

    ## Specificity (i.e., True Negative Rate) = TN/(TN+FP)
    spec <- c(spec, conf_matrix[1,1]/(conf_matrix[1,1] + conf_matrix[2,1]))

    ## False Positive Rate (FPR) = FP/(FP + TN)
    fpr <- c(fpr, conf_matrix[2,1]/(conf_matrix[2,1] + conf_matrix[1,1]))
  }
}

## Create a DF with all this info
thrs_df <- as.data.frame(cbind(thrs, sens, spec, fpr))
rm(sens, spec, thrs, fpr, thr, conf_matrix)


## iv) Plot results and find intercept
## First, plot sensitivity and specificity vs. threshold, to get a visual idea of where
## the cutoff point is
thrs_df %>% ggplot(aes(thrs, sens)) +
  geom_line(colour = "steelblue", linewidth = 2) +
  geom_line(aes(thrs, spec, col = "coral"), linewidth = 2) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Specificity")) +
  # xlim(-4.5, -4) + ylim(0.875, 0.925) +
  labs(x = "Threshold", y = "Sensitivity") +
  theme(axis.title.y = element_text(colour = "steelblue"),
        axis.title.y.right = element_text(colour = "coral"),
        axis.title = element_text(size = 12), # face = "bold"
        legend.position = "none")

## Now, find the intersect (I'm sure there's a better, more precise and sophisticated
## way to go about this...)
## We saw above that the cutoff point is between -4.5 and -4.0 (as this is where the
## lines intersect), so now we can compute linear models for both the sensitivity and
## the specificity lines between these two points
lm_sens <- lm(thrs_df[thrs_df$thrs == -4.5 | thrs_df$thrs == -4.0, c(1)] ~
                thrs_df[thrs_df$thrs == -4.5 | thrs_df$thrs == -4.0, c(2)])
lm_spec <- lm(thrs_df[thrs_df$thrs == -4.5 | thrs_df$thrs == -4.0, c(1)] ~
                thrs_df[thrs_df$thrs == -4.5 | thrs_df$thrs == -4.0, c(3)])

## Compute the intercept
## See https://stackoverflow.com/questions/7114703/finding-where-two-linear-fits-inter
## sect-in-r
cm <- rbind(coef(lm_sens),coef(lm_spec))    ## Coefficient matrix
cutoff <- c(-solve(cbind(cm[,2],-1))
            %*% cm[,1])
cutoff                                      ## y (sens|spec) = 0.905; x (thr) = -4.154


## And finally, plot properly with cutoff points
# svg("results/cutoff.svg")
thrs_df %>% ggplot(aes(thrs, sens)) +
  geom_line(colour = "steelblue", linewidth = 2) +
  geom_line(aes(thrs, spec, col = "coral"), linewidth = 2) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Specificity")) +
  geom_hline(yintercept = cutoff[1], linetype="dashed", color = "darkgrey") +
  geom_vline(xintercept = cutoff[2], linetype="dashed", color = "darkgrey") +
  labs(x = "Threshold", y = "Sensitivity") +
  theme(axis.title.y = element_text(colour = "steelblue"),
        axis.title.y.right = element_text(colour = "coral"),
        axis.title = element_text(size = 12), # face = "bold"
        legend.position = "none")
# dev.off()

## To end with, also have a look at the AUC plot, which can be used to assess the
## ability of the model to differentiate between known hosts and alleged non-hosts
## under different probability thresholds. "A value of 0.8 for the AUC means that
## for 80% of the time a random selection from the positive group will have a score
## greater than a random selection from the negative class"
## ROC: plot(thrs_df$sens~c(1-thrs_df$spec))
thrs_df %>% ggplot(aes(fpr, sens)) +
  geom_line() +
  labs(x = 'FPR', y = "TPR") +
  theme(axis.title.y.right = element_text(colour = "red"), legend.position="none")

## Compute AUC value
DescTools::AUC(thrs_df$fpr, thrs_df$sens)    ## 0.963
rm(thrs_df, lm_spec, lm_sens, cm)



## B) PLOT RESULTS
## i) Violin plot
# svg("results/preditions_violin_general_oak_mod043.svg")
ggplot(data = predictions,
       aes(y = prediction.logodds, x = host.status, fill = host.status)) +

  ## Threshold using 'baseline' log-odds of sth being a host
  ## NB, 7552 = no. observations & 116 == no. 'positive' observations
  ## logit(p) = log(p/(1-p)), with p = prob
  # geom_hline(yintercept = log((116/7552)/(1-(116/7552))),
  #            linetype = "dashed", col = "black") +
  ## New inferred threshold (very close to the one above!)
  geom_hline(yintercept = cutoff[2], linewidth = 1.5,
             linetype = "dotted", col = "maroon") +

  geom_violin(width = 1) +
  geom_boxplot(width = 0.1, fill = "grey") +

  scale_fill_manual(values = c("aquamarine2", "coral")) +
  ggtitle("Predicted interactions") +
  xlab("Host status") + ylab("Log-odds predicted host status") +

  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))
# dev.off()

median(predictions[predictions$host.status == 1, ]$prediction.logodds)
median(predictions[predictions$host.status == 0, ]$prediction.logodds)


## ii) Binary plots using threshold
## Prepare data
predictions <- data.frame(predictions,
                          ifelse(predictions$prediction.logodds > cutoff[2], 1, 0))
colnames(predictions) <- c(colnames(predictions)[1:5], "prediction.cutoff")
predictions$prediction.cutoff <- as.factor(predictions$prediction.cutoff)

## Check confidence matrix (predictions = rows, observations = columns)
conf_matrix <- table(ifelse(predictions$prediction.logodds > cutoff[2], 1, 0),
                     predictions$host.status)
## % true positives: 90.52%
round(conf_matrix[2,2]/(conf_matrix[1,2] + conf_matrix[2,2])*100, 2)
## % true negatives: 90.64%
round(conf_matrix[1,1]/(conf_matrix[1,1] + conf_matrix[2,1])*100, 2)

## Plot confusion matrix as a fourfold plot
## See https://stackoverflow.com/questions/23891140/r-how-to-visualize-confusion-ma
## trix-using-the-caret-package
##            0 (Obs)   1 (Obs)
## 0 (Pred)   TN        FN
## 1 (pred)   FP        TP
# svg("results/preditions_binary_general_oak_mod043.svg")
fourfoldplot(conf_matrix, color = c("#CC6666", "#99CC99"),
             conf.level = 0, main = "Confusion Matrix",
             margin = c(2))       ## standardizing the col
# dev.off()

# write.table(subset(predictions,
#                    prediction.logodds >= cutoff[2] & host.status == 0)[, c(1,2,5)],
#             "results/potential_future_hosts.tsv",
#             sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# ## Plot confusion matrix as a heatmap
# cm <- yardstick::conf_mat(predictions, host.status, prediction.cutoff)
# autoplot(cm, type = "heatmap") +
#   scale_fill_gradient(low = "#769df5", high = "#1e3e7d") +
#   xlab("Reported status") +
#   ylab("Predicted status")
# rm(cm)



##### X)    PLOT PREDICTIONS BY OAK/AGRILUS SPECIES                                  #####
## A) PER OAK SPECIES
## i) All oak species
# svg("results/predictions_oaks.svg",width = 20, height = 20)
ggplot(predictions,
       aes(x = host.status, y = prediction.logodds)) +
  geom_point() +
  geom_hline(yintercept = cutoff[2], linewidth = 1.5,
             linetype = "dotted", col = "maroon") +
  facet_wrap(~ quercus.sp, ncol = 16) +
  theme(strip.text = element_text(face = "italic"))
# dev.off()


## ii) Only hosts
hosts <- unique(subset(predictions, host.status == 1)$quercus.sp)
hosts <- subset(predictions, quercus.sp %in% hosts)

ggplot(hosts, aes(x = host.status, y = prediction.logodds)) +
  geom_point() +
  geom_hline(yintercept = cutoff[2], linewidth = 1.5,
             linetype = "dotted", col = "maroon") +
  facet_wrap(~ quercus.sp, ncol = 10)
which.max(coefficients(oak_mod043)$quercus.sp[1:236])        ## Q. acutissima, -3.97



## B) PER AGRILUS SPECIES
## i) For all individual Agrilus spp.
# svg("results/predictions_agrilus.svg",width = 20, height = 10)
ggplot(predictions,
       aes(x = host.status, y = prediction.logodds)) +
  geom_hline(yintercept = cutoff[2], linewidth = 1.5,
             linetype = "dotted", col = "maroon") +
  geom_point() + facet_wrap(~ agrilus.sp, ncol = 8) +
  theme(strip.text = element_text(face = "italic"))
# dev.off()

## Check whether Q. robur and Q. rubra are always on top
# predictions$colour <- "black"
# predictions[predictions$quercus.sp %in%
#             c("Q. robur", "Q. rubra"),]$colour <- rep("maroon", 64)
# ggplot(predictions, aes(x = host.status, y = prediction.logodds, col = colour)) +
#   geom_point() + facet_wrap(~ agrilus.sp, ncol = 8) + theme(legend.position = "none")
# predictions$colour <- NULL


## ii) For a particular Agrilus sp.
## Continuous plot
q <- subset(predictions, agrilus.sp == "A. angustulus" & host.status == 0 &
              prediction.logodds >= cutoff[2])$quercus.sp

# svg(file = "a_angustulus_preds.svg", width = 40, height = 25)
ggplot(subset(predictions, agrilus.sp == "A. angustulus"),
       aes(x = host.status, y = prediction.logodds)) + ## , fill = host.status))
  # geom_violin(width = 1) +
  # geom_boxplot(width = 0.1, fill = "grey") +
  # scale_fill_manual(values = c("aquamarine2", "coral")) +
  geom_point(size = 5) +
  facet_wrap(~ agrilus.sp, ncol = 8) +
  geom_hline(yintercept = cutoff[2],
             linetype = "dashed", color = "red") +
  ggrepel::geom_label_repel(aes(label = ifelse(quercus.sp %in% q,
                                               as.character(quercus.sp), ''),),
                            hjust = 0, vjust = 0, max.overlaps = 100, size = 4)
# dev.off()

## Binary plot
## Check confidence matrix (predictions = rows, observations = columns)
conf_matrix <- table(ifelse(subset(predictions,
                                   agrilus.sp == "A. biguttatus")$prediction.logodds >
                              cutoff[2], 1, 0),
                     subset(predictions,
                            agrilus.sp == "A. biguttatus")$host.status)
## % true positives: 100%
round(conf_matrix[2,2]/(conf_matrix[1,2] + conf_matrix[2,2])*100, 2)
## % true negatives: 89.91%
round(conf_matrix[1,1]/(conf_matrix[1,1] + conf_matrix[2,1])*100, 2)

## Plot confusion matrix as a fourfold plot
## See https://stackoverflow.com/questions/23891140/r-how-to-visualize-confusion-ma
## trix-using-the-caret-package
##            0 (Obs)   1 (Obs)
## 0 (Pred)   TN        FN
## 1 (pred)   FP        TP
# svg(file = "a_biguttatus_preds_binary.svg", width = 40, height = 25)
fourfoldplot(conf_matrix, color = c("#CC6666", "#99CC99"),
             conf.level = 0, main = "Confusion Matrix",
             margin = c(2))
# dev.off()


## iii) Agrilus auroguttatus (native to AZ in USA, introduced to CA in USA)
q <- c("Q. kelloggii", "Q. agrifolia", "Q. engelmannii", "Q. chrysolepis")
predictions$colour <- "black"
predictions[predictions$quercus.sp %in% q &
            predictions$agrilus.sp == "A. auroguttatus",]$colour <- rep("maroon",
                                                                        length(q))

# svg(file = paste0(pw, "results/a_aur_preds.svg"))
ggplot(subset(predictions, agrilus.sp == "A. auroguttatus"),
       aes(x = host.status, y = prediction.logodds, col = colour)) +
  scale_color_manual(values = c("black", "maroon")) +
  geom_point(size = 2) +
  facet_wrap(~ agrilus.sp, ncol = 8) +
  geom_hline(yintercept = cutoff[2],
             linetype = "dashed", color = "maroon") +
  ggrepel::geom_label_repel(aes(label = ifelse(quercus.sp %in% q,
                                               as.character(quercus.sp), ''),),
                            col ="black", size = 5,
                            hjust = 1.5, vjust = 0, max.overlaps = 100) +
  theme_bw() +
  theme(legend.position="none")
# dev.off()

predictions$colour <- NULL

cutoff[2]
subset(predictions, agrilus.sp == "A. auroguttatus" &
                    quercus.sp %in% c("Q. kelloggii", "Q. agrifolia",
                                      "Q. engelmannii"))[,c(1,2,5,6)]


## iii) Agrilus bilineatus
q <- c("Q. robur")
predictions$colour <- "black"
predictions[predictions$quercus.sp %in% q &
              predictions$agrilus.sp == "A. bilineatus",]$colour <- rep("maroon",
                                                                        length(q))

# svg(file = "results/a_bil_preds.svg")
ggplot(subset(predictions, agrilus.sp == "A. bilineatus"),
       aes(x = host.status, y = prediction.logodds, col = colour)) +
  scale_color_manual(values = c("black", "maroon")) +
  geom_point(size = 2) +
  facet_wrap(~ agrilus.sp, ncol = 8) +
  geom_hline(yintercept = cutoff[2],
             linetype = "dashed", color = "maroon") +
  ggrepel::geom_label_repel(aes(label = ifelse(quercus.sp %in% q,
                                               as.character(quercus.sp), ''),),
                            col ="black", size = 5,
                            hjust = 1.5, vjust = 0, max.overlaps = 100) +
  theme_bw() +
  theme(legend.position="none")
# dev.off()

predictions$colour <- NULL


## v) Agrilus sulcicollis
q <- c("Q. rubra")
predictions$colour <- "black"
predictions[predictions$quercus.sp %in% q &
              predictions$agrilus.sp == "A. sulcicollis",]$colour <- rep("maroon",
                                                                        length(q))

# svg(file = "results/a_sul_preds.svg")
ggplot(subset(predictions, agrilus.sp == "A. sulcicollis"),
       aes(x = host.status, y = prediction.logodds, col = colour)) +
  scale_color_manual(values = c("black", "maroon")) +
  geom_point(size = 2) +
  facet_wrap(~ agrilus.sp, ncol = 8) +
  geom_hline(yintercept = cutoff[2],
             linetype = "dashed", color = "maroon") +
  ggrepel::geom_label_repel(aes(label = ifelse(quercus.sp %in% q,
                                               as.character(quercus.sp), ''),),
                            col ="black", size = 5,
                            hjust = 1.5, vjust = 0, max.overlaps = 100) +
  theme_bw() +
  theme(legend.position="none")
# dev.off()

predictions$colour <- NULL



#
##### XI)   LEAVE-ONE-OUT CHECKS                                                     #####
## A) FIRST WAY - CODE ONE INTERACTION AS A '0' AT A TIME
# ## This modified (using 'fitted' instead of 'predict') piece of code comes from:
# ## https://www.r-bloggers.com/2015/12/calculate-leave-one-out-prediction-for-glm/
# ## Ask for 240 h, 1 core [it only does 1 chain at a time], 15 GB, -t 1-116; takes ca.
# ## 2h to run.
# n <- as.numeric(commandArgs(trailingOnly = TRUE))
# obs <- which(interaction_data$interaction == 1)
#
# loo_fitted <- function(mod, obs, n) {
#   ndat <- mod$data; i <- obs[n]
#   print(interaction_data[i, 1:3])
#
#   ## If the current interaction == 1, then update the model by making it == 0,
#   ## and predict the value for the current observation
#   ndat$interaction[i] <- 0
#   return(data.frame(row.no = i,
#                     prediction.logodds.loo = fitted(update(mod, newdata = ndat),
#                                                     mod$data[i, ],
#                                                     scale = "linear")[, 1],
#                     row.names = NULL))
# }
#
# predictions_loo0 <- loo_fitted(oak_mod043, obs, n)
#
# # write.table(predictions_loo0,
# #             paste("results/loo/predictions_loo0", n, ".tsv", sep = ""),
# #             row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
#
# ## Make everything into one file (once everything has finished running)
# loo0_pred_loods <- lapply(paste0("results/loo/", list.files(path = "results/loo")),
#                           read.table, header = TRUE, sep = "\t")
# loo0_pred_loods <- do.call(rbind.data.frame, loo0_pred_loods)
# loo0_pred_loods <- loo0_pred_loods[order(loo0_pred_loods$row.no),]
#
# # write.table(loo0_pred_loods, "results/predictions_loo0.tsv",
# #             row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
#
#
## ii) Compare results against original predictions
loo0_pred_loods <- read.table("results/predictions_loo0.tsv",
                              sep = "\t", header = T)$prediction.logodds.loo
predictions1 <- subset(predictions, host.status == 1)
interaction_data1 <- subset(interaction_data, interaction == 1)
predictions_loo1 <- data.frame(agrilus.sp = interaction_data1$agrilus.sp,
                               quercus.sp = interaction_data1$quercus.sp,
                               interaction = as.factor(interaction_data1$interaction),
                               prediction.logodds.orig = predictions1$prediction.logodds,
                               prediction.logodds.loo0 = loo0_pred_loods)

## Check the differences between the "loo0" and the "original" predictions
plot(loo0_pred_loods ~ predictions1$prediction.logodds)
abline(a = 0, b = 1, col = "steelblue")
# identify(y = loo0_pred_loods, x = predictions$prediction.logodds)


## iii) Plot results (only for '1's)
## LOO0 vs. original predicitons
predictions_loo1$thr <- "no"
for (i in 1:nrow(predictions_loo1)) {
  if (predictions_loo1$prediction.logodds.orig[i] >= cutoff[2] &
      predictions_loo1$prediction.logodds.loo0[i] < cutoff[2]) {
    predictions_loo1$thr[i] <- "yes"
  }
}

# svg("results/loo0vsorigpreds.svg")
ggplot(predictions_loo1,
       aes(x = prediction.logodds.orig, y = prediction.logodds.loo0,
           col = thr)) +
  geom_abline(linetype = "dashed", col = "darkgrey", linewidth = 1.5) +
  geom_point() +
  scale_color_manual(values = c("black", "coral")) +
  # geom_hline(yintercept = cutoff[2], linewidth = 1.5,
  #            linetype = "dashed", col = "coral") +
  # geom_vline(xintercept = cutoff[2], linewidth = 1.5,
  #            linetype = "dashed", col = "coral") +
  xlab("Original predictions") + ylab("'One 0 in turn' predictions") +
  ggtitle("Quercus sp. - Agrilus sp.") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")
# dev.off()

predictions_loo1[predictions_loo1$thr == "yes",]
predictions_loo1$thr <- NULL

## Violin plot
ggplot(data = predictions_loo1,
       aes(y = prediction.logodds.loo0, x = interaction, fill = interaction)) +

  # ## Threshold
  # geom_hline(yintercept = cutoff[2], linewidth = 1.5,
  #            linetype = "dotted", col = "maroon") +

  geom_violin(width = 1) +
  geom_boxplot(width = 0.1, fill = "grey") +

  scale_fill_manual(values = c("coral")) +
  ggtitle("Predicted interactions") +
  xlab("Host status") + ylab("Log-odds predicted host status") +

  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))


## iv) Check new rate of true positives
# ## Prepare data
# predictions_loo1 <- data.frame(predictions_loo1,
#                           as.factor(ifelse(predictions_loo1$prediction.logodds.loo0 >
#                                             cutoff[2],
#                                  1, 0)))
# colnames(predictions_loo1) <- c(colnames(predictions_loo1)[1:5],
#                                 "prediction.logodds.loo0.cutoff")
## Create confidence matrix
conf_matrix <- table(ifelse(predictions_loo1$prediction.logodds.loo0 > cutoff[2], 1, 0),
                     predictions_loo1$interaction)

## Rate of true positives = 83.62%
round(conf_matrix[2]/(conf_matrix[1] + conf_matrix[2])*100, 2)
rm(conf_matrix, predictions1, interaction_data1)



## B) SECOND WAY - LEAVE ONE INTERACTION OUT AT A TIME
## i) Compute loo predictions
## I believe loo_predict() is the equivalent to predict(), and loo_linpred() to
## fitted(method = "linear")
## Each one takes <5 min to run
loo_pred_lodds <- loo_linpred(oak_mod043, type = "mean")
loo_pred_prob <- loo_predict(oak_mod043, type = "mean")


## ii) Compare against other method
predictions_loo1$prediction.logodds.loo <- loo_pred_lodds[which(predictions$host.status
                                                                == 1)]
predictions_loo1$thr <- "no"
for (i in 1:nrow(predictions_loo1)) {
  if (predictions_loo1$prediction.logodds.loo[i] >= cutoff[2] &
      predictions_loo1$prediction.logodds.loo0[i] < cutoff[2]) {
    predictions_loo1$thr[i] <- "loo"
  } else if (predictions_loo1$prediction.logodds.loo0[i] >= cutoff[2] &
        predictions_loo1$prediction.logodds.loo[i] < cutoff[2]) {
    predictions_loo1$thr[i] <- "loo0"
  }
}

# svg("results/loovsloo0.svg")
ggplot(predictions_loo1,
       aes(x = prediction.logodds.loo0, y = prediction.logodds.loo,
           col = thr)) +
  geom_abline(linetype = "dashed", col = "darkgrey", linewidth = 1.5) +
  geom_point() +
  scale_color_manual(values = c("black", "coral")) +
  # geom_hline(yintercept = cutoff[2], linewidth = 1.5,
  #            linetype = "dashed", col = "coral") +
  # geom_vline(xintercept = cutoff[2], linewidth = 1.5,
  #            linetype = "dashed", col = "coral") +
  xlab("'One 0 in turn' predictions") + ylab("LOO predictions") +
  ggtitle("Quercus sp. - Agrilus sp.") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")
# dev.off()

predictions_loo1$thr <- NULL

# plot(x = predictions_loo1$prediction.logodds.loo0,
#      y = predictions_loo1$prediction.logodds.loo)
# identify(x = predictions_loo1$prediction.logodds.loo0,
#          y = predictions_loo1$prediction.logodds.loo)
# which(predictions$host.status == 1)[35]
# interaction_data[6275, 1:3] ## Q. durata only hosts A. angelicus (6 hosts, 5 are oaks)


## iii) Compare against original predictions
predictions_loo <- data.frame(agrilus.sp = interaction_data$agrilus.sp,
                              quercus.sp = interaction_data$quercus.sp,
                              interaction = as.factor(interaction_data$interaction),
                              prediction.logodds.orig = predictions$prediction.logodds,
                              prediction.logodds.loo = loo_pred_lodds)

## Check the differences between the "loo" and the "original" predictions
plot(loo_pred_lodds ~ predictions$prediction.logodds)
abline(a = 0, b = 1, col = "steelblue")
# identify(y = loo_pred_lodds, x = predictions$prediction.logodds)

# Some of the points under the y = x line:
## They are all '1's, which makes sense: if they are taken out, it makes sense that the
## prediction would be slightly lower.
interaction_data[c(734,  753,  755, 1750, 2348, 4459, 4648, 4760, 5084, 5512, 5698,
                   6275, 6368, 6369, 6499, 6587, 6659, 6668, 6736, 6772, 6799, 6822,
                   6858, 7146, 7163, 7227), 1:3]
interaction_data[6369, 1:3]


## iv) Plot results
## Colour by interaction
ggplot(predictions_loo,
       aes(x = prediction.logodds.orig, y = prediction.logodds.loo, col = interaction)) +
  geom_abline(linetype = "dashed", col = "darkgrey") +
  geom_point() +
  xlab("Original model") + ylab("Leave One Out ('automatic')") +
  scale_color_manual(values = c("black", "coral")) +
  ggtitle("Quercus sp. - Agrilus sp.") +
  theme(plot.title = element_text(hjust = 0.5))

## Only plot 1s, colouring those points that have moved under the threshold
predictions_loo$thr <- "no"
for (i in 1:nrow(predictions_loo)) {
  if (predictions_loo$prediction.logodds.orig[i] >= cutoff[2] &
      predictions_loo$prediction.logodds.loo[i] < cutoff[2]) {
    predictions_loo$thr[i] <- "yes"
  }
}

# svg("results/loovsorigpreds.svg")
ggplot(predictions_loo[predictions_loo$interaction == 1, ],
       aes(x = prediction.logodds.orig, y = prediction.logodds.loo,
           col = thr)) +
  geom_abline(linetype = "dashed", col = "darkgrey", linewidth = 1.5) +
  geom_point() +
  scale_color_manual(values = c("black", "coral")) +
  # geom_hline(yintercept = cutoff[2], linewidth = 1.5,
  #            linetype = "dashed", col = "coral") +
  # geom_vline(xintercept = cutoff[2], linewidth = 1.5,
  #            linetype = "dashed", col = "coral") +
  xlab("Original model") + ylab("Leave One Out ('automatic')") +
  ggtitle("Quercus sp. - Agrilus sp.") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")
# dev.off()

predictions_loo[predictions_loo$interaction == 1 & predictions_loo$thr == "yes",]
predictions_loo$thr <- NULL

## Violin plot
ggplot(data = predictions_loo,
       aes(y = prediction.logodds.loo, x = interaction, fill = interaction)) +

  # ## Threshold
  # geom_hline(yintercept = cutoff[2], linewidth = 1.5,
  #            linetype = "dotted", col = "maroon") +

  geom_violin(width = 1) +
  geom_boxplot(width = 0.1, fill = "grey") +

  scale_fill_manual(values = c("aquamarine2", "coral")) +
  ggtitle("Predicted interactions") +
  xlab("Host status") + ylab("Log-odds predicted host status") +

  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))


## v) Check new rate of true positives
# ## Prepare data
# predictions_loo1 <- data.frame(predictions_loo1,
#                                as.factor(ifelse(predictions_loo1$prediction.logodds.loo >
#                                                   cutoff[2],
#                                                 1, 0)))
# colnames(predictions_loo1) <- c(colnames(predictions_loo1)[1:7],
#                                 "prediction.logodds.loo.cutoff")

## Create confidence matrix
conf_matrix <- table(ifelse(predictions_loo1$prediction.logodds.loo > cutoff[2], 1, 0),
                     predictions_loo1$interaction)

## Rate of true positives = 83.62%
round(conf_matrix[2]/(conf_matrix[1] + conf_matrix[2])*100, 2)

rm(conf_matrix, i)

# save.image(paste("/data/SBCS-NicholsLab/elvirahg/oak_analyses/01_original_models/",
#                  "/results/oak_models_final.RData", sep = ""))



#
##### XII)  MISCELANEOUS EXPLORATORY THINGS                                          #####
## A) RELATIONSHIP BETWEEN GEOGRAPHIC DISTANCE AND NUMBER OF INTERACTIONS
## i) Do oak species that are closer in space to Agrilus hosts tend to host more
## Agrilus species?
## Create DF
interaction_q <- data.frame(interaction_data[, c(1,3,5)] %>%
                 group_by(quercus.sp) %>%
                 summarise(mean.dist = mean(min.dist.mean.norm)))
interaction_q$interaction <- data.frame(interaction_data[, c(1,3,5)] %>%
                             group_by(quercus.sp) %>%
                             summarise(interaction = sum(interaction)))$interaction
str(interaction_q)

## Plot
boxplot(interaction_q$mean.dist ~ interaction_q$interaction,
        xlab = "No. Agrilus spp. hosted", ylab = "norm(mean(min.dist.geo))")
summary(lm(interaction_q$mean.dist ~ interaction_q$interaction))
abline(a = 0.866063, b = -0.050027)

t.test(interaction_q$mean.dist[interaction_q$interaction == 0],
       interaction_q$mean.dist[interaction_q$interaction > 0])
rm(interaction_q)

## ii) Do Agrilus species that are hosted by oaks that are closer in space to Agrilus
## hosts tend to attack more oak species?
## Create DF
interaction_a <- data.frame(interaction_data[, c(2,3,5)] %>%
                 group_by(agrilus.sp) %>%
                 summarise(mean.dist = mean(min.dist.mean.norm)))
interaction_a$interaction <- data.frame(interaction_data[, c(2,3,5)] %>%
                             group_by(agrilus.sp) %>%
                             summarise(interaction = sum(interaction)))$interaction
str(interaction_a)

## Plot
boxplot(interaction_a$mean.dist ~ interaction_a$interaction)
plot(interaction_a$mean.dist ~ interaction_a$interaction)
summary(lm(interaction_a$mean.dist ~ interaction_a$interaction))
abline(a = 0.883096, b = -0.011482)

rm(interaction_a)



## B) PHYLOGENETIC AND GEOGRAPHIC DATA FOR SOME AGRILUS SPECIES TO UNDERSTAND PREDICTIONS
## i) Load data
## Read nexus tree created with phylocom nodesig (first tree: host status info)
oak_nodesig <- read.nexus(paste0(pw, "input/quercus_nodesig_result.nex"))[[1]]

## Repace "Quercus" for "Q." (tip labels) and "_" for " "
oak_nodesig$tip.label <- gsub(pattern = "Quercus_", replacement = "Q. ",
                              oak_nodesig$tip.label)
oak_nodesig$tip.label <- gsub(pattern = "_", replacement = " ",
                              oak_nodesig$tip.label)

## Read file with info on number of Agrilus species hosted per Quercus species
hosts <- read.table(paste0(pw, "tmp/agrilus_quercus_hosts.txt"),
                    sep = "\t", header = F)
hosts <- data.frame(cbind(hosts$V2, hosts$V1))
colnames(hosts) <- c("quercus.sp", "no.agrilus")
hosts$quercus.sp <- as.character(hosts$quercus.sp)
hosts <- as.data.frame(table(hosts$quercus.sp))
colnames(hosts) <- c("quercus.sp", "no.agrilus")
hosts$quercus.sp <- gsub(pattern = "Quercus", replacement = "Q.", fixed = T,
                         x = hosts$quercus.sp)
str(hosts)


## ii) Plot data for A. auroguttatus
## Phylogenetic tree
## Group oaks according to host status: 0 = oaks not known to be hosts, 1 = native
## known hosts in CA, 2 = observed new hosts in introduced range (AZ), 3 = confirmed
## non-host in introduced range
hosts_aur <- list("1" = c("Q. conzattii", "Q. emoryi",
                          "Q. hypoleucoides", "Q. peduncularis"),
                  "2" = c("Q. agrifolia", "Q. kelloggii",
                        "Q. chrysolepis", "Q. engelmannii"),
                  "3" = "Q. arizonica")
oak_nodesig <-groupOTU(oak_nodesig, hosts_aur)

oak_tree_aur <- ggtree(oak_nodesig, col = "gray26", size = 1, layout = "circular") +
                geom_tiplab(aes(subset = group == 0, color = group), offset = 10,
                size = 1.7, fontface = "italic") +
                geom_tiplab(aes(subset = group != 0, color = group), offset = 10,
                            size = 1.7, fontface = "bold.italic") +
                scale_color_manual(values = c("black",
                                              brewer.pal(n = 3, name = "Dark2"))) +

                ggnewscale::new_scale_color() +
                geom_tippoint(aes(color = group)) +
                scale_color_manual(values = c("white",
                                   brewer.pal(n = 3, name = "Dark2"))) +

                ggnewscale::new_scale_fill() +
                geom_fruit(data = hosts, geom = geom_bar,
                           mapping = aes(y = quercus.sp, x =  no.agrilus, size = 100),
                           offset = 0.05, pwidth = 0.1,
                           stat = "identity", orientation = "y") +
                theme(legend.position = "none")

plot(oak_tree_aur)

pdists_aur <- subset(interaction_data, agrilus.sp == "A. auroguttatus" &
                       quercus.sp %in% c("Q. conzattii",
                                         "Q. emoryi",
                                         "Q. hypoleucoides",
                                         "Q. peduncularis",
                                         "Q. agrifolia",
                                         "Q. kelloggii",
                                         "Q. chrysolepis",
                                         "Q. engelmannii",
                                         "Q. arizonica"))[,c(1,3,12)]
pdists_aur <- pdists_aur[c(1, 3, 4:9, 2),]

# svg("results/A_auroguttatus_phylo_min_dist.svg")
barplot(pdists_aur$phylo.dist.min, names.arg = pdists_aur$quercus.sp,
        ylab = "min phylo dist", las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(1,1,1,1,2,2,2,2,3)])
# dev.off()

## Ranef values
## High values make it more likely for sth to be a host
ranef_aur <- data.frame(rbind(ranef(oak_mod043)$quercus.sp[,,1]['Q. conzattii',],     ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. emoryi',],        ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. hypoleucoides',], ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. peduncularis',],  ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. agrifolia',],     ## 2
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. kelloggii',],     ## 2
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. engelmannii',],   ## 2
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. chrysolepis',],   ## 2
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. arizonica',]      ## 3
                              ))
ranef_aur$quercus.sp <- c('Q. conzattii', 'Q. emoryi', 'Q. hypoleucoides',
                          'Q. peduncularis', 'Q. agrifolia', 'Q. kelloggii',
                          'Q. engelmannii', 'Q. chrysolepis', 'Q. arizonica')

# svg("results/A_auroguttatus_ranef.svg")
barplot(ranef_aur$Estimate, names.arg = ranef_aur$quercus.sp,
        ylab = "ranef estimate", las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(1,1,1,1,2,2,2,2,3)])
# dev.off()

## Geographic distances
gdists_aur <- subset(interaction_data, agrilus.sp == "A. auroguttatus" &
                                       quercus.sp %in% c("Q. conzattii",
                                                         "Q. emoryi",
                                                         "Q. hypoleucoides",
                                                         "Q. peduncularis",
                                                         "Q. agrifolia",
                                                         "Q. kelloggii",
                                                         "Q. chrysolepis",
                                                         "Q. engelmannii",
                                                         "Q. arizonica"))[,c(1,3,5)]
gdists_aur <- gdists_aur[c(1, 3, 4:9, 2),]

# svg("results/A_auroguttatus_geo_min_dist_mean_norm.svg")
barplot(gdists_aur$min.dist.mean.norm, names.arg = gdists_aur$quercus.sp,
        ylab = "min.dist.mean (geo)", las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(1,1,1,1,2,2,2,2,3)])
# dev.off()


## iii) Plot data for A. bilineatus
## Phylogenetic tree
## Group oaks according to host status: 0 = oaks not known to be hosts, 1 = native
## known hosts in Europe, 2 = observed new hosts in introduced range (AZ), 3 = confirmed
## non-host in introduced range
hosts_bil <- list("1" = c("Q. alba", "Q. coccinea", "Q. lyrata", "Q. macrocarpa",
                          "Q. michauxii", "Q. muehlenbergii", "Q. rubra",
                          "Q. velutina"),
                  "2" = c("Q. robur"))
oak_nodesig <-groupOTU(oak_nodesig, hosts_bil)

oak_tree_bil <- ggtree(oak_nodesig, col = "gray26", size = 1, layout = "circular") +
                geom_tiplab(aes(subset = group == 0, color = group), offset = 10,
                            size = 1.7, fontface = "italic") +
                geom_tiplab(aes(subset = group != 0, color = group), offset = 10,
                            size = 1.7, fontface = "bold.italic") +
                scale_color_manual(values = c("black",
                                              brewer.pal(n = 3, name = "Dark2"))) +

                ggnewscale::new_scale_color() +
                geom_tippoint(aes(color = group)) +
                scale_color_manual(values = c("white",
                                              brewer.pal(n = 3, name = "Dark2"))) +

                ggnewscale::new_scale_fill() +
                geom_fruit(data = hosts, geom = geom_bar,
                           mapping = aes(y = quercus.sp, x =  no.agrilus, size = 100),
                           offset = 0.05, pwidth = 0.1,
                           stat = "identity", orientation = "y") +
                theme(legend.position = "none")
plot(oak_tree_bil)

pdists_bil <- subset(interaction_data, agrilus.sp == "A. bilineatus" &
                       quercus.sp %in% c("Q. alba",
                                         "Q. coccinea",
                                         "Q. lyrata",
                                         "Q. macrocarpa",
                                         "Q. michauxii",
                                         "Q. muehlenbergii",
                                         "Q. rubra",
                                         "Q. velutina",
                                         "Q. robur"))[,c(1,3,12)]

# svg("A_bilineatus_phylo_min_dist.svg")
barplot(pdists_bil$phylo.dist.min, names.arg = pdists_bil$quercus.sp,
        ylab = "min phylo dist", las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(rep(1, 8), 2)])
# dev.off()

## Ranef values
## High values make it more likely for sth to be a host
ranef_bil <- data.frame(rbind(ranef(oak_mod043)$quercus.sp[,,1]['Q. alba',],          ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. coccinea',],      ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. lyrata',],        ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. macrocarpa',],    ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. michauxii',],     ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. muehlenbergii',], ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. rubra',],         ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. velutina',],      ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. robur',]          ## 2
                              ))
ranef_bil$quercus.sp <- c('Q. alba', 'Q. coccinea', 'Q. lyrata',
                          'Q. macrocarpa', 'Q. michauxii', 'Q. muehlenbergii',
                          'Q. rubra', 'Q. velutina', 'Q. robur')

# svg("A_bilineatus_raneff.svg")
barplot(ranef_bil$Estimate, names.arg = ranef_bil$quercus.sp,
        ylab = "ranef estimate", las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(1,1,1,1,1,1,1,1,2)])
# dev.off()

## Geographic distances
gdists_bil <- subset(interaction_data, agrilus.sp == "A. bilineatus" &
                                       quercus.sp %in% c("Q. alba",
                                                         "Q. coccinea",
                                                         "Q. lyrata",
                                                         "Q. macrocarpa",
                                                         "Q. michauxii",
                                                         "Q. muehlenbergii",
                                                         "Q. rubra",
                                                         "Q. velutina",
                                                         "Q. robur"))[,c(1,3,5)]

# svg("A_bilineatus_geo_min_dist_mean_norm.svg")
barplot(gdists_bil$min.dist.mean.norm, names.arg = gdists_bil$quercus.sp,
        ylab = "min.dist.mean (geo)", las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(rep(1, 8), 2)])
# dev.off()


## iv) Plot data for A. sulcicollis
## Phylogenetic tree
## Group oaks according to host status: 0 = oaks not known to be hosts, 1 = native
## known hosts in CA, 2 = observed new hosts in introduced range (AZ), 3 = confirmed
## non-host in introduced range
hosts_sul <- list("1" = c("Q. frainetto", "Q. petraea", "Q. pubescens",
                          "Q. robur", "Q. macrocarpa"),
                  "2" = c("Q. rubra"))
oak_nodesig <-groupOTU(oak_nodesig, hosts_sul)

oak_tree_sul <- ggtree(oak_nodesig, col = "gray26", size = 1, layout = "circular") +
                geom_tiplab(aes(subset = group == 0, color = group), offset = 10,
                            size = 1.7, fontface = "italic") +
                geom_tiplab(aes(subset = group != 0, color = group), offset = 10,
                            size = 1.7, fontface = "bold.italic") +
                scale_color_manual(values = c("black",
                                              brewer.pal(n = 3, name = "Dark2"))) +

                ggnewscale::new_scale_color() +
                geom_tippoint(aes(color = group)) +
                scale_color_manual(values = c("white",
                                              brewer.pal(n = 3, name = "Dark2"))) +

                ggnewscale::new_scale_fill() +
                geom_fruit(data = hosts, geom = geom_bar,
                           mapping = aes(y = quercus.sp, x =  no.agrilus, size = 100),
                           offset = 0.05, pwidth = 0.1,
                           stat = "identity", orientation = "y") +
                theme(legend.position = "none")

plot(oak_tree_sul)

pdists_sul <- subset(interaction_data, agrilus.sp == "A. sulcicollis" &
                       quercus.sp %in% c("Q. frainetto", "Q. petraea",
                                         "Q. pubescens", "Q. robur",
                                         "Q. macrocarpa",
                                         "Q. rubra"))[,c(1,3,12)]
pdists_sul <- pdists_sul[c(1:4,6,5),]

# svg("A_sulcicollis_phylo_min_dist.svg")
barplot(pdists_sul$phylo.dist.min, names.arg = pdists_sul$quercus.sp,
        ylab = "min phylo mean",
        las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(rep(1, 5), 2)])
# dev.off()

## Ranef values
## High values make it more likely for sth to be a host
ranef_sul <- data.frame(rbind(ranef(oak_mod043)$quercus.sp[,,1]['Q. frainetto',],   ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. petraea',],     ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. pubescens',],   ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. robur',],       ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. macrocarpa',],  ## 1
                              ranef(oak_mod043)$quercus.sp[,,1]['Q. rubra',]        ## 2
                              ))
ranef_sul$quercus.sp <- c('Q. frainetto', 'Q. petraea', 'Q. pubescens',
                          'Q. robur', 'Q. macrocarpa', 'Q. rubra')

# svg("A_sulcicollis_ranef.svg")
barplot(ranef_sul$Estimate, names.arg = ranef_sul$quercus.sp,
        ylab = "ranef estimate", las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(1,1,1,1,1,2)])
# dev.off()

## Geographic distances
gdists_sul <- subset(interaction_data, agrilus.sp == "A. sulcicollis" &
                                       quercus.sp %in% c("Q. frainetto", "Q. petraea",
                                                         "Q. pubescens", "Q. robur",
                                                         "Q. macrocarpa",
                                                         "Q. rubra"))[,c(1,3,5)]
gdists_sul <- gdists_sul[c(1:4,6,5),]

# svg("A_sulcicollis_geo_min_dist_mean_norm.svg")
barplot(gdists_sul$min.dist.mean.norm, names.arg = gdists_sul$quercus.sp,
        ylab = "min.dist.mean (geo)",
        las = 2, cex.names = 0.5,
        col = brewer.pal(n = 3, name = "Dark2")[c(rep(1, 5), 2)])
# dev.off()



## C) LOOK FOR JUMPS BETWEEN RED AND WHITE OAKS
## Reminder:
## A. auroguttatus: AZ (USA) to CA (USA), new hosts: Q. kelloggii Q. chrysolepis & Q. agrifolia
## A. bilineatus: introduced from N America to Turkey (no new reported hosts so far)
## A. sulcicollis: Europe to N America, novel hosts: Q. macrocarpa

str(agrilus_hosts)
subset(agrilus_hosts, agrilus.sp == "A. sulcicollis")
subset(agrilus_hosts, agrilus.sp == "A. auroguttatus")

## i) Divide opak spp. into red and white
## Red oaks = sec. Lobatae
## White oaks s.s. = sec. Quercus, Virentes & Ponticae

red_oaks <- c("Quercus crassifolia",       ## Quercus brachystachys is a syn
              "Quercus mcvaughii",
              "Quercus scytophylla",
              "Quercus hypoleucoides",
              "Quercus sideroxyla",
              "Quercus jonesii",
              "Quercus emoryi",
              "Quercus durifolia",
              "Quercus conzattii",
              "Quercus radiata",
              "Quercus calophylla",
              "Quercus fulva",
              "Quercus urbani",
              "Quercus viminea",
              "Quercus eduardi",
              "Quercus crassipes",
              "Quercus confertifolia",     ## Quercus gentryi is a syn
              "Quercus castanea",
              "Quercus acutifolia",
              "Quercus grahamii",
              "Quercus humboldtii",
              "Quercus costaricensis",
              "Quercus seemannii",         ## Quercus eugeniifolia is a syn
              "Quercus cortesii",
              "Quercus benthamii",         ## Quercus lowilliamsii is a syn
              "Quercus delgadoana",
              "Quercus sapotifolia",
              "Quercus iltisii",
              "Quercus crispifolia",
              "Quercus elliptica",
              "Quercus uxoris",
              "Quercus aristata",
              "Quercus planipocula",
              "Quercus laurina",
              "Quercus pinnativenulosa",
              "Quercus sartorii",
              "Quercus mexicana",
              "Quercus affinis",
              "Quercus gravesii",
              "Quercus canbyi",
              "Quercus incana",
              "Quercus hemisphaerica",
              "Quercus inopina",
              "Quercus myrtifolia",
              "Quercus laevis",
              "Quercus nigra",
              "Quercus arkansana",
              "Quercus pumila",            ## Quercus elliottii is a syn
              "Quercus phellos",
              "Quercus laurifolia",
              "Quercus falcata",
              "Quercus marilandica",
              "Quercus pagoda",
              "Quercus imbricaria",
              "Quercus georgiana",
              "Quercus ilicifolia",
              "Quercus acerifolia",
              "Quercus shumardii",
              "Quercus buckleyi",
              "Quercus rubra",
              "Quercus ellipsoidalis",
              "Quercus velutina",
              "Quercus coccinea",
              "Quercus texana",
              "Quercus palustris",
              "Quercus wislizeni",
              "Quercus parvula",
              "Quercus agrifolia",
              "Quercus kelloggii")

red_oaks[!(red_oaks %in% gsub("Q.", "Quercus", interaction_data$quercus.sp))]

white_oaks <- c("Quercus turbinella",
                "Quercus ajoensis",
                "Quercus toumeyi",
                "Quercus striatula",
                "Quercus grisea",
                "Quercus arizonica",
                "Quercus engelmannii",
                "Quercus oblongifolia",
                "Quercus rugosa",
                "Quercus greggii",
                "Quercus diversifolia",
                "Quercus obtusata",
                "Quercus potosina",
                "Quercus peduncularis",
                "Quercus laeta",
                "Quercus chihuahuensis",
                "Quercus deserticola",
                "Quercus glaucoides",
                "Quercus resinosa",
                "Quercus magnoliifolia",         ## Quercus nudinervis is a syn
                "Quercus subspathulata",
                "Quercus liebmannii",
                "Quercus segoviensis",
                "Quercus purulhana",
                "Quercus glaucescens",
                "Quercus lancifolia",
                "Quercus copeyensis",
                "Quercus insignis",
                "Quercus corrugata",
                "Quercus glabrescens",
                "Quercus martinezii",
                "Quercus germana",
                "Quercus vaseyana",
                "Quercus pungens",
                "Quercus hinckleyi",
                "Quercus polymorpha",
                "Quercus mohriana",
                "Quercus laceyi",
                "Quercus margarettae",
                "Quercus austrina",
                "Quercus havardii",
                "Quercus chapmanii",
                "Quercus similis",
                "Quercus stellata",
                "Quercus boyntonii",
                "Quercus oglethorpensis",
                "Quercus sinuata",
                "Quercus infectoria",            ## Quercus boissieri is a syn
                "Quercus kotschyana",
                "Quercus pubescens",
                "Quercus vulcanica",
                "Quercus faginea",
                "Quercus macranthera",
                "Quercus frainetto",
                "Quercus dalechampii",
                "Quercus lusitanica",
                "Quercus petraea",              ## Quercus cedrorum is a syn
                "Quercus pyrenaica",
                "Quercus canariensis",
                "Quercus robur",
                "Quercus hartwissiana",
                "Quercus aliena",
                "Quercus griffithii",
                "Quercus fabrei",
                "Quercus serrata",
                "Quercus mongolica",
                "Quercus dentata",              ## Quercus yunnanensis is a syn
                "Quercus alba",
                "Quercus michauxii",
                "Quercus montana",
                "Quercus bicolor",
                "Quercus lyrata",
                "Quercus macrocarpa",
                "Quercus muehlenbergii",
                "Quercus prinoides",
                "Quercus pacifica",
                "Quercus dumosa",
                "Quercus douglasii",
                "Quercus johntuckeri",
                "Quercus corneliusmulleri",
                "Quercus berberidifolia",
                "Quercus durata",
                "Quercus garryana",
                "Quercus lobata",
                "Quercus geminata",
                "Quercus virginiana",
                "Quercus minima",
                "Quercus oleoides",
                # "Quercus sagraeana",          ## not in my DS (not in GBIF)
                "Quercus brandegeei",
                "Quercus fusiformis",
                "Quercus pontica",
                "Quercus sadleriana")

white_oaks[!(white_oaks %in% gsub("Q.", "Quercus", interaction_data$quercus.sp))]
which(white_oaks %in% red_oaks)

## Divide oaks into red and white
oaks <- data.frame(matrix(nrow = 236, ncol = 2))
colnames(oaks) <- c("quercus.sp", "type")
oaks$quercus.sp <- gsub("Q.", "Quercus",
                        unique(grep("Q.", interaction_data$quercus.sp, value = TRUE)))
str(oaks)

oaks[oaks$quercus.sp %in% red_oaks, ]$type <- "red"
oaks[oaks$quercus.sp %in% white_oaks, ]$type <- "white"

## 236 spp in total, out of which 75 are not white or red
length(which(!(is.na(oaks$type))))
length(which(is.na(oaks$type)))


## ii) Find which Agrilus spp. use red/white oaks
agrilus_redwhite <- data.frame(matrix(ncol = 4,
                                      nrow = length(unique(interaction_data$agrilus.sp))))
colnames(agrilus_redwhite) <- c("agrilus.sp", "red.hosts", "white.hosts", "na.hosts")
agrilus_redwhite$agrilus.sp <- unique(interaction_data$agrilus.sp)
agrilus_redwhite$colour <- NA

for (i in 1:nrow(agrilus_redwhite)) {
  agrilus <- agrilus_redwhite$agrilus.sp[i]
  hosts <- gsub("Q.", "Quercus",
                subset(interaction_data, agrilus.sp == agrilus & interaction == 1)$quercus.sp)
  colours <- subset(oaks, quercus.sp %in% hosts)$type
  colours[is.na(colours)] <- 0
  colours <- round(table(colours)/sum(table(colours)), 2)
  agrilus_redwhite$red.hosts[i] <- colours[2]
  agrilus_redwhite$white.hosts[i] <- colours[3]
  agrilus_redwhite$na.hosts[i] <- colours[1]

  cols <- which(agrilus_redwhite[i, 2:4] > 0)
  if (length(cols > 0) == 1) {
    if (cols == 1) { agrilus_redwhite[i, ]$colour <- "red" }
    if (cols == 2) { agrilus_redwhite[i, ]$colour <- "white" }
    if (cols == 3) { agrilus_redwhite[i, ]$colour <- "na" }
  } else if (length(cols) == 2) {
    if (setequal(cols, c(1, 2))) { agrilus_redwhite[i, ]$colour <- "pink" }
    if (setequal(cols, c(1, 3))) { agrilus_redwhite[i, ]$colour <- "red.na" }
    if (setequal(cols, c(2, 3))) { agrilus_redwhite[i, ]$colour <- "white.na" }
  } else if (length(cols) == 3) {
    agrilus_redwhite[i, ]$colour <- "pink.na"
  }
}; rm(i, agrilus, hosts, colours, cols)
agrilus_redwhite[is.na(agrilus_redwhite)] <- 0

str(agrilus_redwhite)
unique(agrilus_redwhite$colour)


## iii) Find whether there are jumps in the predictions
agrilus_redwhite$red.hosts.pred <- NA
agrilus_redwhite$white.hosts.pred <- NA
agrilus_redwhite$na.hosts.pred <- NA
agrilus_redwhite$colour.pred <- NA

for (i in 1:nrow(agrilus_redwhite)) {
  agrilus <- agrilus_redwhite$agrilus.sp[i]
  hosts <- gsub("Q.", "Quercus",
                subset(predictions, agrilus.sp == agrilus & prediction.cutoff == 1)$quercus.sp)
  known <- gsub("Q.", "Quercus",
                subset(interaction_data, agrilus.sp == agrilus & interaction == 1)$quercus.sp)
  hosts <- hosts[!(hosts %in% known)]
  colours <- subset(oaks, quercus.sp %in% hosts)$type
  colours[is.na(colours)] <- 0
  colours <- round(table(colours)/sum(table(colours)), 2)
  agrilus_redwhite$red.hosts.pred[i] <- colours[2]
  agrilus_redwhite$white.hosts.pred[i] <- colours[3]
  agrilus_redwhite$na.hosts.pred[i] <- colours[1]

  cols <- which(agrilus_redwhite[i, 6:8] > 0)
  if (length(cols > 0) == 1) {
    if (cols == 1) { agrilus_redwhite[i, ]$colour.pred <- "red" }
    if (cols == 2) { agrilus_redwhite[i, ]$colour.pred <- "white" }
    if (cols == 3) { agrilus_redwhite[i, ]$colour.pred <- "na" }
  } else if (length(cols) == 2) {
    if (setequal(cols, c(1, 2))) { agrilus_redwhite[i, ]$colour.pred <- "pink" }
    if (setequal(cols, c(1, 3))) { agrilus_redwhite[i, ]$colour.pred <- "red.na" }
    if (setequal(cols, c(2, 3))) { agrilus_redwhite[i, ]$colour.pred <- "white.na" }
  } else if (length(cols) == 3) {
    agrilus_redwhite[i, ]$colour.pred <- "pink.na"
  }
}; rm(i, agrilus, hosts, known, colours, cols)
agrilus_redwhite[is.na(agrilus_redwhite)] <- 0

str(agrilus_redwhite)
unique(agrilus_redwhite$colour)



## D) FIND WHETHER THERE IS A CORRELATION BETWEEN NO. HOSTS & PHYLO SIGNAL BTW. HOSTS
agrilus_no <-  as.data.frame(table(agrilus_hosts[grep("Quercus",
                                                      agrilus_hosts$plant.sp),]$agrilus.sp))
colnames(agrilus_no) <- c("agrilus.sp", "no.hosts")
agrilus_no <- agrilus_no[agrilus_no$no.hosts > 1, ]
head(agrilus_no)

## Phylo cov matrix
oak_phylo <- paste0("/data/SBCS-NicholsLab/elvirahg/oak_analyses/",
                    "01_original_models/input/quercus_nodesig_result.nex")
oak_phylo <- read.nexus(oak_phylo)[[1]]

## Edit node label information
oak_phylo$node.label <- gsub("'", "", oak_phylo$node.label)
oak_phylo$node.label <- gsub("N", "", oak_phylo$node.label)
oak_phylo$tip.label <- gsub(pattern = "_", replacement = " ", oak_phylo$tip.label)
oak_phylo$tip.label <- gsub(pattern = "Quercus", replacement = "Q.",
                            oak_phylo$tip.label, fixed = TRUE)

## Create phylo matrix
oak_phylo_cov <- vcv.phylo(oak_phylo)
rm(oak_phylo)

## For each Agrilus species, extract the name of its host species and then look for the mean
## phylo distance between them
agrilus_no$mean.phylo.dist.hosts <- NA

for (i in 1:nrow(agrilus_no)) {
  agrilus <- agrilus_no$agrilus.sp[i]
  hosts <- subset(interaction_data, agrilus.sp == agrilus & interaction == 1)$quercus.sp

  # if (length(hosts) == 1) {
  #   mean_phylo_dist <- max(oak_phylo_cov)
  # }
  # else {
  nos <- which(colnames(oak_phylo_cov) %in% hosts)
  mean_phylo_dist <- c()

  for (j in 1:length(hosts)) {

    dist_j <- mean(oak_phylo_cov[nos[j], nos[-j]])
    mean_phylo_dist <- c(mean_phylo_dist, dist_j)

  }
  mean_phylo_dist <- mean(mean_phylo_dist)
  # }
  agrilus_no$mean.phylo.dist.hosts[i] <- mean_phylo_dist
}; rm(dist_j, mean_phylo_dist, agrilus, hosts, i, j, nos)

head(agrilus_no)

mod1 <- lm(no.hosts ~ mean.phylo.dist.hosts, data = agrilus_no)
summary(mod1)

plot(no.hosts ~ mean.phylo.dist.hosts, data = agrilus_no,
     pch = 16)
abline(a = 5.76034, b = -0.07327)

mod2 <- glm(no.hosts ~ mean.phylo.dist.hosts,
            family = poisson(link = log),
            data = agrilus_no)
summary(mod2)

## https://thestatsgeek.com/2014/04/26/deviance-goodness-of-fit-test-for-poisson-regression/
# pchisq(mod2$deviance, df = mod2$df.residual,
#        lower.tail = FALSE)                      ## H0: model is correctly specified

ggplot(data =  agrilus_no, aes(x = mean.phylo.dist.hosts, y = no.hosts)) +
  geom_point() +
  geom_smooth()  # method = "glm", method.args = list(family = "poisson")




#
#####



#
