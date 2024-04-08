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
#####       COMPUTE FRITZ AND PURVIS' D (PHYLOGENETIC SIGNAL)                        #####
# A) READ IN AND PREPARE PHYLOGENETIC TREE OF OAKS
# Source: doi.org/10.1111/nph.16162, github.com/andrew-hipp/global-oaks-2019

# Read tree (NB, with tmp/tr.*Stem_accepted*.tre results are very similar)
# NB, synonmyms have been substituted to align with WCVP accepted names
oak_tree <- read.tree("input/tr.singletons.correlated.1.taxaGrepCrown_accepted_names.tre")

# Remove duplicated tips that form a sister clade
# oak_tree$tip.label[duplicated(oak_tree$tip.label)]
oak_tree <- drop.tip(oak_tree,
                     c("Quercus_magnoliifolia|Mexico|Durango|2017124|QUE002013",
                       "Quercus_crassifolia|Mexico|Durango|2017108|QUE001997",
                       "Quercus_agrifolia|USA|CA|CA-DAV-MH55|QUE000328",
                       "Quercus_coccifera|Italy|Apulia|TUS13-013|QUE001988"))

# Rename tips (only keep species info)
for (i in 1:length(oak_tree$tip.label)) {
  # Extract tip label
  tip <- oak_tree$tip.label[i]
  # Only keep species name info, and change "Quercus" for "Q."
  oak_spp <- gsub(oak_tree$tip.label[i], pattern = "Quercus_([^|]*).*",
                  replacement = "Q._\\1")
  # Substitute any ".", "._" or "_" for ". ", unless "." is at the end of the string
  oak_spp <- gsub("_", " ", oak_spp)
  # Replace old label with new one
  oak_tree$tip.label[i] <- oak_spp
}
# Relabel Q. litoralis (accepted name acc. to WCVP v3 is Atuna excelsa subsp. excelsa)
oak_tree$tip.label[220] <- "Q. litoralis (Atuna excelsa)"

# Remove one tip if name is duplicated (even if tips from the same sp. are not sister
# clades). Results are very similar when removing the other tip
oak_tree$tip.label[duplicated(oak_tree$tip.label)]
grep("Q. arizonica", oak_tree$tip.label)   # 6, 14
oak_tree$tip.label[6] <- "Q. arizonica2"

grep("Q. laeta", oak_tree$tip.label)       # 16, 18
oak_tree$tip.label[16] <- "Q. laeta2"

grep("Q. conzattii", oak_tree$tip.label)   # 114, 116
oak_tree$tip.label[114] <- "Q. conzattii2"

grep("Q. benthamii", oak_tree$tip.label)   # 131, 132
oak_tree$tip.label[131] <- "Q. benthamii2"

grep("Q. rehderiana", oak_tree$tip.label)  # 179, 182
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

# Plot tree and make sure it looks OK
# plot(oak_tree)



# B) READ IN HOST DATA
# Read in host data for oak trees
oak_hosts <- read.table("input/no_agrilus_quercus_hosts.txt", sep = "\t")
colnames(oak_hosts) <- c("quercus.sp", "no.agrilus.sp")
oak_hosts <- as.data.frame(lapply(oak_hosts, gsub, pattern = "Quercus",
                                  replacement = "Q.", fixed = TRUE))
str(oak_hosts)

# Create a 0/1 host status DF
oak_host_status <- as.data.frame(matrix(nrow = length(oak_tree$tip.label), ncol = 2))
colnames(oak_host_status) <- c("quercus.sp", "host.status")
oak_host_status$quercus.sp <- oak_tree$tip.label
oak_host_status$host.status <- 0

for (i in 1:nrow(oak_host_status)) {
  if (oak_host_status$quercus.sp[i] %in% oak_hosts$quercus.sp) {
    oak_host_status$host.status[i] <- 1
  }
}



# C) COMPUTE SIGNAL
# Create a comparative data object to match the spp. names in the phylo tree to the DS
oak_tree$node.label <- NULL
oak_comparative <- comparative.data(phy = oak_tree,
                                    data = oak_host_status,
                                    # Name of the column with spp. names
                                    names.col = quercus.sp,
                                    # Store a variance covariance matrix of your tree
                                    vcv = TRUE,
                                    # Don't remove spp. w/o data for certain variables
                                    na.omit = FALSE,
                                    # Shows dropped spp.
                                    warn.dropped = TRUE)

# Now compute strength of the phylogenetic signal (Fritz & Purvis' D metric)
d_oaks <- phylo.d(data = oak_comparative,
                  binvar = host.status, permut = 1000)
d_oaks

# D = 0.55496, which suggests clumping is moderately strong (according to F&P 2010,
# 0.56 is "moderately strong").
# Phylogenetic pattern differs significantly from the Brownian expectation (prob. of
# estimating D under Brownian evolution = 0.001)
# Also differs significantly from 1 (prob of estimating D with random phylogenetic
# structure = 0).

# See this graphically
plot(d_oaks, lwd = 5, cex.axis = 2.5, cex.lab = 2)
axis(side = 2, lwd = 3, lab = F)
axis(side = 1, lwd = 3, lab = F)
box(lwd=3)



#
#####       PLOT HOT CLADES                                                          #####
# A) EXTRACT NO. ENTRIES IN GBIF FOR OAK SPP. IN PHYLOGENY
# i) Use no. of entries according to GBIF
# Read the phylogenetic tree of oaks and extract spp. names (also from Hipp et al., 2019)
oaks <- read.tree("input/quercus_hipp19_singleton_crown_sp_level.tre")
oaks <- gsub(oaks$tip.label, pattern = "Q.", replacement = "Quercus", fixed = T)
oaks <- gsub(oaks, pattern = "_", replacement = " ", fixed = T)
oaks <- gsub(oaks, pattern = " -Atuna excelsa", replacement = "", fixed = T)
oaks <- gsub(oaks, pattern = "?? ", replacement = "", fixed = T)

# Extract their GBIF info (takes ca. 1 min to run).
oak_keys_info <- lapply(oaks, name_lookup, rank = "SPECIES",
                        higherTaxonKey = "220", status = "Accepted", limit = 1)

# Create a dataframe to store oak species name and their GBIF key. If there is no key
# info for the current sp., save its index in manual_keys. Takes < 1 min to run.
oak_keys <- as.data.frame(matrix(ncol = 2, nrow = length(oak_keys_info)))
colnames(oak_keys) <- c("oak", "key")
manual_keys <- c()
for (i in 1:length(oak_keys_info)) {
  if (!is.null(oak_keys_info[[i]][4]$hierarchies)) {
    sp <- oak_keys_info[[i]][2]$data$species
    # So that it doesn't include ""Quercus faginea-mirbeckii" (sp. in "oaks" is
    # "Quercus faginea")
    if (sp %in% oaks) {
      oak_keys[i,1] <- sp
      oak_keys[i,2] <- oak_keys_info[[i]][2]$data$key
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
oaks[manual_keys]
no <- grep(oaks[manual_keys][1], x = oaks)    # "Quercus new species"
oak_keys[no,1] <- oaks[manual_keys][1]        # "Quercus new species"
oak_keys[no,2] <- NA                          # "Quercus new species"
no <- grep(oaks[manual_keys][2], x = oaks)    # "Quercus corrugata"
oak_keys[no,1] <- oaks[manual_keys][2]        # "Quercus corrugata"
oak_keys[no,2] <- 2880406                     # "Quercus corrugata"*
no <- grep(oaks[manual_keys][3], x = oaks)    # "Quercus faginea"
oak_keys[no,1] <- oaks[manual_keys][3]        # "Quercus faginea"
oak_keys[no,2] <- 2881480                     # "Quercus faginea"***
no <- grep(oaks[manual_keys][4], x = oaks)    # "Quercus calophylla"
oak_keys[no,1] <- oaks[manual_keys][4]        # "Quercus calophylla"
oak_keys[no,2] <- 2880691                     # "Quercus calophylla"*
no <- grep(oaks[manual_keys][5], x = oaks)    # "Quercus confertifolia"
oak_keys[no,1] <- oaks[manual_keys][5]        # "Quercus confertifolia"
oak_keys[no,2] <- 2878551                     # "Quercus confertifolia"*
no <- grep(oaks[manual_keys][6], x = oaks)    # "Quercus sartorii"
oak_keys[no,1] <- oaks[manual_keys][6]        # "Quercus sartorii"
oak_keys[no,2] <- 2880720                     # "Quercus sartorii"*
no <- grep(oaks[manual_keys][7], x = oaks)    # "Quercus palustris"
oak_keys[no,1] <- oaks[manual_keys][7]        # "Quercus palustris"
oak_keys[no,2] <- 8313153                     # "Quercus palustris"
no <- grep(oaks[manual_keys][8], x = oaks)    # "Quercus sp. nov."
oak_keys[no,1] <- oaks[manual_keys][8]        # "Quercus sp. nov."
oak_keys[no,2] <- NA                          # "Quercus sp. nov."
no <- grep(oaks[manual_keys][9], x = oaks)    # "Quercus litoralis"
oak_keys[no,1] <- oaks[manual_keys][9]        # "Quercus litoralis"
oak_keys[no,2] <- 10839442                    # "Quercus litoralis"**
no <- grep(oaks[manual_keys][10], x = oaks)   # "Quercus acuta"
oak_keys[no,1] <- oaks[manual_keys][10]       # "Quercus acuta"
oak_keys[no,2] <- 2876863                     # "Quercus acuta"***

nrow(oak_keys) == length(oaks)

entries <- c()
for (i in 1:nrow(oak_keys)) {
  entries <- c(entries, occ_count(taxonKey = oak_keys$key[i],
                                  georeferenced = TRUE, from = 1950, type = "count"))
}

entries <- data.frame(cbind(oak_keys$oak, entries))
colnames(entries) <- c("oak.sp", "no.entries")

entries$no.entries <- as.numeric(entries$no.entries)
entries$oak.sp <- gsub("Quercus", "Q.", entries$oak.sp, fixed = T)
entries$oak.sp <- gsub(pattern = "-", replacement = "",
                       entries$oak.sp, fixed = T)

entries[entries$oak.sp == "Q. sp.nov.", ]$no.entries <- 0
entries[entries$oak.sp == "Q. new species", ]$no.entries <- 0


# B) LOAD  REST OF DATA
# Read nexus tree created with phylocom nodesig (first tree: host status info)
oak_nodesig <- read.nexus("results/quercus_nodesig_result.nex")[[1]]

# Repace "Quercus" for "Q." (tip labels) and "_" for " "
oak_nodesig$tip.label <- gsub(pattern = "Quercus_", replacement = "Q. ",
                              oak_nodesig$tip.label)
oak_nodesig$tip.label <- gsub(pattern = "_", replacement = " ",
                              oak_nodesig$tip.label)

# Read file with info on number of Agrilus species hosted per Quercus species
hosts <- read.table("input/agrilus_quercus_hosts.txt", sep = "\t", header = F)
hosts <- data.frame(cbind(hosts$V2, hosts$V1))
colnames(hosts) <- c("quercus.sp", "no.agrilus")

hosts$quercus.sp <- as.character(hosts$quercus.sp)
hosts <- as.data.frame(table(hosts$quercus.sp))
colnames(hosts) <- c("quercus.sp", "no.agrilus")
str(hosts)

# Replace "Quercus" with "Q."
hosts$quercus.sp <- gsub(pattern = "Quercus", replacement = "Q.", fixed = T,
                         x = hosts$quercus.sp)



# C) PLOT TREE WITH HOT CLADES
# Check with Mesquite (in Windows) which clades are hot/cold, and colour them
# accordingly with R
oak_nodesig$node.label <- gsub(pattern = "'", replacement = "",
                               oak_nodesig$node.label)

oak_nodesig <- groupClade(.data = oak_nodesig, .node = c(
                          "N246", "N247",
                          "N293", "N294", "N295", "N296", "N297", "N298",
                          "N316", "N317",
                          "N342", "N343", "N344",
                          "N381",
                          "N404",
                          "N419",
                          "N431",
                          "N443", "N444", "N445", "N446"))

cols <- c("grey10",
          brewer.pal(n = 9, name = "Blues")[4:5],
          brewer.pal(n = 9, name = "YlOrRd")[4:9],
          brewer.pal(n = 9, name = "YlOrRd")[5:6],
          brewer.pal(n = 9, name = "Blues")[4:6],
          brewer.pal(n = 9, name = "YlOrRd")[4],
          brewer.pal(n = 9, name = "YlOrRd")[4],
          brewer.pal(n = 9, name = "YlOrRd")[4],
          brewer.pal(n = 9, name = "YlOrRd")[4],
          brewer.pal(n = 9, name = "Blues")[4:7])

oak_tree <- ggtree(oak_nodesig,
                   aes(color = group), size = 1.3, layout = "circular") +
            scale_color_manual(values = cols) +
            geom_tiplab(color = "black", offset = 15, size = 1.7, fontface = "italic") +
            # theme(legend.position = "none") +
            geom_point2(aes(subset = (label == "N246")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N247")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N293")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N294")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N295")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N296")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N297")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N298")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N316")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N317")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N342")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N343")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N344")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N381")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N404")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N419")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N431")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N443")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N444")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N445")), shape = 20, size = 5) +
            geom_point2(aes(subset = (label == "N446")), shape = 20, size = 5)
# plot(oak_tree)

# NB, re stat = "identity": For geom_bar(), the default (stat = "count") is to count
# the rows for each x value. With stat = "identity", you're telling ggplot2 to skip
# the aggregation and that you'll provide the y values. From: stackoverflow.com/
# questions/59008974/why-is-stat-identity-necessary-in-geom-bar-in-ggplot
oak_tree2 <- oak_tree + ggnewscale::new_scale_fill() +
             geom_fruit(data = hosts,
                        geom = geom_bar,
                        mapping = aes(y = quercus.sp, x =  no.agrilus, size = 100),
                        offset = 0.15,
                        pwidth = 0.1,
                        stat = "identity",                   # see above
                        orientation = "y"                    # axis orientation
                        # axis.params = list(axis = "x",     # add axis text
                        #               text.angle = -45,    # text size of axis
                        #               hjust = 0            # adjust text horiz. pos.
                        #               )
                        # grid.params = list()               # add grid line
                        )

plot(oak_tree2)



# D) ADD INFO ON NO. GBIF ENTRIES
# oak_tree2$data$label[1:238][!(oak_tree2$data$label[1:238] %in% entries$oak.sp)]
oak_tree3  <- oak_tree + ggnewscale::new_scale_fill()
# entries[entries$no.entries >= 100000,]$no.entries <- 100000

oak_tree3 <- oak_tree2 %<+% entries +
  ggnewscale::new_scale_fill() +
  geom_tippoint(aes(x = x+4, size = no.entries), # , fill = no.entries
                 shape = 21, fill = alpha("#86929C", .5), colour = "grey") +
  scale_size_continuous(range = c(0.1, 10),
                        limits = c(1, 500000),
                        breaks = c(1, 1000, 10000,
                                   100000, 500000)) +
  guides(color = "none") +
  theme(legend.position = c(0.9, 0.95))
  # scale_fill_gradient(low = "lightgoldenrod1",
  #                     high = "red3",
  #                     na.value = NA,
  #                     limits = c(0.1, 8))

plot(oak_tree3)



#