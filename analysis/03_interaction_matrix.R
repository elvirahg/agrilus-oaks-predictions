##########################################################################################
####                           AGRILUS OAK HOSTS ANALYSES                          ######
####                           ELVIRA HERNANDEZ GUTIERREZ                          ######
##########################################################################################

####       SET ENVIRONMENT                                                          #####
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
####       CREATE A DF WITH INFO ON OAK - INTERACTION DATA                          #####
# A) READ PLANT GEO DATA FROM GBIF
# Read in file
geo_info <- read.table("/tmp/gbif_plants_clean.tsv", header = T, sep = "\t")
geo_info <- geo_info[, c("plant.sp", "lon", "lat")]
str(geo_info)

# Check that all plant species that should be there are present
plants <- read.tree("/input/quercus_hipp19_singleton_crown_sp_level.tre")

plants <- gsub(plants$tip.label, pattern = "Q.", replacement = "Quercus",
                   fixed = T)
plants <- gsub(plants, pattern = "_", replacement = " ",
                   fixed = T)
plants <- gsub(plants, pattern = " -Atuna excelsa", replacement = "",
                   fixed = T)

plants <- gsub(plants, pattern = "?? ", replacement = "", fixed = T)
plants <- c(plants,
                unique(read.table("input/non_oak_hosts.tsv",
                                  header = T, sep = "\t")$hosts))
length(plants) == (238 + 19)  # 238 Quercus spp. + 19 extra hosts



# B) EXTRACT INFO ON AGRILUS SPP. THAT EXPLOIT OAKS AND THEIR HOSTS
agrilus_hosts <- read.table("input/agrilus_quercus_hosts.txt",
                            header = F, sep = "\t")
colnames(agrilus_hosts) <- c("agrilus", "hosts")

agrilus_hosts <- rbind(agrilus_hosts,
                       read.table("input/non_oak_hosts.tsv",
                                  header = T, sep = "\t"))
agrilus_hosts <- agrilus_hosts[,c(2,1)]
colnames(agrilus_hosts) <- c("plant.sp", "agrilus.sp")
agrilus_hosts$agrilus.sp <- gsub("Agrilus", "A.", agrilus_hosts$agrilus.sp)

# Remove any hosts that are not in the phylogeny (only Q. gambelii) from agrilus_hosts
agrilus_hosts <- agrilus_hosts[!grepl("gambelii", agrilus_hosts$plant.sp),]


# IMPORTANT: I have also decided to remove any reported hosts that are non-native
# for those Agrilus species that have been introduced into new areas
# A. auroguttatus: introduced from AZ (USA) to CA (USA). In CA, reported novel larval
# hosts are Q. kelloggii, Q. chrysolepis & Q. agrifolia
# A. bilineatus: introduced from N America to Turkey (no new reported hosts so far)
# A. sulcicollis: introduced from Europe to N America, where Q. macrocarpa (in my DS)
# is a new reported larval host
head(agrilus_hosts[agrilus_hosts$plant.sp == "Quercus macrocarpa" &
                   agrilus_hosts$agrilus.sp == "A. sulcicollis",])

agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Quercus agrifolia"
                               & agrilus_hosts$agrilus.sp == "A. auroguttatus")), ]
agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Quercus macrocarpa"
                               & agrilus_hosts$agrilus.sp == "A. sulcicollis")), ]

# IMPORTANT: Let's also remove any larval hosts according to Jendek & Polakova (2014)
# (i.e., there's at least one reference with '!') but with a  confidence index < 3.
# A. graminis: Corylus avellana, Ostrya carpinifolia, Euonymus europaeus and Castanea
# sativa.
# A. obscuricollis: Ficus carica
# A. relegatus (reported as A. dualis)	Quercus robur
agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Ficus carica"
                                 & agrilus_hosts$agrilus.sp == "A. obscuricollis")), ]

agrilus_hosts <- agrilus_hosts[-(which(agrilus_hosts$plant.sp == "Quercus robur"
                                 & agrilus_hosts$agrilus.sp == "A. relegatus")), ]

agrilus_hosts <- agrilus_hosts[-(which((agrilus_hosts$plant.sp == "Corylus avellana" |
                                 agrilus_hosts$plant.sp == "Ostrya carpinifolia" |
                                 agrilus_hosts$plant.sp == "Euonymus europaeus" |
                                 agrilus_hosts$plant.sp == "Castanea sativa") &
                                 agrilus_hosts$agrilus.sp == "A. graminis")), ]



# C) CREATE INTERACTION DF
# Extract all plant and Agrilus spp.
oaks <- grep("Quercus", unique(geo_info$plant.sp), value = TRUE)
agrilus <- unique(agrilus_hosts$agrilus.sp)

# Check which oak spp from the phylogeny & which non-oak hosts are missing in "oaks"
# (they won't be present in the final dataset)
setdiff(plants, oaks)
rm(plants)

# Create a DF with all Agrilus sp. - oak sp. combinations (interaction_data)
interaction_data <- expand.grid(oaks, agrilus)
colnames(interaction_data) <- c("plant.sp", "agrilus.sp")
nrow(interaction_data)

# Append info on Agrilus hosts (i.e., oak spp. and the beetle sp. they host) at the
# end of interaction_data
interaction_data <- rbind(interaction_data,
                          agrilus_hosts[grep("Quercus", agrilus_hosts$plant.sp),])

# Make it so that hosts = 1, non-hosts = 0 (i.e., count number of entries for each
# plant - Agrilus pair and if == 1 => 0, 2 => 1)
interaction_data <- ddply(interaction_data, .(plant.sp, agrilus.sp),  nrow)
colnames(interaction_data) <- c("quercus.sp", "agrilus.sp", "interaction")
interaction_data$interaction <- interaction_data$interaction - 1
nrow(interaction_data) == length(agrilus)*length(oaks)



#
####       CREATE A DF WITH INFO ON OAK - ADD GEO OVERLAP DATA                      #####
# A) ADD DISTANCE INFORMATION
# For each plant sp. - Agrilus sp. pair, find the mean min. distance to the nearest
# host (from a different plant sp.)  of the current Agrilus sp.

# NB: "geodist() accepts only one or two primary arguments, which must be rectagular
# objects with unambiguously labelled longitude and latitude columns (i.e., some
# variant of lon/lat, or x/y)." Input(s) "can be in arbitrary rectangular format."
# See cran.r-project.org/web/packages/geodist/vignettes/geodist.html

# NB: Geodesic distance = "a curve representing the shortest path between two points"
# measure = "haversine" ("geodesic" is too computationally demanding)

# Add future columns to interaction_data
interaction_data$min.dist.mean <- rep(NA, nrow(interaction_data))
interaction_data$min.dist.mean.norm <- rep(NA, nrow(interaction_data))
interaction_data$min.dist.mean.log <- rep(NA, nrow(interaction_data))
interaction_data$min.dist.median <- rep(NA, nrow(interaction_data))
interaction_data$min.dist.median.norm <- rep(NA, nrow(interaction_data))
interaction_data$min.dist.median.log <- rep(NA, nrow(interaction_data))
str(interaction_data)

# Loop is not the most efficient: takes ca. 24 h to run (ask for 5 GB/iteration)
# Most iterations take < 1 h but some (30 or so??) take longer
# THIS THIS TO BE RUN INSIDE A JOBSCRIPT WITH 236 ARRAYS
oak <- as.numeric(commandArgs(trailingOnly=TRUE))
oak <- as.character(unique(interaction_data$quercus.sp)[oak])
first <- grep(oak, interaction_data$quercus.sp)[1]
print(oak)
print(nrow(interaction_data))

# For each oak - Agrilus pair (one oak at a time)
options(warn = 2)
for(i in first:(first+31)) {  # 1:nrow(interaction_data)
  # Get current plant sp
  oak <- interaction_data$quercus.sp[i]

  # Get current agrilus sp.
  agrilus <- interaction_data$agrilus.sp[i]

  # Find all host spp. (that are not the current oak sp., if it happens to be a host)
  hosts <- subset(agrilus_hosts, agrilus.sp == agrilus &
                    plant.sp != oak)$plant.sp
  print(paste0(i, ": ", oak, " - ", agrilus), )

  # If there are no hosts for the current Agrilus sp. other than the current oak sp.,
  # let's set min.dist.mean to 30k km (Earth's circumference = ca. 40k km)
  if (length(hosts) == 0) {
    min.dist.mean <- 3e4
    min.dist.median <- 3e4
  }
  # Otherwise
  else {
    # Extract geo info
    # i) For all individuals from the current oak sp.
    oak_geo <- subset(geo_info, plant.sp == oak)

    # ii) For all [other] plant host spp. of the current Agrilus sp.
    hosts_geo <- subset(geo_info, plant.sp %in% hosts)

    # Now, for all individuals from this oak sp. with coord. info from GBIF, find the
    # closest host in hosts_geo
    min.dists <- c()

    for (j in 1:nrow(oak_geo)){
      # Compute min dist (in km) to find the closest host in hosts_geo
      # Suppressed message when using "cheap" measure: "Maximum distance is > 100km.
      # The 'cheap' measure is inaccurate over such large distances, you'd likely be
      # better using a different 'measure'."
      suppressMessages(
        dists <- geodist(x = oak_geo[j,], y = hosts_geo, measure = "cheap")/1e3
      )
      min.dist <- min(dists)
      min.dists <- c(min.dists, min.dist)
    }

    # For this oak sp., find the mean min.dist and median min.dist
    min.dist.mean <- mean(min.dists)
    min.dist.median <- median(min.dists)
  }

  # Append the info for the current oak.sp to interaction.data
  interaction_data[i,]$min.dist.mean <- min.dist.mean
  interaction_data[i,]$min.dist.median <- min.dist.median
}
options(warn=1)

# Now transform distances
# First, multiply values * normal distr so that things that are close have a stronger
# weight. Use dnorm(), with mean = 0 and SD = 50 so that pairs that are > 100 km
# km away have very low values. The distr is rescaled so that scores look prettier.
# NB, dnorm(0, 0, 50) [max val] and dnorm(1e5, 0, 50) [min val] are just to make
# sure that all pairs are rescaled equally
# Mean
l <- length(interaction_data$min.dist.mean)
min.dist.mean.norm <- scales::rescale(-c(dnorm(0, 0, 50),
                                         dnorm(1e5, 0, 50),
                                         dnorm(interaction_data$min.dist.mean,
                                               0, 50)))[3:(l+2)]
interaction_data$min.dist.mean.norm <- min.dist.mean.norm; rm(min.dist.mean.norm, l)

# Median
l <- length(interaction_data$min.dist.median)
min.dist.median.norm <- scales::rescale(-c(dnorm(0, 0, 50),
                                           dnorm(1e5, 0, 50),
                                           dnorm(interaction_data$min.dist.median,
                                                 0, 50)))[3:(l+2)]
interaction_data$min.dist.median.norm <- min.dist.median.norm
rm(min.dist.median.norm, l)

# Now, log distances to check if a log distribution fits the data better
# Mean
interaction_data$min.dist.mean.log <- scales::rescale(log(
  interaction_data$min.dist.mean + 0.1))
# Median
interaction_data$min.dist.median.log <- scales::rescale(log(
  interaction_data$min.dist.median + 0.1))

# Finally, change "Quercus" for "Q.":
interaction_data$quercus.sp <-  gsub("Quercus", "Q.",
                                     interaction_data$quercus.sp,
                                     fixed = TRUE)


# Write table by appending (only add header if file doesn't already exist)
interaction_data <- interaction_data[first:(first+31), ]
nrow(interaction_data)

# write.table(x = interaction_data,
#             file = "tmp/interaction_data_oak_hosts.txt",
#             row.names = FALSE,
#             col.names = !file.exists("tmp/interaction_data_oak_hosts.txt"),
#             sep = "\t", quote = FALSE, append = TRUE)
print("done")                                                       # END OF JOBSCRIPT

# # NB: Note that the distance calculated is, inf fact, very different for large dists
# # when using the 'cheap' metric.
# # E.g., in our DF, the biggest imputed mean min. dist is ca. 23k km for Q. litoralis
# # and A. chiricahuae (only host = Q. hypoleucoides). However, if we use "harvestine"
# # or "geodesic" as our measure, we get ca. 14k km. Still, the transformed distances
# # performed better in our models (see analyses below), so this is not likely to be
# # something we'd need to worry about...
# ggplot() +
#   borders("world", colour = "gray50", fill = "gray50") +
#   # coord_cartesian(xlim = c(-10, 50), ylim = c(37, 70)) +
#   geom_point(aes(lon, lat),
#              data = geo_info[geo_info$plant.sp == "Quercus litoralis",]) +
#   geom_point(aes(lon, lat),
#              data = geo_info[geo_info$plant.sp == "Quercus hypoleucoides",],
#              col = "red")
#
# dists <- geodist(x = geo_info[geo_info$plant.sp == "Quercus litoralis",],
#                  y = geo_info[geo_info$plant.sp == "Quercus hypoleucoides",],
#                  measure = "haversine")/1e3
# m <- c()
# for (i in 1:nrow(dists)) {
#   m <- c(m, min(dists[i, ]))
# }
# mean(m); rm(dists, m)

# # NB2: Also note that, when using measure = "cheap", the distance values depend on
# # the number of locations that the distances are being computed for, e.g.: [note that
# # this doesn't happen if measure is set to "haversine" or "geodesic"]
# geodist(x = geo_info[geo_info$plant.sp == "Quercus uxoris",][1,],
#         y = geo_info[geo_info$plant.sp  == "Quercus imbricaria", ][1,],
#         measure = "cheap")/1e3
#
# geodist(x = geo_info[geo_info$plant.sp == "Quercus uxoris",][1,],
#         y = geo_info[geo_info$plant.sp  == "Quercus imbricaria", ][1:5,],
#         measure = "cheap")/1e3


# B) EXPLORE THE DATA
interaction_data <- read.table("tmp/interaction_data_oak_hosts.txt",
                               header = T, sep = "\t")

# i) Have a look at min.dist.mean and min.dist.mean
# min.dist.mean
min(interaction_data$min.dist.mean)
subset(interaction_data, min.dist.mean == min(interaction_data$min.dist.mean))
hist(interaction_data$min.dist.mean[interaction_data$interaction == 1],
     seq(0, 3e04, length.out = 30), main = "mean(min.dist) hosts - another host sp.
     of the current Agrilus sp.", xlab = "mean(dist.min)")
hist(interaction_data$min.dist.mean[interaction_data$interaction == 0],
     seq(0, 3e04, length.out = 30), main = "mean(min.dist) non-hosts - host sp.
     of the current Agrilus sp.", xlab = "mean(dist.min)")

# min.dist.median
min(interaction_data$min.dist.median)
subset(interaction_data, min.dist.median == min(interaction_data$min.dist.median))
hist(interaction_data$min.dist.median[interaction_data$interaction == 1],
     seq(0, 3e04, length.out = 30), main = "median(min.dist) hosts - another host sp.
     of the current Agrilus sp.", xlab = "median(dist.min)")
hist(interaction_data$min.dist.median[interaction_data$interaction == 0],
     seq(0, 3e04, length.out = 30), main = "median(min.dist) non-hosts - host sp.
     of the current Agrilus sp.", xlab = "median(dist.min)")


# ii) Explore max. difference between 'mean' and 'median' distances
# Plot histogram for min.dist.mean - min.dist.median
hist(interaction_data$min.dist.mean - interaction_data$min.dist.median,
     breaks = 100)

# Plot histogram for $min.dist.mean.norm - min.dist.median.norm
hist(interaction_data$min.dist.mean.norm - interaction_data$min.dist.median.norm,
     breaks = 100)

# Check max. difference between mean and median min. distances [could also do min.]
# Turns out it's for Q. acutissima [mostly present in NE America and SE Asia] - A.
# albocomus (hosts = Q. emoryi & Q. grisea [both NW America]) [in this case, preds
# using median perform better than w/ mean]
interaction_data[which.max(interaction_data$min.dist.mean -
                             interaction_data$min.dist.median), ]

# Plot distribution of hosts
ggplot(geo_info[geo_info$plant.sp == "Quercus emoryi" |
                  geo_info$plant.sp == "Quercus grisea",],
       aes(x = lon, y = lat, colour = plant.sp)) +
  coord_fixed() + borders("world", colour = "lightgrey", fill = "lightgrey") +
  geom_point() +
  scale_color_manual(values = c("Quercus emoryi" = scales::alpha('darkblue', 1),
                                'Quercus grisea' = scales::alpha('coral', 0.3))) +
  # geom_point(data = poe_geo, aes(x = lon, y = lat),
  # colour = scales::alpha('darkblue', 1)) +
  labs(colour = "Quercus spp.") +
  theme_bw() +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))

# Have a look at the min. distances
dists <- geodist(x = subset(geo_info, plant.sp == "Quercus acutissima"),
                 y = subset(geo_info, plant.sp == "Quercus emoryi" |
                              plant.sp == "Quercus grisea"),
                 measure = "haversine")/1e3
min_d <- c()
for (i in 1:nrow(dists)) {
  min_d <- c(min_d, min(dists[i,]))
}

# Plot a histogram to see what's going on
hist(min_d, breaks = 50) # makes sense
abline(v = interaction_data[which.max(interaction_data$min.dist.mean -
                                        interaction_data$min.dist.median),
                            ]$min.dist.mean,
       col = "red", lwd = 2, lty = "dashed")
abline(v = interaction_data[which.max(interaction_data$min.dist.mean -
                                        interaction_data$min.dist.median),
                            ]$min.dist.median,
       col = "blue", lwd = 2, lty = "dashed")


# iii) Explore max. difference between norm 'mean' and norm 'median' distances
# Check max. diff. between norm mean and median min. dists [could also do min.]
# Turns out it's for Q. acutissima [mostly present in NE America and SE Asia] - A.
# osburni (hosts = Q. virginiana [NE America]  & Q. phellos [W Europe & NE America])
interaction_data[which.max(interaction_data$min.dist.mean.norm -
                             interaction_data$min.dist.median.norm), ]

# Plot geo location of Q. acutissima (black) and hosts (blue)
plot(geo_info[geo_info$plant.sp == "Quercus acutissima", 2:3],
     xlim = c(-180, 180), ylim = c(-90, 90), pch = 16)
points(geo_info[geo_info$plant.sp == "Quercus virginiana", 2:3],
       col = scales::alpha("steelblue", 0.8), pch = 8)
points(geo_info[geo_info$plant.sp == "Quercus phellos", 2:3],
       col = scales::alpha("steelblue", 0.8), pch = 8)

# Have a look at the min. distances
dists <- geodist(x = subset(geo_info, plant.sp == "Quercus acutissima"),
                 y = subset(geo_info, plant.sp == "Quercus virginiana" |
                              plant.sp == "Quercus phellos"),
                 measure = "haversine")/1e3
min_d <- c()
for (i in 1:nrow(dists)) {
  min_d <- c(min_d, min(dists[i,]))
}

# Plot a histogram to see what's going on
hist(min_d, breaks = 50) # makes sense
abline(v = interaction_data[which.max(interaction_data$min.dist.mean.norm -
                                        interaction_data$min.dist.median.norm),
                            ]$min.dist.mean,
       col = "red", lwd = 2, lty = "dashed")
abline(v = interaction_data[which.max(interaction_data$min.dist.mean.norm -
                                        interaction_data$min.dist.median.norm),
                            ]$min.dist.median,
       col = "blue", lwd = 2, lty = "dashed")


# iv) Have a look at A. ussuricola (prediction is better using the mean instead of
# the median)
model_predictions[model_predictions$agrilus.sp == "A. ussuricola" &
                    model_predictions$host.status == 1,]

min.dists <- c()
for (j in 1:nrow(geo_info[geo_info$plant.sp == "Quercus acutissima",])) {
  suppressMessages(
    dists <- geodist(x = geo_info[geo_info$plant.sp == "Quercus acutissima",][j,],
                     y = geo_info[geo_info$plant.sp == "Quercus serrata",],
                     measure = "cheap")/1e3
  )
  min.dist <- min(dists)
  min.dists <- c(min.dists, min.dist)
}
hist(min.dists, breaks = 50)     # bimodal
abline(v = mean(min.dists), col = "red")
abline(v = median(min.dists), col = "blue")



# C) PLOT THE 3 GEO DIST METRCIS TO COMPARE THEM
geo_dists <- data.frame(matrix(c(
                        sort(interaction_data$min.dist.mean, decreasing = F),
                        sort(interaction_data$min.dist.mean.norm, decreasing = T),
                        sort(interaction_data$min.dist.mean.log, decreasing = F)),
                        ncol = 3))
colnames(geo_dists) <- c("raw.dist", "norm.dist", "log.dist")


# i) min.dist.mean
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


# ii) min.dist.mean.norm
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



# iii) min.dist.mean.log
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



#
