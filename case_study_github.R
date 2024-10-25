rm(list = ls())

library(MASS)
library(mvtnorm)
library(salso)
library(mcclust)
library(maps)
library(progress)
library(shiny)
library(shinythemes)
library(leaflet)

setwd("C:/Dati/Dottorato/Aiello_Argiento_Legramanti_Paci")

# MCMC algorithm
source("code/tseriesclust_similarity.R")
# Auxiliary functions necessary for the algorithm
source("code/auxiliary_functions/scaleandperiods.R")
source("code/auxiliary_functions/designmatrices.R")
source("code/auxiliary_functions/comp11.R")
source("code/auxiliary_functions/similarity.R")

# Time series data
pollutant <- read.csv('code/data/timeSeriesData.csv') 
pollutant <- pollutant[2191:2290,]

stationInfo <- read.csv('code/data/stationsInfo.csv')
##get latitude and longitude
latitude <- c()
longitude <- c()
for (i in 1:dim(pollutant)[2]) {
  site <- colnames(pollutant)[i]
  latitude <- c(latitude,stationInfo$latitude[which(stationInfo$site==site)])
  longitude <- c(longitude,stationInfo$longitude[which(stationInfo$site==site)])
}

locations <- data.frame(lat=latitude,lon=longitude)

plot(longitude,latitude,pch=19,cex=0.9,xlab='Longitude',ylab='Latitude')
map("world",add=T)

library(tictoc)

# 43896.79 seconds
tic()
pb=progress_bar$new(total=5000)
invisible(pb$tick(0))
tseriesc_similarity.out <- tseriesclust_similarity(pollutant,locations,maxiter=5000,burnin=2000,thinning=1,
                                                   frequency = 365,seasonfreq = 4,seasondelay = 79)
toc()

load("C:/Dati/Dottorato/Aiello_Argiento_Legramanti_Paci/code/results/result_fix_similarity_25_10_24.RData")

library(leaflegend)
library(htmlwidgets)
library(mcclust.ext)

cts <- scaleandperiods(pollutant,TRUE)$cts  # were removed from the original data set because they were constant.
locations <- locations[-cts,]

# Binder loss function

gnstar_bin <- as.numeric(salso(tseriesc_similarity.out$memorygn, maxNClusters=4, loss="binder", nRuns=50, nCores=4))

d <- dim(tseriesc_similarity.out$rhosample)[1]

rho <- tseriesc_similarity.out$rhosample[d,]
rho_m <- c()
tau <- tseriesc_similarity.out$sig2thesample[d,]
tau_m <- c()
for (i in 1:length(levels(as.factor(gnstar_bin)))) {
  rho_m <- c(rho_m, round(mean(rho[which(gnstar_bin==i)]),2))
  tau_m <- c(tau_m, round(mean(tau[which(gnstar_bin==i)]),2))
}

locations.clust <- data.frame(lon = locations$lon, lat = locations$lat, cluster = gnstar_bin, rho = rho_m[gnstar_bin], tau = tau_m[gnstar_bin])
locations.clust$gamma <- factor(paste("(",locations.clust$rho,",", locations.clust$tau,")", sep = ""))

pal <- colorFactor(c("blue", "green", "red","black"), domain = levels(as.factor(locations.clust$gamma)),ordered = FALSE)

# Create leaflet map
leaflet(data = locations.clust, options = list(zoomControl = FALSE,attributionControl=FALSE)) %>%
  addProviderTiles("Stadia.AlidadeSmooth") %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, radius = 5, stroke = FALSE,
                   fillColor = ~pal(gamma), fillOpacity = 0.9) %>%
  addLegendFactor(pal = pal, values = ~factor(gamma),
                  opacity = 0.9, shape = 'circle', title = HTML("<span style='font-size: 24px;'>&gamma;</span>"),  #HTML("&gamma;"), 
                  position ="topleft", width = 10, height = 10,
                  labelStyle = 'font-size: 18px;') 

table(gnstar_bin)

# similarity matrix

library(RColorBrewer)

simm <- comp.psm(tseriesc_similarity.out$memorygn)

# Assuming simm_bin and z_bin are already created:
simm_bin <- simm[order(gnstar_bin), order(gnstar_bin)]
z_bin <- gnstar_bin[order(gnstar_bin)]

library(reshape2)
# Convert the matrix into a long-format data frame
simm_bin_df <- melt(simm_bin)

# Define the color palette
grey_palette <- brewer.pal(8, "Greys")

# Create a data frame for the row side colors
row_side_colors_df <- data.frame(Var1 = 1:length(z_bin), cluster = factor(z_bin))

# Create a data frame for the column side colors
col_side_colors_df <- data.frame(Var2 = 1:length(z_bin), cluster = factor(z_bin))

# Create a color vector for row and column side colors based on z_bin
side_colors <- c("green", "red", "black", "blue")[z_bin]

# Create the heatmap with side color bars
ggplot(simm_bin_df, aes(Var1, Var2, fill = value)) + 
  # Heatmap tiles
  geom_tile() +
  # Reverse the Y axis for the 'revC' effect
  scale_y_reverse() +
  # Use grayscale for the heatmap colors (continuous scale)
  scale_fill_gradientn(colors = grey_palette, limits = c(min(simm_bin_df$value), max(simm_bin_df$value))) +
  # Add row side colors using `fill_side` aesthetic for the side bars
  geom_tile(data = row_side_colors_df, aes(x = -2.5, y = Var1), fill = side_colors, width = 5) +
  # Add column side colors using `fill_side` aesthetic for the side bars
  geom_tile(data = col_side_colors_df, aes(x = Var2, y = -2.5), fill = side_colors, height = 5) +
  # Ensure equal aspect ratio for the tiles
  coord_fixed() +
  # Customize the theme and remove unnecessary elements
  theme_void() +
  labs(x = "", y = "", fill = "Co-clustering \nprobability") +
  theme(axis.text = element_blank(),
        legend.title = element_text(size=20), legend.text = element_text(size=15), legend.key.size = unit(1, "cm"))

# VI loss function

gnstar_VI <- as.numeric(salso(tseriesc_similarity.out$memorygn, maxNClusters=4, loss=salso::VI(a=0.5), nRuns=50, nCores=4))

d <- dim(tseriesc_similarity.out$rhosample)[1]

rho <- tseriesc_similarity.out$rhosample[d,]
rho_m <- c()
tau <- tseriesc_similarity.out$sig2thesample[d,]
tau_m <- c()
for (i in 1:length(levels(as.factor(gnstar_VI)))) {
  rho_m <- c(rho_m, round(mean(rho[which(gnstar_VI==i)]),2))
  tau_m <- c(tau_m, round(mean(tau[which(gnstar_VI==i)]),2))
}

locations.clust <- data.frame(lon = locations$lon, lat = locations$lat, cluster = gnstar_VI, rho = rho_m[gnstar_VI], tau = tau_m[gnstar_VI])
locations.clust$gamma <- factor(paste("(",locations.clust$rho,",", locations.clust$tau,")", sep = ""))

pal <- colorFactor(c("blue", "green", "red","black"), domain = levels(as.factor(locations.clust$gamma)),ordered = FALSE)

# Create leaflet map
leaflet(data = locations.clust, options = list(zoomControl = FALSE,attributionControl=FALSE)) %>%
  addProviderTiles("Stadia.AlidadeSmooth") %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, radius = 5, stroke = FALSE,
                   fillColor = ~pal(gamma), fillOpacity = 0.9) %>%
  addLegendFactor(pal = pal, values = ~factor(gamma),
                  opacity = 0.9, shape = 'circle', title = HTML("<span style='font-size: 24px;'>&gamma;</span>"),  #HTML("&gamma;"), 
                  position ="topleft", width = 10, height = 10,
                  labelStyle = 'font-size: 18px;') 

table(gnstar_VI)


# similarity matrix

# Assuming simm_VI and z_VI are already created:
simm_VI <- simm[order(gnstar_VI), order(gnstar_VI)]
z_VI <- gnstar_VI[order(gnstar_VI)]

library(reshape2)
# Convert the matrix into a long-format data frame
simm_VI_df <- melt(simm_VI)

# Define the color palette
grey_palette <- brewer.pal(8, "Greys")

# Create a data frame for the row side colors
row_side_colors_df <- data.frame(Var1 = 1:length(z_VI), cluster = factor(z_VI))

# Create a data frame for the column side colors
col_side_colors_df <- data.frame(Var2 = 1:length(z_VI), cluster = factor(z_VI))

# Create a color vector for row and column side colors based on z_VI
side_colors <- c("red", "black", "green", "blue")[z_VI]

# Create the heatmap with side color bars
ggplot(simm_VI_df, aes(Var1, Var2, fill = value)) + 
  # Heatmap tiles
  geom_tile() +
  # Reverse the Y axis for the 'revC' effect
  scale_y_reverse() +
  # Use grayscale for the heatmap colors (continuous scale)
  scale_fill_gradientn(colors = grey_palette, limits = c(min(simm_VI_df$value), max(simm_VI_df$value))) +
  # Add row side colors using `fill_side` aesthetic for the side bars
  geom_tile(data = row_side_colors_df, aes(x = -2.5, y = Var1), fill = side_colors, width = 5) +
  # Add column side colors using `fill_side` aesthetic for the side bars
  geom_tile(data = col_side_colors_df, aes(x = Var2, y = -2.5), fill = side_colors, height = 5) +
  # Ensure equal aspect ratio for the tiles
  coord_fixed() +
  # Customize the theme and remove unnecessary elements
  theme_void() +
  labs(x = "", y = "", fill = "Co-clustering \nprobability") +
  theme(axis.text = element_blank(),
        legend.title = element_text(size=20), legend.text = element_text(size=15), legend.key.size = unit(1, "cm"))
