rm(list = ls())

pollutant_all <- read.csv('data/timeSeriesData.csv') 
pollutant <- pollutant_all[2191:2555,]

# Define the day ranges for each season
winter_days <- c(1:79, 356:365)   # Jan 1 - Mar 20, Dec 22 - Dec 31
spring_days <- 80:171             # Mar 21 - Jun 20
summer_days <- 172:263            # Jun 21 - Sep 20
autumn_days <- 264:355            # Sep 21 - Dec 21

# Create a vector for the seasons
seasons <- character(365)
seasons[winter_days] <- "Winter"
seasons[spring_days] <- "Spring"
seasons[summer_days] <- "Summer"
seasons[autumn_days] <- "Fall"

# Convert to a factor with "Winter" as the reference level
season_factor <- factor(seasons, levels = c("Winter", "Spring", "Summer", "Fall"))

# Generate dummy variables (excluding Winter as baseline)
Z <- model.matrix(~ season_factor)[, -1]

Z <- cbind(rep(1,365), Z)

# Rename columns for clarity
colnames(Z) <- c("Intercept", "Spring", "Summer", "Fall")

# Print first few rows
head(Z)

stationInfo <- read.csv('data/stationsInfo.csv')

library(Rcpp)
library(RcppArmadillo)

sourceCpp('samplers/sampler_sPPM.cpp')

# initial values

set.seed(123)

n <- dim(pollutant)[2]
p <- dim(Z)[2]

init_sig2beta <- 1/rgamma(p, shape = 3, rate = 2) 

init_beta <- matrix(rnorm(n * p, 0, sqrt(init_sig2beta)), p, n) 
init_sigma2 <- 1/rgamma(n, shape = 3, rate = 2)

init_phi <- runif(n) # 
init_tau2 <- 1/rgamma(n, shape = 3, rate = 2) # 

init_M <- 2.212124 # posterior mean of the same parameter using PPM

pollutant_scaled <- scale(pollutant, scale = FALSE)

library(tictoc)

tic()
tseriesc_similarity.out <- tseriesclust_similarity_cpp(as.matrix(pollutant_scaled), as.matrix(locations), Z,
                                                       runs = 15000, burn = 10000, thin = 1, similarity = 3, 
                                                       init_beta = init_beta,
                                                       init_sigma2 = init_sigma2, init_sig2beta = init_sig2beta,
                                                       init_phi = init_phi, init_tau2 = init_tau2, init_M = init_M,
                                                       a_eps = 2, b_eps = 1,
                                                       a_beta = 2, b_beta = 1,
                                                       a_tau = 2, b_tau = 1,
                                                       a_phi = 1, b_phi = 1, eta_phi = 0.05)
toc()

library(leaflet)
library(leaflegend)
library(htmlwidgets)
library(htmltools)
library(mcclust.ext)
library(salso)

sPPM_clusters <- tseriesc_similarity.out$cluster

pollutant.out <- pollutant[ ,-tseriesc_similarity.out$rem_obs]
pollutant_scaled.out <- pollutant_scaled[ ,-tseriesc_similarity.out$rem_obs]
locations.out <- locations[-tseriesc_similarity.out$rem_obs, ]

similarity_clust_VI <- as.numeric(salso(tseriesc_similarity.out$cluster, maxNClusters=0, loss=salso::VI(a=1), nRuns=100, nCores=6))

init_phi <- runif(length(unique(similarity_clust_VI)))[similarity_clust_VI] # 
init_tau2 <- 1/rgamma(length(unique(similarity_clust_VI)), shape = 3, rate = 2)[similarity_clust_VI] #

init_phi <- runif(length(unique(similarity_clust_VI)))[similarity_clust_VI] # 
init_tau2 <- 1/rgamma(length(unique(similarity_clust_VI)), shape = 3, rate = 2)[similarity_clust_VI] #

# rerunning the algorithm with partition obtained theough the VI for cluster specific estimates 

tic()
tseriesc_similarity.out <- tseriesclust_similarity_cpp(as.matrix(pollutant_scaled), as.matrix(locations), Z,
                                                   runs = 2000, burn = 1000, thin = 1, similarity = 3, clust = FALSE,
                                                   init_beta = init_beta,
                                                   init_sigma2 = init_sigma2, init_sig2beta = init_sig2beta,
                                                   init_phi = init_phi, init_tau2 = init_tau2, init_M = init_M,
                                                   a_eps = 2, b_eps = 1,
                                                   a_beta = 2, b_beta = 1,
                                                   a_tau = 2, b_tau = 1,
                                                   a_phi = 1, b_phi = 1, eta_phi = 0.05,
                                                   a_M = 2, b_M = 0.5)
toc()

phi <- apply(tseriesc_similarity.out$phi, 2, median)
phi_m <- c()
tau2 <- apply(tseriesc_similarity.out$tau2, 2, median)
tau2_m <- c()
for (i in 1:length(levels(as.factor(similarity_clust_VI)))) {
  phi_m <- c(phi_m, round(mean(phi[which(similarity_clust_VI==i)]),2))
  tau2_m <- c(tau2_m, round(mean(tau2[which(similarity_clust_VI==i)]),2))
}

locations.clust <- data.frame(lon = locations.out$lon, lat = locations.out$lat, cluster = similarity_clust_VI, phi = phi_m[similarity_clust_VI], tau2 = tau2_m[similarity_clust_VI])
locations.clust$gamma <- factor(paste("(",locations.clust$phi,",", locations.clust$tau2,")", sep = ""))

table(similarity_clust_VI)

# Create a color factor using the generated palette
pal <- colorFactor(c("green4", "green2", "yellow2", "orange2", "red2"), domain = levels(as.factor(locations.clust$gamma)),ordered = FALSE)

# Create leaflet map
leaflet(data = locations.clust, options = list(zoomControl = FALSE,attributionControl=FALSE)) %>%
  addProviderTiles("Stadia.AlidadeSmooth") %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, radius = 9, stroke = FALSE,
                   fillColor = ~pal(gamma), fillOpacity = 0.9) %>%
  addLegendFactor(pal = pal, values = ~factor(gamma),
                  opacity = 0.9, shape = 'circle', 
                  title = HTML("<span style='font-size: 30px;'>&theta;</span>"),
                  position ="topleft", width = 20, height = 20,
                  labelStyle = 'font-size: 20px;')

# coclustering matrix -----------------------------------------------------

simm <- comp.psm(sPPM_clusters)

heatmap(simm)

# Assuming simm_VI and z_VI are already created:
simm_VI <- simm[order(similarity_clust_VI), order(similarity_clust_VI)]
z_VI <- similarity_clust_VI[order(similarity_clust_VI)]

library(reshape2)
# Convert the matrix into a long-format data frame
simm_VI_df <- melt(simm_VI)

library(RColorBrewer)
# Define the color palette
grey_palette <- brewer.pal(8, "Greys")

# Create a data frame for the row side colors
row_side_colors_df <- data.frame(Var1 = 1:length(z_VI), cluster = factor(z_VI))

# Create a data frame for the column side colors
col_side_colors_df <- data.frame(Var2 = 1:length(z_VI), cluster = factor(z_VI))

# Create a color vector for row and column side colors based on z_VI
side_colors <- pal(factor(locations.clust$gamma))[order(similarity_clust_VI)]

library(ggplot2)
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
        legend.position = "none")


# ts comparison with bands ------------------------------------------------

library(ggplot2)
library(tidyr)
library(dplyr)

# Convert `pollutant.out` from wide to long format
df_long <- pollutant.out %>%
  mutate(time = 1:nrow(pollutant.out)) %>%  # Add a time variable
  pivot_longer(cols = -time, names_to = "Series", values_to = "Value")

# Create a data frame for cluster memberships
cluster_df <- data.frame(Series = colnames(pollutant.out), 
                         Cluster = factor(similarity_clust_VI),
                         gamma = factor(paste("(",phi_m[similarity_clust_VI],",", tau2_m[similarity_clust_VI],")", sep = "")))  # Ensure it's a factor

df_long <- df_long %>%
  left_join(cluster_df, by = "Series")

# Compute summary statistics per cluster and time point
df_summary <- df_long %>%
  group_by(time, gamma) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    lower_90 = quantile(Value, 0.05, na.rm = TRUE),
    upper_90 = quantile(Value, 0.95, na.rm = TRUE),
    lower_95 = quantile(Value, 0.025, na.rm = TRUE),
    upper_95 = quantile(Value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
ggplot(df_summary, aes(x = time, y = mean_value, color = gamma)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = gamma), alpha = 0.2, color = NA) +  # 95% interval
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = gamma), alpha = 0.4, color = NA) +  # 90% interval
  geom_line(linewidth = 0.66) +  # Mean line
  scale_color_manual(values = pal(df_summary$gamma)) +
  scale_fill_manual(values = pal(df_summary$gamma)) +  # Ensure fill colors match line colors
  facet_wrap(~ gamma) + 
  labs(x = "Day",
       y = "PM10",
       color = expression(gamma),
       fill = expression(gamma)
  ) +
  theme_minimal() + 
  theme(
    axis.title.x = element_text(size = 16),                 # X-axis label size
    axis.title.y = element_text(size = 16),                 # Y-axis label size
    axis.text.x = element_text(size = 14),                  # X-axis tick size
    axis.text.y = element_text(size = 14),                  # Y-axis tick size
    strip.text = element_text(size = 16),
    legend.position = "none",
    legend.key.width = unit(0.5, "cm")
  )

