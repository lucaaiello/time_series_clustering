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

sourceCpp('samplers/sampler_noclust.cpp')

# initial values

set.seed(123)

n <- dim(pollutant)[2]
p <- dim(Z)[2]

init_sig2beta <- 1/rgamma(p, shape = 3, rate = 2) # 1

init_beta <- matrix(rnorm(n * p, 0, sqrt(init_sig2beta)), p, n) # matrix(0, p, n) 
init_sigma2 <- 1/rgamma(n, shape = 3, rate = 2) # rep(1, n) 

init_phi <- runif(n) # 
init_tau2 <- 1/rgamma(n, shape = 3, rate = 2) # 

pollutant_scaled <- scale(pollutant, scale = FALSE)

library(tictoc)

tic()
tseries.out <- tseries_cpp(as.matrix(pollutant_scaled), Z,
                           runs = 15000, burn = 10000, thin = 1, 
                           init_beta = init_beta,
                           init_sigma2 = init_sigma2, init_sig2beta = init_sig2beta,
                           init_phi = init_phi, init_tau2 = init_tau2, 
                           a_eps = 2, b_eps = 1,
                           a_beta = 2, b_beta = 1,
                           a_tau = 2, b_tau = 1,
                           a_phi = 1, b_phi = 1, eta_phi = 0.1)
toc()

# Output maps

pollutant.out <- pollutant[ ,-tseries.out$rem_obs]
pollutant_scaled.out <- pollutant_scaled[ ,-tseries.out$rem_obs]
locations.out <- locations[-tseries.out$rem_obs, ]

library(coda)
HPD_phi <- HPDinterval(mcmc(tseries.out$phi),prob = 0.01)
HPD_tau2 <- HPDinterval(mcmc(tseries.out$tau2),prob = 0.01)

phi <- round(apply(HPD_phi, 1, mean),2)
tau2 <- round(apply(HPD_tau2, 1, mean),2)

locations.clust <- data.frame(lon = locations.out$lon, lat = locations.out$lat, phi = phi, tau2 = tau2)
locations.clust$gamma <- factor(paste("(",locations.clust$phi,",", locations.clust$tau2,")", sep = ""))

pal_phi <- colorNumeric(
  palette = colorRampPalette(c("green4", "green2", "yellow2", "orange2", "red2"))(5), 
  domain = locations.clust$phi[-37])

leaflet(data = locations.clust[-37,], options = list(zoomControl = FALSE, attributionControl = FALSE)) %>%
  addProviderTiles("Stadia.AlidadeSmooth") %>%
  addCircleMarkers(
    lng = ~lon, lat = ~lat, radius = 9, stroke = FALSE,
    fillColor = ~pal_phi(phi), fillOpacity = 0.9
  ) %>%
  addLegend(
    pal = pal_phi, values = locations.clust$phi[-37],  # Ensure `phi` is properly referenced
    opacity = 0.9, # title = HTML("<span style='font-size: 20px;'>&phi;</span>"),
    title = HTML("<span style='font-size: 30px; font-family: Times New Roman, serif;'>φ</span>"),
    position = "topleft", labFormat = labelFormat(digits = 2)  # Format labels correctly
  )

pal_tau2 <- colorNumeric(
  palette = colorRampPalette(c("green4", "green2", "yellow2", "orange2", "red2"))(5), 
  domain = locations.clust$tau2)

leaflet(data = locations.clust, options = list(zoomControl = FALSE, attributionControl = FALSE)) %>%
  addProviderTiles("Stadia.AlidadeSmooth") %>%
  addCircleMarkers(
    lng = ~lon, lat = ~lat, radius = 9, stroke = FALSE,
    fillColor = ~pal_tau2(tau2), fillOpacity = 0.9
  ) %>%
  addLegend(
    pal = pal_tau2, values = locations.clust$tau2,  # Ensure `phi` is properly referenced
    opacity = 0.9, #title = HTML("<span style='font-size: 20px;'>&tau;</span>"),
    title = HTML("<span style='font-size: 30px; font-family: Times New Roman, serif;'>τ²</span>"),
    position = "topleft", labFormat = labelFormat(digits = 2)  # Format labels correctly
  )


