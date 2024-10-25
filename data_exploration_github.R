rm(list = ls())

setwd("C:/Dati/Dottorato/Aiello_Argiento_Legramanti_Paci/code")

pollutant_all <- read.csv("data/timeSeriesData.csv") 
pollutant <- pollutant_all[2191:2290,]

stationInfo <- read.csv("data/stationsInfo.csv")

##get latitude and longitude
latitude <- c()
longitude <- c()
region <- c()
for (i in 1:dim(pollutant)[2]) {
  site <- colnames(pollutant)[i]
  latitude <- c(latitude,stationInfo$latitude[which(stationInfo$site==site)])
  longitude <- c(longitude,stationInfo$longitude[which(stationInfo$site==site)])
}

locations <- data.frame(lat=latitude,lon=longitude)

plot(longitude,latitude,pch=19,cex=0.9,xlab='Longitude',ylab='Latitude')
map("world",add=T)

# plotting all the stations mean and sd ----------------------------------------

PM10 <- pollutant

PM10_summary <- cbind(stationInfo$longitude,stationInfo$latitude,colMeans(PM10),sapply(PM10, sd))
colnames(PM10_summary) <- c("Longitude","Latitude","avg","sd")
PM10_summary <- as.data.frame(PM10_summary)

library(leaflet)

pal_avg <- colorNumeric(
  palette = colorRampPalette(c('green4', 'yellow', 'orange', 'red3'))(dim(PM10_summary)[1]), 
  domain = PM10_summary$avg)

leaflet(PM10_summary, options = list(zoomControl = FALSE,attributionControl=FALSE)) %>% 
  addProviderTiles("Stadia.AlidadeSmooth") %>%
  addCircleMarkers(lng = ~Longitude, lat = ~Latitude,
                   stroke = FALSE, radius = 5, fillOpacity = 0.9, color = ~pal_avg(avg)) %>% 
  addLegend('topleft', title = "PM10 mean", pal = pal_avg, values = PM10_summary$avg, opacity = 0.9)

pal_sd <- colorNumeric(
  palette = colorRampPalette(c('green4', 'yellow', 'orange', 'red3'))(dim(PM10_summary)[1]), 
  domain = PM10$sd)

leaflet(PM10_summary, options = list(zoomControl = FALSE,attributionControl=FALSE)) %>% 
  addProviderTiles("Stadia.AlidadeSmooth") %>%
  addCircleMarkers(lng = ~Longitude, lat = ~Latitude,
                   stroke = FALSE, radius = 5, fillOpacity = 0.9, color = ~pal_sd(sd)) %>% 
  addLegend('topleft', title = "PM10 sd", pal = pal_sd, values = PM10_summary$sd, opacity = 0.9)

################################################################################

PM10_boxplot <- data.frame(q10 = apply(PM10,1,quantile, 0.10),
                           q25 = apply(PM10,1,quantile, 0.25),
                           q50 = apply(PM10,1,quantile, 0.50),
                           q75 = apply(PM10,1,quantile, 0.75),
                           q90 = apply(PM10,1,quantile, 0.90))

days <- seq(as.Date("2019-01-01"), as.Date("2019-04-10"), by = "days")

ggplot(PM10_boxplot) +
  geom_ribbon(aes(x = days, ymin = q10, ymax = q90), fill = "yellow2", alpha= 0.5) +
  geom_ribbon(aes(x = days, ymin = q25, ymax = q75), fill = "orange", alpha= 0.75) +
  geom_line(aes(x = days, y = q50), color = "red", linewidth = 1.5) +
  labs(y = expression("PM10 ("*mu*g/m^3*")"),
       color = "Quantile") +
  theme_minimal() +
  theme(axis.title = element_text(size=45),
        axis.text.y = element_text(size=30),
        axis.text.x = element_text(size=30, angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        strip.text = element_text(size=45)) +
  scale_x_date(date_labels = "%d %b", 
               breaks = seq(as.Date("2019-01-01"), as.Date("2019-04-10"), by = "1 week"))

################################################################################

# Define the maximum daily value and maximum dangerous days
max_daily_val <- 50  # critical daily level
max_dang_days <- 35  # threshold for the number of dangerous days

# Store the number of "dangerous" days per station
n_dang_days <- colSums(pollutant_all[2191:2555,] > max_daily_val)

# Create a data frame for ggplot
df <- data.frame(n_dang_days = n_dang_days)

# Plot the histogram using ggplot2
ggplot(df, aes(x = n_dang_days)) +
  geom_histogram(binwidth = 5, center = 2.5, fill = "skyblue", color = "black") +
  geom_vline(xintercept = max_dang_days, color = "red", linetype = "dashed", linewidth = 4) +
  labs(
    x = "Number of days above the limit",
    y = "Number of stations"
  ) +
  theme_bw() +
  theme(title = element_text(size=45),
        axis.title = element_text(size=45), 
        axis.text = element_text(size=45))

################################################################################

PM10_arima <- PM10[,-c(54,64)] # removing constant time series
locations_arima <- locations[-c(54,64),]

# Construction of the design matrix, the same used in  the final model 
source("auxiliary_functions/designmatrices.R")
DM <- designmatrices(deg = 2, 100, frequency = 365, seasonfreq = 2, seasondelay = 79)     
Z <- DM$Z

b0 <- c()
rho <- c()
sigma <- c()

for(r in 1:dim(PM10_arima)[2]){
  a <- arima(as.numeric(PM10_arima[,r]), order=c(1,0,0), xreg = Z, include.mean = FALSE)
  b0[r] <- as.numeric(a$coef[2])
  rho[r] <- as.numeric(a$coef[1])
  sigma[r] <- sqrt(as.numeric(a$sigma2))
}

arima_output <- data.frame(Longitude = locations_arima$lon, Latitude = locations_arima$lat,
                           b0 = b0, rho = rho, sigma = sigma)

library(leaflet)

pal_b0 <- colorNumeric(
  palette = colorRampPalette(c('green4', 'yellow', 'orange', 'red3'))(dim(arima_output)[1]), 
  domain = arima_output$b0)

leaflet(arima_output, options = list(zoomControl = FALSE,attributionControl=FALSE)) %>% 
  addProviderTiles("Stadia.AlidadeSmooth") %>%
  addCircleMarkers(lng = ~Longitude, lat = ~Latitude,
                   stroke = FALSE, radius = 5, fillOpacity = 0.9, color = ~pal_b0(b0)) %>% 
  addLegend('topleft', title = "Mean", pal = pal_b0, values = arima_output$b0, opacity = 0.9)

pal_rho <- colorNumeric(
  palette = colorRampPalette(c('green4', 'yellow', 'orange', 'red3'))(dim(arima_output)[1]), 
  domain = arima_output$rho)

leaflet(arima_output, options = list(zoomControl = FALSE,attributionControl=FALSE)) %>% 
  addProviderTiles("Stadia.AlidadeSmooth") %>%
  addCircleMarkers(lng = ~Longitude, lat = ~Latitude,
                   stroke = FALSE, radius = 5, fillOpacity = 0.9, color = ~pal_rho(rho)) %>% 
  addLegend('topleft', title = "Autocorrelation", pal = pal_rho, values = arima_output$rho, opacity = 0.9)

pal_sigma <- colorNumeric(
  palette = colorRampPalette(c('green4', 'yellow', 'orange', 'red3'))(dim(arima_output)[1]), 
  domain = arima_output$sigma)

leaflet(arima_output, options = list(zoomControl = FALSE,attributionControl=FALSE)) %>% 
  addProviderTiles("Stadia.AlidadeSmooth") %>%
  addCircleMarkers(lng = ~Longitude, lat = ~Latitude,
                   stroke = FALSE, radius = 5, fillOpacity = 0.9, color = ~pal_sigma(sigma)) %>% 
  addLegend('topleft', title = "Residual sd", pal = pal_sigma, values = arima_output$sigma, opacity = 0.9)
