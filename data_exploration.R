rm(list = ls())

pollutant_all <- read.csv("data/timeSeriesData.csv") 
PM10 <- pollutant_all[2191:2555,] # year 2019

stationInfo <- read.csv("data/stationsInfo.csv")

# plotting all the stations -----------------------------------------------

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
                   stroke = FALSE, radius = 9, fillOpacity = 0.9, color = ~pal_avg(avg)) %>% 
  addLegend('topleft', title = "PM10 mean", pal = pal_avg, values = PM10_summary$avg, opacity = 0.9)

pal_sd <- colorNumeric(
  palette = colorRampPalette(c('green4', 'yellow', 'orange', 'red3'))(dim(PM10_summary)[1]), 
  domain = PM10$sd)

leaflet(PM10_summary, options = list(zoomControl = FALSE,attributionControl=FALSE)) %>% 
  addProviderTiles("Stadia.AlidadeSmooth") %>%
  addCircleMarkers(lng = ~Longitude, lat = ~Latitude,
                   stroke = FALSE, radius = 9, fillOpacity = 0.9, color = ~pal_sd(sd)) %>% 
  addLegend('topleft', title = "PM10 sd", pal = pal_sd, values = PM10_summary$sd, opacity = 0.9)

# Time series with bands --------------------------------------------------

PM10_boxplot <- data.frame(q05 = apply(PM10,1,quantile, 0.05),
                           q25 = apply(PM10,1,quantile, 0.25),
                           q50 = apply(PM10,1,quantile, 0.50),
                           q75 = apply(PM10,1,quantile, 0.75),
                           q95 = apply(PM10,1,quantile, 0.95))

days <- seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days")

ggplot(PM10_boxplot) +
  geom_ribbon(aes(x = days, ymin = q05, ymax = q95), fill = "yellow2", alpha= 0.5) +
  geom_ribbon(aes(x = days, ymin = q25, ymax = q75), fill = "orange", alpha= 0.75) +
  geom_line(aes(x = days, y = q50), color = "red", linewidth = 1.5) +
  labs(y = expression("PM10 ("*mu*g/m^3*")")) +
  theme_minimal() +
  theme(axis.title = element_text(size=30),
        axis.text.y = element_text(size=30),
        axis.text.x = element_text(size=30, angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        strip.text = element_text(size=45)) +
  scale_x_date(date_labels = "%d %b", 
               breaks = seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "2 weeks"))


# Maximum days above threshold --------------------------------------------

# Define the maximum daily value and maximum dangerous days

max_daily_val <- 50  # critical daily level
max_dang_days <- 35  # threshold for the number of dangerous days

# Store the number of "dangerous" days per station
n_dang_days <- colSums(PM10 > max_daily_val)

# Create a data frame for ggplot
df <- data.frame(n_dang_days = n_dang_days)

# Plot the histogram using ggplot2
ggplot(df, aes(x = n_dang_days)) +
  geom_histogram(binwidth = 5, center = 2.5, fill = "skyblue", color = "black") +
  geom_vline(xintercept = max_dang_days, color = "red", linetype = "dashed", linewidth = 1.5) +
  labs(
    x = "Number of days above the limit",
    y = "Number of stations"
  ) +
  theme_minimal() +
  theme(axis.title = element_text(size=30),
        axis.text = element_text(size=30))
