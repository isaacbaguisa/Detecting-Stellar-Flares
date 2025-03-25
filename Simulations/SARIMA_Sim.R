# This script simulates light curves for three stars by fitting a SARIMA model 
# to the flux data, generating simulated data, and injecting synthetic flares. 
# It also visualizes the actual and simulated light curves, highlighting the 
# synthetic flares.

# @param starData: A data frame containing the time series data of a star, 
# including flux values (with missing values imputed).
# @param starName: A string representing the name of the star (used in the plot title).
# @return A data frame containing the simulated light curve with injected 
# synthetic flares, including the time, simulated flux, and flare labels.

# Load required libraries
library(forecast)
library(ggplot2)
library(imputeTS)
library(tidyverse)

set.seed(12345)

# Read in data
s1 <- read_csv("031381302.csv")
s2 <- read_csv("129646813.csv")
s3 <- read_csv("0131799991.csv")

# Impute missing flux values using Kalman filter
s1$pdcsap_flux_imputed <- na_kalman(s1$pdcsap_flux, model = "StructTS")
s2$pdcsap_flux_imputed <- na_kalman(s2$pdcsap_flux, model = "StructTS")
s3$pdcsap_flux_imputed <- na_kalman(s3$pdcsap_flux, model = "StructTS")

# Function to simulate light curves and inject synthetic flares
simulateLightCurve <- function(starData, starName) {
  fit <- auto.arima(starData$pdcsap_flux_imputed, seasonal = TRUE)
  simulatedData <- simulate(fit, nsim = length(starData$pdcsap_flux_imputed))
  
  # Convert to a data frame and keep the original columns 
  simulatedDf <- data.frame(
    time = 1:length(starData$pdcsap_flux_imputed),
    pdcsap_flux_imputed = simulatedData,
    flare_label = rep(0, length(starData$pdcsap_flux_imputed)) # No flares
  )
  
  # Inject 10 synthetic flares into the simulated data
  flareIndices <- sample(1:length(starData$pdcsap_flux_imputed), 
                         size = 10, 
                         replace = FALSE)
  
  for (index in flareIndices) {
    simulatedDf$pdcsap_flux_imputed[index] <- simulatedDf$pdcsap_flux_imputed[index] + 
      runif(1, min = 50, max = 350) # Add flare
    simulatedDf$flare_label[index] <- 1 # Mark flare
  }
  
  # Plot the actual data and the simulated data together
  p <- ggplot() +
    geom_line(data = starData, aes(x = 1:length(starData$pdcsap_flux_imputed), 
                                   y = starData$pdcsap_flux_imputed), 
              color = "blue", alpha = 0.7) +  # Actual data
    geom_line(data = simulatedDf, aes(x = time, y = pdcsap_flux_imputed), 
              color = "green", alpha = 0.7) +  # Simulated data
    geom_point(data = simulatedDf[simulatedDf$flare_label == 1, ], 
               aes(x = time, y = pdcsap_flux_imputed), color = "red") + 
    labs(title = paste("Actual vs. Simulated Light Curve for", starName),
         x = "Time", y = "PDCSAP Flux") +
    theme_minimal() +
    scale_color_manual(values = c("blue", "green"))
  
  print(p)
  
  return(simulatedDf)
}

# Simulate light curves for all stars and visualize actual vs. simulated data
png("simS1.png", width = 800, height = 600, res = 150)
simulatedS1 <- simulateLightCurve(s1, "Star 1")
dev.off()

png("simS2.png", width = 800, height = 600, res = 150)
simulatedS2 <- simulateLightCurve(s2, "Star 2")
dev.off()

png("simS3.png", width = 800, height = 600, res = 150)
simulatedS3 <- simulateLightCurve(s3, "Star 3")
dev.off()






