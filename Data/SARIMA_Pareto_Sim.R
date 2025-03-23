# This script simulates light curves for three stars by fitting a SARIMA model 
# to the flux data, generating simulated data, and injecting synthetic flares. 
# It also visualizes the actual and simulated light curves, highlighting the 
# synthetic flares injected following a Pareto distribution 

# @params starData: A data frame containing the time series data of a star, 
# including flux values (with missing values imputed).
# @params starName: A string representing the name of the star (used in the plot title).
# @returns A data frame containing the simulated light curve with injected 
# synthetic flares, including the time, simulated flux, and flare labels.

# Load required libraries
library(forecast)
library(ggplot2)
library(imputeTS)
library(tidyverse)
library(extraDistr)

set.seed(12345)

# Read in data
s1 <- read_csv("031381302.csv")
s2 <- read_csv("129646813.csv")
s3 <- read_csv("0131799991.csv")

# Impute missing flux values using Kalman filter
s1$pdcsap_flux_imputed <- na_kalman(s1$pdcsap_flux, model = "StructTS")
s2$pdcsap_flux_imputed <- na_kalman(s2$pdcsap_flux, model = "StructTS")
s3$pdcsap_flux_imputed <- na_kalman(s3$pdcsap_flux, model = "StructTS")

# Function to inject flares based on Pareto distribution
simulate_light_curve_pareto <- function(star_data, star_name) {
  
  # Fit the SARIMA model (or another model if needed)
  fit_sarima <- auto.arima(star_data$pdcsap_flux_imputed, seasonal = TRUE)
  
  # Simulate the light curve using the SARIMA model
  simulated_data <- simulate(fit_sarima, nsim = 
                               length(star_data$pdcsap_flux_imputed))
  
  # Initialize simulated dataframe
  simulated_df <- data.frame(
    time = 1:length(star_data$pdcsap_flux_imputed),
    pdcsap_flux_imputed = simulated_data,
    flare_label = rep(0, length(star_data$pdcsap_flux_imputed))  # Initially no flares
  )
  
  # Inject 10 synthetic flares using a Pareto distribution
  flare_indices <- sample(1:length(star_data$pdcsap_flux_imputed), 
                          size = 10, replace = FALSE)
  
  for (index in flare_indices) {
    # Use the Pareto distribution to determine flare amplitude
    flare_amplitude <- rpareto(1, a = 1, b = 50)  # Using parameters for large flares
    simulated_df$pdcsap_flux_imputed[index] <- simulated_df$pdcsap_flux_imputed[index] + 
      flare_amplitude
    simulated_df$flare_label[index] <- 1  # Mark flare
  }
  
  # Plot the actual data and simulated data
  p <- ggplot() +
    geom_line(data = star_data, aes(x = 1:length(star_data$pdcsap_flux_imputed), 
                                    y = star_data$pdcsap_flux_imputed), 
              color = "blue", alpha = 0.7) +  # Actual data
    geom_line(data = simulated_df, aes(x = time, y = pdcsap_flux_imputed), 
              color = "green", alpha = 0.7) +  # Simulated data
    geom_point(data = simulated_df[simulated_df$flare_label == 1, ], 
               aes(x = time, y = pdcsap_flux_imputed), color = "red") +  # Flares
    labs(title = paste("Actual vs. Simulated Light Curve for", star_name),
         x = "Time", y = "PDCSAP Flux") +
    theme_minimal() +
    scale_color_manual(values = c("blue", "green"))
  
  print(p)
  
  return(simulated_df)
}

# Simulate light curves for all stars and visualize actual vs. simulated data
simulatedS1 <- simulate_light_curve_pareto(s1, "Star 1")
simulatedS2 <- simulate_light_curve_pareto(s2, "Star 2")
simulatedS3 <- simulate_light_curve_pareto(s3, "Star 3")


# Simulate light curves for all stars and visualize actual vs. simulated data
png("simS1Pareto.png", width = 800, height = 600, res = 150)
simulatedS1 <- simulate_light_curve_pareto(s1, "Star 1")
dev.off()

png("simS2Pareto.png", width = 800, height = 600, res = 150)
simulatedS2 <- simulate_light_curve_pareto(s2, "Star 2")
dev.off()

png("simS3Pareto.png", width = 800, height = 600, res = 150)
simulatedS3 <- simulate_light_curve_pareto(s3, "Star 3")
dev.off()
















