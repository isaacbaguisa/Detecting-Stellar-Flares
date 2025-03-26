# This script processes flux data for three stars, imputes missing flux values,
# fits a Poisson Process Regression (PPR) model to predict flux values, 
# calculates residuals, flags anomalies based on a defined threshold, 
# and visualizes anomalies on a time series plot.

# Load required libraries
library(tidyverse)
library(glmnet)

# Read in data
s1 <- read_csv("Data/031381302.csv")
s2 <- read_csv("Data/129646813.csv")
s3 <- read_csv("Data/0131799991.csv")

# Impute missing flux values using Kalman filter
s1$pdcsap_flux_imputed <- na_kalman(s1$pdcsap_flux, model = "StructTS")
s2$pdcsap_flux_imputed <- na_kalman(s2$pdcsap_flux, model = "StructTS")
s3$pdcsap_flux_imputed <- na_kalman(s3$pdcsap_flux, model = "StructTS")

# Fit Poisson Process Regression model
pprModel <- glm(pdcsap_flux_imputed ~ time, family = poisson(link = "log"), 
                data = s1)

# Predicting the flux using the fitted model
pprPred <- predict(pprModel, newdata = s1, type = "response")

# Calculate residuals (difference between observed and predicted flux)
pprResiduals <- s1$pdcsap_flux_imputed - pprPred

# Define an anomaly threshold (residuals greater than 3 sd from the mean)
threshold <- mean(pprResiduals, na.rm = TRUE) + 3 * sd(pprResiduals, 
                                                       na.rm = TRUE)

# Flag anomalies based on the residuals
s1$ppr_anomalies <- ifelse(pprResiduals > threshold, "Flare", "Normal")

# Create a data frame for plotting
plotData <- data.frame(time = s1$time, flux = s1$pdcsap_flux_imputed, 
                       flare = s1$ppr_anomalies)

# Visualize anomalies detected by PPR using ggplot2
ggplot(plotData, aes(x = time, y = flux)) +
  geom_point(aes(color = flare), size = 1) +
  scale_color_manual(values = c("Normal" = "black", "Flare" = "red")) +
  labs(title = "Poisson Process Regression Anomalies", 
       x = "Time", y = "PDCSAP Flux")
