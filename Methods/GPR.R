# This script processes flux data for three stars, imputes missing flux values,
# fits a Gaussian Process Regression (GPR) model to predict flux values, 
# calculates residuals, flags anomalies based on a defined threshold, 
# and visualizes anomalies on a time series plot.

# Load required libraries
library(tidyverse)
library(kernlab)

# Read in data
s1 <- read_csv("Data/031381302.csv")
s2 <- read_csv("Data/129646813.csv")
s3 <- read_csv("Data/0131799991.csv")

# Impute missing flux values using Kalman filter
s1$pdcsap_flux_imputed <- na_kalman(s1$pdcsap_flux, model = "StructTS")
s2$pdcsap_flux_imputed <- na_kalman(s2$pdcsap_flux, model = "StructTS")
s3$pdcsap_flux_imputed <- na_kalman(s3$pdcsap_flux, model = "StructTS")

# Commented out because it takes long to run:
# gpr_model <- gausspr(pdcsap_flux_imputed ~ time, data = s1, kernel = "rbfdot")
# save(gprModel, file = "gprModel.RData")

# Load the pre-trained GPR model
load("Data/gprModel.RData")
gprModel <- gpr_model

# Predict flux values using the GPR model
gprPred <- predict(gprModel, newdata = s1)

# Calculate residuals (difference between observed and predicted flux)
gprResiduals <- s1$pdcsap_flux_imputed - gprPred

# Define an anomaly threshold
threshold <- mean(gprResiduals, na.rm = TRUE) + 3 * sd(gprResiduals, na.rm = TRUE)

# Flag anomalies based on the residuals
s1$gpr_anomalies <- ifelse(gprResiduals > threshold, "Flare", "Normal")

# Create a data frame for plotting
plotData <- data.frame(time = s1$time, flux = s1$pdcsap_flux_imputed, flare = s1$gpr_anomalies)

# Visualize anomalies detected by GPR
ggplot(plotData, aes(x = time, y = flux)) +
  geom_point(aes(color = flare), size = 1) +
  scale_color_manual(values = c("Normal" = "black", "Flare" = "red")) +
  labs(title = "Gaussian Process Regression Anomalies",
       x = "Time", y = "PDCSAP Flux")
