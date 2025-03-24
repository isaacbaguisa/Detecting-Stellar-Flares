# This script processes flux data for three stars, imputes missing flux values,
# fits a One-Class Support Vector Machine (OCSVM) model to detect anomalies, 
# calculates decision function scores, flags anomalies based on a defined 
# threshold, and visualizes anomalies on a time series plot.

# Load required libraries
library(tidyverse)
library(e1071)

# Read in data
s1 <- read_csv("Data/031381302.csv")
s2 <- read_csv("Data/129646813.csv")
s3 <- read_csv("Data/0131799991.csv")

# Impute missing flux values using Kalman filter
s1$pdcsap_flux_imputed <- na_kalman(s1$pdcsap_flux, model = "StructTS")
s2$pdcsap_flux_imputed <- na_kalman(s2$pdcsap_flux, model = "StructTS")
s3$pdcsap_flux_imputed <- na_kalman(s3$pdcsap_flux, model = "StructTS")

# Fit the One-Class SVM model with decision function enabled
ocsvmModel <- svm(s1$pdcsap_flux_imputed, type = "one-classification", nu = 0.01,
                  kernel = "radial", scale = TRUE, 
                  gamma = 0.05, decision.values = TRUE)

# Get the decision function scores (distance from the boundary)
decisionScores <- ocsvmModel$decision.values

# Define a manual threshold for anomalies
threshold <- -2000

# Flag anomalies based on the threshold
s1$ocsvm_anomalies <- ifelse(decisionScores < threshold, "Anomaly", "Normal")

# Create a data frame for plotting
plotData <- data.frame(time = s1$time, flux = s1$pdcsap_flux_imputed, 
                       anomalies = s1$ocsvm_anomalies)
colnames(plotData) <- c("time", "flux", "anomalies")

# Visualize anomalies detected by OCSVM with manual threshold using ggplot2
ggplot(plotData, aes(x = time, y = flux)) +
  geom_point(aes(color = anomalies), size = 1) +
  scale_color_manual(values = c("Normal" = "black", "Anomaly" = "red")) +
  labs(title = "One-Class SVM Anomalies with Manual Threshold", 
       x = "Time", y = "PDCSAP Flux")
