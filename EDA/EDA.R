# This code does an exploratory data analysis and processes stellar flare data 
# for three stars(TIC 031381302, TIC 129646813, and TIC 0131799991).
# It involves the following steps:
# 1. Reading and cleaning data, including the imputation of missing values 
# in the PDCSAP_FLUX.
# 2. Generating plots comparing original and imputed flux data.
# 3. Analyzing the distribution of PDCSAP_FLUX and identifying potential flares 
# using statistical thresholds.
# 4. Performing DBSCAN clustering on the flux data to identify clusters of 
# brightness events.
# 5. Using ARMA models to analyze residuals and assess model fit.

# Load libraries
library(knitr)
library(dplyr)
library(ggplot2)
library(dbscan)
library(tidyverse)
library(zoo)
library(gridExtra)
library(imputeTS)
library(forecast)

# Read in data
s1 <- read_csv("Data/031381302.csv")
s2 <- read_csv("Data/129646813.csv")
s3 <- read_csv("Data/0131799991.csv")

s1$star_id <- "TIC 031381302"
s2$star_id <- "TIC 129646813"
s3$star_id <- "TIC 0131799991"

# Summary of missing values
missingS1 <- table(is.na(s1$pdcsap_flux))
missingS2 <- table(is.na(s2$pdcsap_flux))
missingS3 <- table(is.na(s3$pdcsap_flux))

missingValues <- data.frame(
  Dataset = c("TIC 031381302", "TIC 129646813", "TIC 0131799991"),
  Available = c(missingS1["FALSE"], missingS2["FALSE"], missingS3["FALSE"]),
  Missing = c(missingS1["TRUE"], missingS2["TRUE"], missingS3["TRUE"])
)
kable(missingValues, caption = "Number of Missing Values in PDCSAP_FLUX by Star")

# Imputation using Kalman filter
s1$pdcsap_flux_imputed <- na_kalman(s1$pdcsap_flux, model = "StructTS")
s2$pdcsap_flux_imputed <- na_kalman(s2$pdcsap_flux, model = "StructTS")
s3$pdcsap_flux_imputed <- na_kalman(s3$pdcsap_flux, model = "StructTS")

# Plot original vs imputed flux for each star
pl1 <- ggplot(s1, aes(x = time)) +
  geom_line(aes(y = pdcsap_flux), color = "blue", alpha = 0.5) +
  geom_line(aes(y = pdcsap_flux_imputed), color = "red", linetype = "dashed", 
            alpha = 0.5) +
  labs(title = "Comparison of Original and Imputed Flux for TIC 031381302", 
       x = "Time", y = "PDCSAP_FLUX")

pl2 <- ggplot(s2, aes(x = time)) +
  geom_line(aes(y = pdcsap_flux), color = "blue", alpha = 0.5) +
  geom_line(aes(y = pdcsap_flux_imputed), color = "red", linetype = "dashed", 
            alpha = 0.5) +
  labs(title = "Comparison of Original and Imputed Flux for TIC 129646813", 
       x = "Time", y = "PDCSAP_FLUX")

pl3 <- ggplot(s3, aes(x = time)) +
  geom_line(aes(y = pdcsap_flux), color = "blue", alpha = 0.5) +
  geom_line(aes(y = pdcsap_flux_imputed), color = "red", linetype = "dashed", 
            alpha = 0.5) +
  labs(title = "Comparison of Original and Imputed Flux for TIC 0131799991", 
       x = "Time", y = "PDCSAP_FLUX")

# Distribution of PDCSAP Flux
data <- bind_rows(s1, s2, s3)
data %>%
  ggplot(aes(x = pdcsap_flux_imputed)) +
  geom_histogram(binwidth = 0.1, color = "blue", alpha = 0.7) +
  facet_wrap(~ star_id, scales = "free") +
  labs(title = "Histogram of PDCSAP Flux", x = "PDCSAP_FLUX", y = "Count")

# Time series of imputed flux for each star
p1 <- s1 %>%
  ggplot(aes(x = time, y = pdcsap_flux_imputed)) +
  geom_line() +
  labs(title = "Time Series of PDCSAP_FLUX for TIC 031381302", 
       x = "Time", y = "Flux")

p2 <- s2 %>%
  ggplot(aes(x = time, y = pdcsap_flux_imputed)) +
  geom_line() +
  labs(title = "Time Series of PDCSAP_FLUX for TIC 129646813", 
       x = "Time", y = "Flux")

p3 <- s3 %>%
  ggplot(aes(x = time, y = pdcsap_flux_imputed)) +
  geom_line() +
  labs(title = "Time Series of PDCSAP_FLUX for TIC 0131799991", 
       x = "Time", y = "Flux")

# Correlation of Flux Errors
c1 <- cor(s1$pdcsap_flux, s1$pdcsap_flux_err, use = "complete.obs")
c2 <- cor(s2$pdcsap_flux, s2$pdcsap_flux_err, use = "complete.obs")
c3 <- cor(s3$pdcsap_flux, s3$pdcsap_flux_err, use = "complete.obs")

corDf <- data.frame(
  Star = c("TIC 031381302", "TIC 129646813", "TIC 0131799991"),
  Correlation = c(c1, c2, c3)
)

kable(corDf, caption = "PDCSAP Flux and Flux Error Correlation by Star", 
      align = 'c', format = "markdown")

# Identifying Potential Flares
thresholdS1_1 <- mean(s1$pdcsap_flux_imputed, na.rm = TRUE) + 4 * 
  sd(s1$pdcsap_flux_imputed, na.rm = TRUE)
thresholdS1_2 <- mean(s1$pdcsap_flux_imputed, na.rm = TRUE) + 3 * 
  sd(s1$pdcsap_flux_imputed, na.rm = TRUE)

s1 <- s1 %>% mutate(
  Classification = if_else(pdcsap_flux_imputed < thresholdS1_2, 
                           "Normal", "Normal"),
  Classification = if_else(pdcsap_flux_imputed > thresholdS1_2, 
                           "Potential Flare", Classification),
  Classification = if_else(pdcsap_flux_imputed > thresholdS1_1, 
                           "Strong Potential Flare", Classification)
)

thresholdS2_1 <- mean(s2$pdcsap_flux_imputed, na.rm = TRUE) + 4 * 
  sd(s2$pdcsap_flux_imputed, na.rm = TRUE)
thresholdS2_2 <- mean(s2$pdcsap_flux_imputed, na.rm = TRUE) + 3 * 
  sd(s2$pdcsap_flux_imputed, na.rm = TRUE)

s2 <- s2 %>% mutate(
  Classification = if_else(pdcsap_flux_imputed < thresholdS2_2, 
                           "Normal", "Normal"),
  Classification = if_else(pdcsap_flux_imputed > thresholdS2_2, 
                           "Potential Flare", Classification),
  Classification = if_else(pdcsap_flux_imputed > thresholdS2_1, 
                           "Strong Potential Flare", Classification)
)

thresholdS3_1 <- mean(s3$pdcsap_flux_imputed, na.rm = TRUE) + 4 * 
  sd(s3$pdcsap_flux_imputed, na.rm = TRUE)
thresholdS3_2 <- mean(s3$pdcsap_flux_imputed, na.rm = TRUE) + 3 * 
  sd(s3$pdcsap_flux_imputed, na.rm = TRUE)

s3 <- s3 %>% mutate(
  Classification = if_else(pdcsap_flux_imputed < thresholdS3_2, 
                           "Normal", "Normal"),
  Classification = if_else(pdcsap_flux_imputed > thresholdS3_2, 
                           "Potential Flare", Classification),
  Classification = if_else(pdcsap_flux_imputed > thresholdS3_1, 
                           "Strong Potential Flare", Classification)
)

# Plot potential flares
q1 <- ggplot(s1, aes(x = time, y = pdcsap_flux_imputed, 
                     color = Classification)) +
  geom_point(alpha = 0.6) +
  labs(title = "Potential Flares for TIC 031381302", 
       x = "Time", y = "PDCSAP_FLUX") +
  scale_color_manual(values = c("Normal" = "black", 
                                "Strong Potential Flare" = "red", 
                                "Potential Flare" = "blue")) 

q2 <- ggplot(s2, aes(x = time, y = pdcsap_flux_imputed, 
                     color = Classification)) +
  geom_point(alpha = 0.6) +
  labs(title = "Potential Flares for TIC 129646813", 
       x = "Time", y = "PDCSAP_FLUX") +
  scale_color_manual(values = c("Normal" = "black", 
                                "Strong Potential Flare" = "red", 
                                "Potential Flare" = "blue")) 

q3 <- ggplot(s3, aes(x = time, y = pdcsap_flux_imputed, 
                     color = Classification)) +
  geom_point(alpha = 0.6) +
  labs(title = "Potential Flares for TIC 0131799991", 
       x = "Time", y = "PDCSAP_FLUX") +
  scale_color_manual(values = c("Normal" = "black", 
                                "Strong Potential Flare" = "red", 
                                "Potential Flare" = "blue")) 

# Clustering Brightness Events with DBSCAN
set.seed(2453)
dataCluster <- data %>%
  group_by(star_id) %>%
  mutate(scaled_flux = scale(pdcsap_flux_imputed)) %>%
  ungroup()

# Apply DBSCAN separately for each star
dataCluster <- dataCluster %>%
  group_by(star_id) %>%
  mutate(cluster = factor(dbscan(scaled_flux, 
                                 eps = 0.5,
                                 minPts = 5)$cluster)) %>%
  ungroup()

# Plot DBSCAN results
ggplot(dataCluster, aes(x = time, y = scaled_flux, color = cluster)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ star_id, scales = "free") +
  labs(title = "DBSCAN Clustered Brightness Events by Star", 
       x = "Time", y = "Scaled Flux") +
  scale_color_discrete(name = "Cluster")

# Cluster Summary
dataClusterSummary <- dataCluster %>% group_by(star_id, cluster) %>% 
  summarise(count = n())
colnames(dataClusterSummary) <- c("Star_ID", "Cluster", "Count")
kable(dataClusterSummary, caption = "Number of Points in Cluster by Star")

# ARMA to Analyze Residuals
pdcsap_flux_ts1 <- ts(s1$pdcsap_flux_imputed, frequency = 365)
arma_model1 <- Arima(pdcsap_flux_ts1, order = c(1, 0, 1))
residuals1 <- residuals(arma_model1)

pdcsap_flux_ts2 <- ts(s2$pdcsap_flux_imputed, frequency = 365)
arma_model2 <- Arima(pdcsap_flux_ts2, order = c(1, 0, 1))
residuals2 <- residuals(arma_model2)

pdcsap_flux_ts3 <- ts(s3$pdcsap_flux_imputed, frequency = 365)
arma_model3 <- Arima(pdcsap_flux_ts3, order = c(1, 0, 1))
residuals3 <- residuals(arma_model3)

acf(residuals2, main = "ACF of Residuals for TIC 129646813")
