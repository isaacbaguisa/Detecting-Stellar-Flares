# This script scores detection methods (GPR, PPR, and OCSVM) 
# based on detection rate, sensitivity (true positive rate), 
# and predictive accuracy for detecting flares in the simulated light curves 
# of Star 3

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

#########################GPR######################################


# Scoring GPR on Star 1

# Predict using the trained GPR model on the simulated data for Star 3
gprPred <- predict(gprModel, newdata = data.frame(time = simulatedS1$time))

# Calculate residuals (difference between observed and predicted flux)
gprResiduals <- simulatedS1$pdcsap_flux_imputed - gprPred

# Define an anomaly threshold (e.g., residuals greater than 3 standard deviations from the mean)
thresholdGpr <- mean(gprResiduals, na.rm = TRUE) + 3 * sd(gprResiduals, 
                                                          na.rm = TRUE)

# Flag anomalies based on the residuals
simulatedS1$predicted_flare_gpr <- ifelse(gprResiduals > thresholdGpr, 1, 0)  
# 1 for flare, 0 for normal

# Calculate Detection Rate for GPR
detectionRateGpr <- mean(simulatedS1$predicted_flare_gpr == 1 & 
                           simulatedS1$flare_label == 1)

# Calculate Sensitivity (True Positive Rate) for GPR
truePositivesGpr <- sum(simulatedS1$predicted_flare_gpr == 1 & 
                          simulatedS1$flare_label == 1)
falseNegativesGpr <- sum(simulatedS1$predicted_flare_gpr == 0 & 
                           simulatedS1$flare_label == 1)
falsePositivesGpr <- sum(simulatedS1$predicted_flare_gpr == 1 & 
                           simulatedS1$flare_label == 0)

sensitivityGpr <- truePositivesGpr / (truePositivesGpr + falseNegativesGpr)

# Calculate Predictive Accuracy for GPR
accuracyGpr <- mean(simulatedS1$predicted_flare_gpr == simulatedS1$flare_label)

# Print the results
cat("GPR Detection Rate: ", detectionRateGpr, "\n")
cat("GPR Sensitivity (True Positive Rate): ", sensitivityGpr, "\n")
cat("GPR Predictive Accuracy: ", accuracyGpr, "\n")


sum(simulatedS1$flare_label==1)
sum(simulatedS1$predicted_flare_gpr==1)

which(simulatedS1$flare_label==1)
which(simulatedS1$predicted_flare_gpr==1)

which(simulatedS1$flare_label & simulatedS1$predicted_flare_gpr == 1)


# Scoring GPR Trained on Star 1, on Star 3
# Predict using the trained GPR model on the simulated data for Star 3
gprPred <- predict(gprModel, newdata = data.frame(time = simulatedS3$time))

# Calculate residuals (difference between observed and predicted flux)
gprResiduals <- simulatedS3$pdcsap_flux_imputed - gprPred

# Define an anomaly threshold (e.g., residuals greater than 3 standard deviations from the mean)
thresholdGpr <- mean(gprResiduals, na.rm = TRUE) + 3 * sd(gprResiduals, 
                                                          na.rm = TRUE)

# Flag anomalies based on the residuals
simulatedS3$predicted_flare_gpr <- ifelse(gprResiduals > thresholdGpr, 1, 0)  
# 1 for flare, 0 for normal

# Calculate Detection Rate for GPR
detectionRateGpr <- mean(simulatedS3$predicted_flare_gpr == 1 & 
                           simulatedS3$flare_label == 1)

# Calculate Sensitivity (True Positive Rate) for GPR
truePositivesGpr <- sum(simulatedS3$predicted_flare_gpr == 1 & 
                          simulatedS3$flare_label == 1)
falseNegativesGpr <- sum(simulatedS3$predicted_flare_gpr == 0 & 
                           simulatedS3$flare_label == 1)
falsePositivesGpr <- sum(simulatedS3$predicted_flare_gpr == 1 & 
                           simulatedS3$flare_label == 0)

sensitivityGpr <- truePositivesGpr / (truePositivesGpr + falseNegativesGpr)

# Calculate Predictive Accuracy for GPR
accuracyGpr <- mean(simulatedS3$predicted_flare_gpr == simulatedS3$flare_label)

# Print the results
cat("GPR Detection Rate: ", detectionRateGpr, "\n")
cat("GPR Sensitivity (True Positive Rate): ", sensitivityGpr, "\n")
cat("GPR Predictive Accuracy: ", accuracyGpr, "\n")


sum(simulatedS3$flare_label==1)
sum(simulatedS3$predicted_flare_gpr==1)

which(simulatedS3$flare_label==1)
which(simulatedS3$predicted_flare_gpr==1)

which(simulatedS3$flare_label & simulatedS3$predicted_flare_gpr == 1)

#########################PPR######################################


# Predict using the trained PPR model on the simulated data for Star 3
pprPred <- predict(pprModel, newdata = data.frame(time = simulatedS3$time))

# Calculate residuals (difference between observed and predicted flux)
pprResiduals <- simulatedS3$pdcsap_flux_imputed - pprPred

# Define an anomaly threshold (residuals greater than 3 sd from the mean)
thresholdPpr <- mean(pprResiduals, na.rm = TRUE) + 3 * sd(pprResiduals, 
                                                          na.rm = TRUE)

# Flag anomalies based on the residuals
simulatedS3$predicted_flare_ppr <- ifelse(pprResiduals > thresholdPpr, 1, 0)  
# 1 for flare, 0 for normal

# Calculate Detection Rate for PPR
detectionRatePpr <- mean(simulatedS3$predicted_flare_ppr == 1 & 
                           simulatedS3$flare_label == 1)

# Calculate Sensitivity (True Positive Rate) for PPR
truePositivesPpr <- sum(simulatedS3$predicted_flare_ppr == 1 & 
                          simulatedS3$flare_label == 1)
falseNegativesPpr <- sum(simulatedS3$predicted_flare_ppr == 0 & 
                           simulatedS3$flare_label == 1)
falsePositivesPpr <- sum(simulatedS3$predicted_flare_ppr == 1 & 
                           simulatedS3$flare_label == 0)

sensitivityPpr <- truePositivesPpr / (truePositivesPpr + falseNegativesPpr)

# Calculate Predictive Accuracy for PPR
accuracyPpr <- mean(simulatedS3$predicted_flare_ppr == simulatedS3$flare_label)

# Print the results
cat("PPR Detection Rate: ", detectionRatePpr, "\n")
cat("PPR Sensitivity (True Positive Rate): ", sensitivityPpr, "\n")
cat("PPR Predictive Accuracy: ", accuracyPpr, "\n")

sum(simulatedS3$flare_label==1)
sum(simulatedS3$predicted_flare_ppr==1)

which(simulatedS3$flare_label==1)
which(simulatedS3$predicted_flare_ppr==1)

which(simulatedS3$flare_label & simulatedS3$predicted_flare_ppr == 1)


#########################OC SVM######################################

# Scoring SVM on Star 1
# Fit the One-Class SVM model with decision function enabled
ocsvmModel <- svm(s1$pdcsap_flux_imputed, type = "one-classification", nu = 0.01,
                  kernel = "radial", scale = TRUE, 
                  gamma = 0.05, decision.values = TRUE)

# Predict using the trained OCSVM model on the simulated data for Star 3
# The OCSVM model returns -1 for outliers (flares) and 1 for normal points
ocsvmPred <- predict(ocsvmModel, 
                     newdata = data.frame(simulatedS1$pdcsap_flux_imputed))

decisionScores <- ocsvmModel$decision.values
threshold <- quantile(decisionScores, 0.95) 
#threshold <- -20

# Flag anomalies based on the threshold
simulatedS1$predicted_flare_ocsvm <- ifelse(decisionScores < threshold, 1, 0)  # 1 for flare, 0 for normal

# Print number of anomalies detected
cat("Number of anomalies detected: ", sum(simulatedS1$predicted_flare_ocsvm == 1), "\n")

# Create a column for the predicted flare labels based on OCSVM results
#simulatedS1$predicted_flare_ocsvm <- ifelse(ocsvmPred == -1, 1, 0)  
# 1 for flare, 0 for normal

# Calculate Detection Rate (the fraction of actual flares detected as anomalies)
detectionRateOcsvm <- mean(simulatedS1$predicted_flare_ocsvm == 1 & 
                             simulatedS1$flare_label == 1)

# Calculate Sensitivity (True Positive Rate) for OCSVM
truePositivesOcsvm <- sum(simulatedS1$predicted_flare_ocsvm == 1 & 
                            simulatedS1$flare_label == 1)
falseNegativesOcsvm <- sum(simulatedS1$predicted_flare_ocsvm == 0 & 
                             simulatedS1$flare_label == 1)
falsePositivesOcsvm <- sum(simulatedS1$predicted_flare_ocsvm == 1 & 
                             simulatedS1$flare_label == 0)
sensitivityOcsvm <- truePositivesOcsvm / (truePositivesOcsvm + falseNegativesOcsvm)

# Calculate Predictive Accuracy for OCSVM
accuracyOcsvm <- mean(simulatedS1$predicted_flare_ocsvm == simulatedS1$flare_label)

# Print the results
cat("OCSVM Detection Rate: ", detectionRateOcsvm, "\n")
cat("OCSVM Sensitivity (True Positive Rate): ", sensitivityOcsvm, "\n")
cat("OCSVM Predictive Accuracy: ", accuracyOcsvm, "\n")


sum(simulatedS1$flare_label==1)
sum(simulatedS1$predicted_flare_ocsvm==1)

which(simulatedS1$flare_label==1)
which(simulatedS1$predicted_flare_ocsvm==1)

which(simulatedS1$flare_label & simulatedS1$predicted_flare_ocsvm == 1)

