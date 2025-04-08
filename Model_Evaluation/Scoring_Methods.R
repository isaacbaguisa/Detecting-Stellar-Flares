# This script scores detection methods (GPR, PPR, and OCSVM) 
# based on detection rate, sensitivity (true positive rate), 
# and predictive accuracy for detecting flares in the simulated light curves 
# of Star 1

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

# Baseline (Always predicting no flares)
simulatedS1$baseline_prediction <- 0
baselineAccuracy <- mean(simulatedS1$baseline_prediction == simulatedS1$flare_label)
cat("Baseline Accuracy: ", baselineAccuracy, "\n")

######################### GPR ######################################

# Scoring GPR on Star 1

# Predict using the trained GPR model on the simulated data
gprPred <- predict(gprModel, newdata = data.frame(time = simulatedS1$time))

# Calculate residuals 
gprResiduals <- simulatedS1$pdcsap_flux_imputed - gprPred

# Define an anomaly threshold (e.g., residuals greater than 3 sd from the mean)
thresholdGpr <- mean(gprResiduals, na.rm = TRUE) + 3 * sd(gprResiduals, 
                                                          na.rm = TRUE)

# Flag anomalies based on the residuals
simulatedS1$predicted_flare_gpr <- ifelse(gprResiduals > thresholdGpr, 1, 0)  

# Calculate TN, TP, FP, FN for GPR
predictedFlaresGpr <- sum(simulatedS1$predicted_flare_gpr==1)
truePositivesGpr <- sum(simulatedS1$predicted_flare_gpr == 1 & 
                          simulatedS1$flare_label == 1)
falseNegativesGpr <- sum(simulatedS1$predicted_flare_gpr == 0 & 
                           simulatedS1$flare_label == 1)
falsePositivesGpr <- sum(simulatedS1$predicted_flare_gpr == 1 & 
                           simulatedS1$flare_label == 0)
trueNegativesGpr <- sum(simulatedS1$predicted_flare_gpr == 0 & 
                          simulatedS1$flare_label == 0)

# Calculate FPR, Sensitivity and Accuracy for GPR
falsePostiveRateGpr <- falsePositivesGpr / (falsePositivesGpr + trueNegativesGpr)
sensitivityGpr <- truePositivesGpr / (truePositivesGpr + falseNegativesGpr)
accuracyGpr <- mean(simulatedS1$predicted_flare_gpr == simulatedS1$flare_label)

# Print the results
cat("GPR Flares Predicted:", predictedFlaresGpr, "\n")
cat("GPR Sensitivity (True Positive Rate): ", sensitivityGpr, "\n")
cat("GPR False Positive Rate: ", falsePostiveRateGpr, "\n")
cat("GPR Predictive Accuracy: ", accuracyGpr, "\n")
#which(simulatedS1$flare_label==1 & simulatedS1$predicted_flare_gpr==0)

# Scoring GPR Trained on Star 1 on Star 3
# Predict using the trained GPR model on the simulated data for Star 3
gprPred <- predict(gprModel, newdata = data.frame(time = simulatedS3$time))

# Calculate residuals 
gprResiduals <- simulatedS3$pdcsap_flux_imputed - gprPred

# Define an anomaly threshold (residuals greater than 3 sd from the mean)
thresholdGpr <- mean(gprResiduals, na.rm = TRUE) + 3 * sd(gprResiduals, 
                                                          na.rm = TRUE)

# Flag anomalies based on the residuals
simulatedS3$predicted_flare_gpr <- ifelse(gprResiduals > thresholdGpr, 1, 0)  

# Calculate TN, TP, FP, FN for GPR
predictedFlaresGpr <- sum(simulatedS3$predicted_flare_gpr==1)
truePositivesGpr <- sum(simulatedS3$predicted_flare_gpr == 1 & 
                          simulatedS3$flare_label == 1)
falseNegativesGpr <- sum(simulatedS3$predicted_flare_gpr == 0 & 
                           simulatedS3$flare_label == 1)
falsePositivesGpr <- sum(simulatedS3$predicted_flare_gpr == 1 & 
                           simulatedS3$flare_label == 0)
trueNegativesGpr <- sum(simulatedS3$predicted_flare_gpr == 0 & 
                         simulatedS3$flare_label == 0)

# Calculate FPR, Sensitivity and Accuracy for GPR
falsePostiveRateGpr <- falsePositivesGpr / (falsePositivesGpr + trueNegativesGpr)
sensitivityGpr <- truePositivesGpr / (truePositivesGpr + falseNegativesGpr)
accuracyGpr <- mean(simulatedS3$predicted_flare_gpr == simulatedS3$flare_label)

# Print the results
cat("GPR Flares Predicted:", predictedFlaresGpr, "\n")
cat("GPR Sensitivity (True Positive Rate): ", sensitivityGpr, "\n")
cat("GPR False Positive Rate: ", falsePostiveRateGpr, "\n")
cat("GPR Predictive Accuracy: ", accuracyGpr, "\n")

#which(simulatedS3$flare_label==1 & simulatedS3$predicted_flare_gpr==0)

ggplot(simulatedS3, aes(x = time, y = pdcsap_flux_imputed)) + 
  geom_line(color = "green") +
  geom_point(data = subset(simulatedS3, predicted_flare_gpr == 1),
             aes(x = time, y = pdcsap_flux_imputed),
             color = "red", size = 2) +
  labs(title = "Flare Detection using GPR on TIC 0131799991 Simulation",
       x = "Time", y = "PDCSAP_FLUX") + 
         theme_minimal()

######################### PPR ######################################

# Predict using the trained PPR model on the simulated data for Star 3
pprPred <- predict(pprModel, newdata = data.frame(time = simulatedS3$time))

# Calculate residuals (difference between observed and predicted flux)
pprResiduals <- simulatedS3$pdcsap_flux_imputed - pprPred

# Define an anomaly threshold (residuals greater than 3 sd from the mean)
thresholdPpr <- mean(pprResiduals, na.rm = TRUE) + 3 * sd(pprResiduals, 
                                                          na.rm = TRUE)

# Flag anomalies based on the residuals
simulatedS3$predicted_flare_ppr <- ifelse(pprResiduals > thresholdPpr, 1, 0)  

# Calculate TN, TP, FP, FN for PPR
truePositivesPpr <- sum(simulatedS3$predicted_flare_ppr == 1 & 
                          simulatedS3$flare_label == 1)
falseNegativesPpr <- sum(simulatedS3$predicted_flare_ppr == 0 & 
                           simulatedS3$flare_label == 1)
falsePositivesPpr <- sum(simulatedS3$predicted_flare_ppr == 1 & 
                           simulatedS3$flare_label == 0)

sensitivityPpr <- truePositivesPpr / (truePositivesPpr + falseNegativesPpr)

# Calculate Accuracy for PPR
accuracyPpr <- mean(simulatedS3$predicted_flare_ppr == simulatedS3$flare_label)

# Print the results
cat("PPR Sensitivity (True Positive Rate): ", sensitivityPpr, "\n")
cat("PPR Predictive Accuracy: ", accuracyPpr, "\n")

######################### OC SVM ######################################

# Scoring SVM on Star 1
ocsvmModel <- svm(s1$pdcsap_flux_imputed, type = "one-classification", nu = 0.01,
                  kernel = "radial", scale = TRUE, decision.values=TRUE,
                  gamma = 0.01)

# Predict using the trained OCSVM model on the simulated data for Star 1
ocsvmPred <- predict(ocsvmModel, 
                     newdata = data.frame(simulatedS1$pdcsap_flux_imputed))

# Define threshold
simulatedS1$predicted_flare_ocsvm <- ifelse(ocsvmPred, 0, 1)
threshold <- quantile(simulatedS1$pdcsap_flux_imputed, 0.95)  
simulatedS1$predicted_flare_ocsvm <- ifelse(simulatedS1$predicted_flare_ocsvm == 1 & 
                                        simulatedS1$pdcsap_flux_imputed > threshold, 1, 0)

# Calculate TN, TP, FN, FP OCSVM
flaresPredictedOcsvm <- sum(simulatedS1$predicted_flare_ocsvm == 1)
truePositivesOcsvm <- sum(simulatedS1$predicted_flare_ocsvm == 1 & 
                            simulatedS1$flare_label == 1)
falseNegativesOcsvm <- sum(simulatedS1$predicted_flare_ocsvm == 0 & 
                             simulatedS1$flare_label == 1)
falsePositivesOcsvm <- sum(simulatedS1$predicted_flare_ocsvm == 1 & 
                             simulatedS1$flare_label == 0)
trueNegativesOcsvm <- sum(simulatedS1$predicted_flare_ocsvm == 0 & 
                          simulatedS1$flare_label == 0)

# Calculate FPR, Sensitivity and Accuracy for OCSVM
falsePostiveRateOcsvm <- falsePositivesOcsvm / (falsePositivesOcsvm + trueNegativesOcsvm)
sensitivityOcsvm <- truePositivesOcsvm / (truePositivesOcsvm + falseNegativesOcsvm)
accuracyOcsvm <- mean(simulatedS1$predicted_flare_ocsvm == simulatedS1$flare_label)

# Print the results
cat("OCSVM Flares Predicted:", flaresPredictedOcsvm, "\n")
cat("OCSVM Sensitivity (True Positive Rate): ", sensitivityOcsvm, "\n")
cat("OCSVM False Positive Rate: ", falsePostiveRateOcsvm, "\n")
cat("OCSVM Predictive Accuracy: ", accuracyOcsvm, "\n")

# Visualize flux over time with anomalies highlighted
ggplot(simulatedS1, aes(x = time, y = pdcsap_flux_imputed)) + 
  geom_line(color = "black") +
  geom_point(data = subset(simulatedS1, predicted_flare_ocsvm == 1),
             aes(x = time, y = pdcsap_flux_imputed),
             color = "red", size = 2) +
  labs(title = "Anomaly Detection using One-Class SVM",
       x = "Time", y = "Flux",
       subtitle = "Red points indicate detected anomalies") +
  theme_minimal()

#which(simulatedS1$flare_label==1 & simulatedS1$predicted_flare_ocsvm==0)
#simulatedS1[12706,]

# Scoring SVM on Star 3
ocsvmModel <- svm(s3$pdcsap_flux_imputed, type = "one-classification", nu = 0.01,
                 kernel = "radial", scale = TRUE, 
                  gamma = 0.01, decision.values = TRUE)

# Predict using the trained OCSVM model on the simulated data for Star 3
ocsvmPred <- predict(ocsvmModel, 
                     newdata = data.frame(simulatedS3$pdcsap_flux_imputed))

simulatedS3$predicted_flare_ocsvm <- ifelse(ocsvmPred, 0, 1)

threshold <- quantile(simulatedS3$pdcsap_flux_imputed, 0.95)  

simulatedS3$predicted_flare_ocsvm <- ifelse(simulatedS3$predicted_flare_ocsvm == 1 & 
                                              simulatedS3$pdcsap_flux_imputed > threshold, 1, 0)

# Calculate TN, TP, FP, FN for OCSVM
flaresPredictedOcsvm <- sum(simulatedS3$predicted_flare_ocsvm == 1)
truePositivesOcsvm <- sum(simulatedS3$predicted_flare_ocsvm == 1 & 
                            simulatedS3$flare_label == 1)
falseNegativesOcsvm <- sum(simulatedS3$predicted_flare_ocsvm == 0 & 
                             simulatedS3$flare_label == 1)
falsePositivesOcsvm <- sum(simulatedS3$predicted_flare_ocsvm == 1 & 
                             simulatedS3$flare_label == 0)
trueNegativesOcsvm <- sum(simulatedS3$predicted_flare_ocsvm == 0 & 
                            simulatedS3$flare_label == 0)

# Calculate FPR, Sensitivity and Accuracy for GPR
falsePostiveRateOcsvm <- falsePositivesOcsvm / (falsePositivesOcsvm + trueNegativesOcsvm)
sensitivityOcsvm <- truePositivesOcsvm / (truePositivesOcsvm + falseNegativesOcsvm)
accuracyOcsvm <- mean(simulatedS3$predicted_flare_ocsvm == simulatedS3$flare_label)

# Print the results
cat("OCSVM Flares Predicted:", flaresPredictedOcsvm, "\n")
cat("OCSVM Sensitivity (True Positive Rate): ", sensitivityOcsvm, "\n")
cat("OCSVM False Positive Rate: ", falsePostiveRateOcsvm, "\n")
cat("OCSVM Predictive Accuracy: ", accuracyOcsvm, "\n")

# which(simulatedS3$flare_label==1 & simulatedS3$predicted_flare_ocsvm==0)
#simulatedS3[2688,]

# Visualize flux over time with anomalies highlighted
ggplot(simulatedS3, aes(x = time, y = pdcsap_flux_imputed)) + 
  geom_line(color = "green") +
  geom_point(data = subset(simulatedS3, predicted_flare_ocsvm == 1),
             aes(x = time, y = pdcsap_flux_imputed),
             color = "red", size = 2) +
  labs(title = "Flare Detection using One-Class SVM on TIC 0131799991 Simulation",
       x = "Time", y = "PDCSAP_FLUX") +
  theme_minimal()

