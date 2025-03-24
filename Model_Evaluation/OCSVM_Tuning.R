# This script tests combinations of nu and gamma in the OC SVM model
# to optimize anomaly detection performance.

# Define possible values for nu and gamma
nu_values <- seq(0.05, 0.5, by = 0.05)  # Range of nu values
gamma_values <- c(0.01, 0.05, 0.1, 0.5)  # Range of gamma values

# Initialize a data frame to store results
results <- data.frame(nu = numeric(), gamma = numeric(), accuracy = numeric())

# Loop through each combination of nu and gamma
for (nu in nu_values) {
  for (gamma in gamma_values) {
    
    # Train the OCSVM model
    ocsvmModel <- svm(s1$pdcsap_flux_imputed, type = "one-classification", 
                      nu = nu, kernel = "radial", scale = TRUE, gamma = gamma, 
                      decision.values = TRUE)
    
    # Get the decision function scores
    decisionScores <- ocsvmModel$decision.values
    
    # Set threshold based on quantiles or another method
    threshold <- quantile(decisionScores, 0.90)
    #threshold <- -1000
    
    # Predict anomalies based on the threshold
    predicted_flare_ocsvm <- ifelse(decisionScores < threshold, 1, 0)
    
    # Calculate accuracy (or other metrics)
    accuracy <- mean(predicted_flare_ocsvm == simulatedS1$flare_label)
    
    # Store the results
    results <- rbind(results, data.frame(nu = nu, gamma = gamma, accuracy = accuracy))
  }
}

# View results to identify best parameters
best_params <- results[which.max(results$accuracy), ]
cat("Best parameters: nu =", best_params$nu, "gamma =", best_params$gamma, "\n")
cat("Best accuracy: ", best_params$accuracy, "\n")

#Best parameters: nu = 0.05 gamma = 0.01 
