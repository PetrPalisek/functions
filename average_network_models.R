average_network_models <- function(data, LV_names, MV_names, draws = 3, 
                                   estimator = "EBICglasso", corMethod = "spearman", tuning = .15) {
  
  library(bootnet)
  library(dplyr)
  
  # List to store weight matrices from each valid draw
  weight_matrices <- list()
  
  # Loop over plausible value draws
  for (i in 1:draws) {
    
    # Select the relevant columns for this draw
    draw_cols <- data %>% select(c(MV_names, 
                                   ends_with(paste0("_draw", i)) & 
                                     starts_with(LV_names)))
    
    # Compute the correlation matrix
    correlation_matrix <- cor(draw_cols, method = corMethod, use = "complete")
    
    # Check if the matrix is positive definite by looking at its eigenvalues
    eigenvalues <- eigen(correlation_matrix)$values
    
    if (any(eigenvalues <= 0)) {
      # Skip this draw if the matrix is not positive definite
      next
    } else {
      # Fit the network model using bootnet
      net <- bootnet::estimateNetwork(draw_cols, 
                                      default = estimator, 
                                      corMethod = corMethod, 
                                      tuning = tuning,
                                      missing = "listwise")
      
      # Store the weight matrix
      weight_matrices[[length(weight_matrices) + 1]] <- as.matrix(net$graph)
    }
  }
  
  # Check if any valid weight matrices were found
  if (length(weight_matrices) == 0) {
    stop("No valid positive definite matrices found among the draws.")
  }
  
  # Stack the weight matrices into a 3D array for easy computation
  weight_array <- array(unlist(weight_matrices), 
                        dim = c(nrow(weight_matrices[[1]]), 
                                ncol(weight_matrices[[1]]), 
                                length(weight_matrices)))
  
  # Compute the average weight matrix
  average_matrix <- apply(weight_array, c(1, 2), mean) %>% data.frame()
  names(average_matrix) <- names(draw_cols)
  
  # Compute variability statistics
  sd_matrix <- apply(weight_array, c(1, 2), sd)   # Standard deviation
  range_matrix <- apply(weight_array, c(1, 2), function(x) diff(range(x))) # Range
  cv_matrix <- sd_matrix / average_matrix         # Coefficient of variation
  
  # Print statistics about variability
  cat("\nStatistics about Variability in Averaged Parameters:\n")
  cat("----------------------------------------------------\n")
  cat("Standard Deviation (Average across parameters):", mean(sd_matrix), "\n")
  cat("Coefficient of Variation (Average across parameters):", mean(cv_matrix, na.rm = TRUE), "\n")
  cat("Range (Average across parameters):", mean(range_matrix), "\n")
  
  # Optionally: Return both the averaged matrix and the variability stats
  return(list(average_matrix = average_matrix, 
              sd_matrix = sd_matrix, 
              cv_matrix = cv_matrix, 
              range_matrix = range_matrix))
}
