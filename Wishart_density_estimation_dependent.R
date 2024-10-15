####################################################################################################
## Density estimation on the space of positive definite matrices using Wishart asymmetric kernels ##
####################################################################################################

## Written by Frederic Ouimet (October 2024)

require("CholWishart")        # for the multivariate gamma function
require("cubature")           # for integrals
require("doFuture")           # for parallel execution with foreach
require("expm")               # for matrix logarithms and exponentials
require("fs")                 # for filesystem operations (dependency)
require("future.batchtools")  # for batchtools integration with future
require("ggplot2")            # for plotting
require("LaplacesDemon")      # for the Wishart and inverse Wishart distributions
require("MASS")               # for the multivariate normal distribution
require("matrixcalc")         # for functions like is.positive.definite
require("matrixsampling")     # for the matrix-type-II Beta distribution
require("optimx")             # for improved bandwidth optimization
require("parallel")           # for parallelization of calculations
require("tidyverse")          # for data manipulation and visualization
require("writexl")            # to write output to Excel files

#####################
## Parallelization ##
#####################

# Define the list of libraries to load on each cluster node

libraries_to_load <- c(
  "CholWishart",
  "cubature",
  "doFuture",
  "expm",
  "fs",
  "future.batchtools",
  "ggplot2",
  "LaplacesDemon",
  "MASS",
  "matrixcalc",
  "matrixsampling",
  "optimx",
  "parallel",
  "tidyverse",
  "writexl"
)

# Define the list of variables/functions to export to the worker nodes

vars_to_export <- c(
  "BB",
  "b_LG_test",
  "b_opt_MC",
  "b_opt_MC_grid",
  "b_test",
  "b_WK_test",
  "CholWishart",
  "compute_ISE",
  "construct_X",
  "cores_per_node",
  "cubature",
  "d",
  "delta",
  "dmatrixbeta_typeII",
  "expm",
  "f",
  "G",
  "generate_VAR1_collection",
  "generate_WAR1",
  "hat_f",
  "II",
  "II_test",
  "integrand",
  "ISE",
  "ISE_MC",
  "ISE_value",
  "JJ",
  "JJ_test",
  "K",
  "K_test",
  "LaplacesDemon",
  "LG",
  "lmvgamma",
  "local_raw_results",
  "logm",
  "LSCV",
  "LSCV_MC",
  "M_list",
  "matrixcalc",
  "matrixsampling",
  "method",
  "MM",
  "NN",
  "n_LG_test",
  "n_test",
  "n_WK_test",
  "optimx",
  "path",
  "raw_results",
  "resources_list",
  "rotation_matrix",
  "RR",
  "S1",
  "S2",
  "S3",
  "S4",
  "setup_parallel_cluster",
  "Sigma_inf_list",
  "Sigma_list",
  "solve_riccati",
  "symmetrize",
  "test_estimator_integral",
  "tol1",
  "tol2",
  "vars_to_export",
  "writexl",
  "XX",
  "XX_data",
  "XX_sample"
)

# Sets up a parallel cluster, loads necessary libraries, and exports required variables globally

setup_parallel_cluster <- function() {
  num_cores <<- detectCores() - 1
  cl <<- makeCluster(num_cores)
  
  # Export the list of libraries to the worker nodes
  clusterExport(cl, varlist = "libraries_to_load")
  
  # Load necessary libraries on each cluster node
  invisible(clusterEvalQ(cl, {
    lapply(libraries_to_load, library, character.only = TRUE)
  }))
  
  # Export all necessary objects, functions, and parameters to the worker nodes
  clusterExport(cl, varlist = vars_to_export)
  
  return(cl) # Return the cluster object
}

# Initialize all variables in the list as NULL except vars_to_export and setup_parallel_cluster

invisible(
  lapply(
    vars_to_export[!(vars_to_export %in% c("vars_to_export", "setup_parallel_cluster"))],
    function(x) assign(x, NULL, envir = .GlobalEnv)
  )
)

##################
## Set the path ##
##################

# path <- file.path(
#   "C://Users//fred1//Desktop//Github_WishartDensityEstimation",
#   fsep = .Platform$file.sep
# )
path <- getwd()
# setwd(path)

################
## Parameters ##
################

d <- 2 # width of the square matrices
delta <- 0.1 # lower bound on the eigenvalues of SPD matrices in LSCV
K <- 4 # fixed degree-of-freedom parameter

h <- function(n) { n ^ 0.25 } # for h-block cross-validation

MM <- list("WK", "LG") # list of density estimation methods
NN <- c(50) # sample sizes
II <- 1:1 # stationary density function indices (related to M)
JJ <- 1:1 # stationary density function indices (related to Sigma)
RR <- 1:1 # replication indices

cores_per_node <- 63 # number of cores for each node in the super-computer

tol1 <- 1e-1
tol2 <- 1e-1

##############################
## Parallelization on nodes ##
##############################

resources_list <- list(
  cpus_per_task = cores_per_node,
  mem = "240G",
  walltime = "28:00:00",
  nodes = 1
  # Omit 'partition' to let SLURM choose
)

#######################
## Special functions ##
#######################

# Symmetrize a matrix

symmetrize <- function(X) {
  return((X + t(X)) / 2)
}

# Generate a 2x2 rotation matrix

rotation_matrix <- function(theta) {
  return(matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2))
}

# Construct X = rotation_matrix(theta) . diag(lambda1, lambda2) . rotation_matrix(theta)^T

construct_X <- function(theta, lambda1, lambda2) {
  R_theta <- rotation_matrix(theta)
  diag_lambda <- diag(c(lambda1, lambda2))
  X <- R_theta %*% diag_lambda %*% t(R_theta)
  return((X + t(X)) / 2)  # Ensure the matrix is symmetric
}

# Bandwidth sequence based on the method

BB <- function(method) {
  if (method == "WK") {
    return(seq(0.01, 0.45, length.out = cores_per_node))  # BB_WK sequence
  } else if (method == "LG") {
    return(seq(0.02, 2, length.out = cores_per_node))  # BB_LG sequence
  } else {
    stop("Invalid method. Should be 'WK' or 'LG'.")
  }
}

##################
## Observations ##
##################

# Autoregressive coefficient matrices
M1 <- matrix(c(0.9, 0, 1, 0), nrow = d)
M2 <- matrix(c(0.3, -0.3, -0.3, 0.3), nrow = d)
M3 <- matrix(c(0.5, 0, 0, 0.5), nrow = d)

# Covariance matrices
Sigma1 <- matrix(c(1, 0, 0, 1), nrow = d)
Sigma2 <- matrix(c(1, 0.5, 0.5, 1), nrow = d)
Sigma3 <- matrix(c(1, 0.9, 0.9, 1), nrow = d)

# List of M and Sigma matrices
M_list <- list(M1, M2, M3)
Sigma_list <- list(Sigma1, Sigma2, Sigma3)

# Function to generate K independent VAR(1) processes
generate_VAR1_collection <- function(i, j, K, n) {
  M <- M_list[[i]]
  Sigma <- Sigma_list[[j]]
  
  # Initialize list to store X[t], each element is a K x d matrix
  X_list <- vector("list", n)
  
  # Generate the first time step X[1] as a K x d matrix
  X_t <- matrix(0, nrow = K, ncol = d)
  for (k in 1:K) {
    X_t[k, ] <- mvrnorm(1, mu = rep(0, d), Sigma = Sigma)
  }
  X_list[[1]] <- X_t
  
  # Generate subsequent time steps using the VAR(1) structure
  for (t in 2:n) {
    X_t <- matrix(0, nrow = K, ncol = d)
    for (k in 1:K) {
      X_t[k, ] <- M %*% X_list[[t - 1]][k, ] + mvrnorm(1, mu = rep(0, d), Sigma = Sigma)
    }
    X_list[[t]] <- X_t
  }
  
  return(X_list)
}

# Function to generate the WAR(1) process as a sum of outer products of VAR(1) processes
generate_WAR1 <- function(i, j, K, n) {
  # Initialize list to store Y[t], each element is a d x d matrix
  Y_list <- vector("list", n)
  
  # Generate the VAR(1) collection as a list of matrices
  X_list <- generate_VAR1_collection(i, j, K, n)
  
  # Compute Y[t] as the sum of outer products of X[k, t, ] over k from 1 to K
  for (t in 1:n) {
    Y_t <- matrix(0, nrow = d, ncol = d)
    for (k in 1:K) {
      X_kt <- X_list[[t]][k, ]
      Y_t <- Y_t + X_kt %*% t(X_kt)
    }
    Y_list[[t]] <- Y_t
  }
  
  return(Y_list)
}

# Example of generating a WAR(1) process
# i is the index of M (1, 2, or 3)
# j is the index of Sigma (1, 2, or 3)
# K is the number of independent VAR(1) processes
# n is the number of time components

XX <- function(i, j, K, n) {
  res <- generate_WAR1(i, j, K, n)
  return(res)
}

# # Example usage:
# XX(i = 1, j = 1, K = 5, n = 10)

##################################
## Stationary density functions ##
##################################

# Function to solve the Riccati equation Sigma_inf = M Sigma_inf M^T + Sigma
solve_riccati <- function(M, Sigma, tol = 1e-8, max_iter = 1000) {
  Sigma_inf <- Sigma  # Initial guess for Sigma_inf
  for (i in 1:max_iter) {
    Sigma_new <- M %*% Sigma_inf %*% t(M) + Sigma
    if (norm(Sigma_new - Sigma_inf, type = "F") < tol) break
    Sigma_inf <- Sigma_new
  }
  # Symmetrization step to ensure the result is symmetric
  Sigma_inf <- (Sigma_inf + t(Sigma_inf)) / 2
  return(Sigma_inf)
}

# Solve the Riccati equation for each combination of M_i and Sigma_j
Sigma_inf_list <- list()

for (i in II) {
  for (j in JJ) {
    # Store the result with key as just "ij"
    Sigma_inf_list[[paste0(i, j)]] <- solve_riccati(M_list[[i]], Sigma_list[[j]])
  }
}

# # Test all 9 combinations of i from 1 to 3, and j from 1 to 3
# for (i in II) {
#   for (j in JJ) {
#     # Get the Riccati solution and the matrices M_i and Sigma_j
#     Sigma_inf_ij <- Sigma_inf_list[[paste0(i, j)]]
#     M_i <- M_list[[i]]
#     Sigma_j <- Sigma_list[[j]]
# 
#     # Test to verify if the Riccati equation is solved for this combination
#     left_side <- Sigma_inf_ij  # Left side: Sigma_inf
#     right_side <- M_i %*% Sigma_inf_ij %*% t(M_i) + Sigma_j  # Right side: M Sigma_inf M^T + Sigma
# 
#     # Check if the two sides are approximately equal
#     test_result <- all.equal(left_side, right_side, tolerance = 1e-8)
# 
#     # Print the result based on whether the equation is satisfied
#     if (isTRUE(test_result)) {
#       print(paste0("Riccati equation satisfied for i = ", i, ", j = ", j))
#     } else {
#       print(paste0("Riccati equation NOT satisfied for i = ", i, ", j = ", j))
#       print("Left side (Sigma_inf):")
#       print(left_side)
#       print("Right side (M Sigma_inf M^T + Sigma):")
#       print(right_side)
#     }
#   }
# }

# Stationary densities
f <- function(i, j, K, X) { # i and j define M_i and Sigma_j, K is degrees of freedom, X is an SPD matrix
  # Retrieve Sigma_inf from the Riccati solutions
  Sigma_inf <- Sigma_inf_list[[paste0(i, j)]]
  
  return(LaplacesDemon::dwishart(X, K, Sigma_inf))
}

# # Tests if the stationary density f(i, j, 2, X) integrates to 1 over all i in II_test, j in JJ_test
# 
# # Set up the parallel environment
# cl <- setup_parallel_cluster()
# 
# # Set the range of indices for i and j
# II_test <- 1:2
# JJ_test <- 1:3
# 
# # Create a list of pairs in the same order as expand.grid(II_test, JJ_test)
# combinations <- lapply(1:nrow(expand.grid(II_test, JJ_test)), function(k) as.numeric(expand.grid(II_test, JJ_test)[k, ]))
# 
# # Start the timer
# start_time <- Sys.time()
# 
# # Run the integral computation in parallel for all combinations of i and j
# results <- parLapply(cl, combinations, function(params) {
#   i <- params[1]
#   j <- params[2]
#   K <- 4  # Fixed K value
#   
#   cat("Calculating integral for i =", i, ", j =", j, "\n")
#   
#   # Integrand function for the triple integral
#   integrand <- function(vars) {
#     density_value <- f(i, j, K, construct_X(vars[1], vars[2], vars[3]))
#     jacobian_value <- abs(vars[2] - vars[3]) / 4
#     
#     return(density_value * jacobian_value)
#   }
#   
#   # Define the limits of integration
#   lower_limit <- c(0, 0, 0)
#   upper_limit <- c(2 * pi, Inf, Inf)
#   
#   # Perform the numerical integration using the cubature package
#   result <- adaptIntegrate(integrand, lowerLimit = lower_limit, upperLimit = upper_limit, tol = tol1)
#   
#   cat("Result for i =", i, ", j =", j, ":", result$integral, "\n")
#   return(result$integral)
# })
# 
# # Stop the cluster after computation
# stopCluster(cl)
# 
# # End the timer
# end_time <- Sys.time()
# 
# # Create a matrix with columns "i", "j", and "Integral"
# results_matrix <- matrix(NA, nrow = length(II_test) * length(JJ_test), ncol = 3)
# results_matrix[, 1:2] <- as.matrix(expand.grid(II_test, JJ_test))
# results_matrix[, 3] <- unlist(results)
# colnames(results_matrix) <- c("i", "j", "Integral")
# 
# # Calculate the elapsed time in minutes
# elapsed_time_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
# 
# # Print the results matrix
# print(results_matrix)
# 
# # Print the elapsed time in minutes
# cat("Elapsed time: ", round(elapsed_time_minutes, 2), "minutes\n")

#########################
## Log-Gaussian kernel ##
#########################

# Centered SMN density

G <- function (Y, b) { # b > 0, Y > 0 is a d x d symmetric matrix
  d <- nrow(Y)
  return (exp(-sum(diag(Y %*% Y)) / (2 * b)) / ((2 * pi * b) ^ (d * (d + 1) / 4)) * 2 ^ (d * (d - 1) / 4))
}

# # Tests if G integrates to 1
# 
# # Start the timer
# start_time <- Sys.time()
# 
# # Define the integrand function
# integrand <- function(vars) {
#   # Construct the 2x2 symmetric matrix
#   Y <- matrix(c(vars[1], vars[2], vars[2], vars[3]), nrow = 2, ncol = 2)
# 
#   # Ensure the matrix is symmetric (this is automatic in our construction)
#   G_value <- G(Y, b = 1)  # Example value of b
# 
#   return(G_value)  # Return the G value as we are integrating in Lebesgue measure over R ^ 3
# }
# 
# # Define the limits of integration
# lower_limit <- c(-Inf, -Inf, -Inf)
# upper_limit <- c(Inf, Inf, Inf)
# 
# # Perform the numerical integration using the cubature package
# result <- adaptIntegrate(integrand, lowerLimit = lower_limit, upperLimit = upper_limit, tol = tol1)
# 
# # End the timer
# end_time <- Sys.time()
# 
# # Calculate the elapsed time in seconds
# elapsed_time_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
# 
# # Print the result of the integral
# print(result$integral)
# 
# # Print the computation time in seconds
# cat("Computation time:", elapsed_time_seconds, "seconds\n")

# Log-Gaussian kernel (introduced by Schwartzman (2016))

LG <- function(X, b, S) {
  # Recenter logm(S)
  Y <- logm(S) - logm(X)
  
  # Calculate the G term
  G_term <- G(Y, b) / det(S)
  
  # Compute the eigenvalues of S
  eigenvalues <- eigen(S)$values
  d <- length(eigenvalues)
  
  # Compute the product factor
  product_factor <- 1
  for (i in 1:d) {
    if (i < d) {
      for (j in (i+1):d) {
        if (abs(eigenvalues[i] - eigenvalues[j]) > .Machine$double.eps) { # > 0
          numerator <- log(eigenvalues[i]) - log(eigenvalues[j])
          denominator <- eigenvalues[i] - eigenvalues[j]
          product_factor <- product_factor * abs(numerator / denominator)
        } else {
          product_factor <- product_factor * (1 / eigenvalues[i])
        }
      }
    }
  }
  
  # Return the product of G_term and product_factor
  return(G_term * product_factor)
}

# # Tests if LG integrates to 1
# 
# # Example matrix X (must be symmetric positive definite) and value for b
# X_test <- matrix(c(2, 0.5, 0.5, 1), nrow = 2)
# b_test <- 1
# 
# # Start the timer
# start_time <- Sys.time()
# 
# # Define the integrand function
# integrand <- function(vars) {
#   # Construct the SPD matrix using the parameters theta, lambda1, and lambda2
#   S <- construct_X(vars[1], vars[2], vars[3])
# 
#   # Compute the value of the LG function with b and X as inputs
#   LG_value <- LG(X_test, b_test, S)
# 
#   # Calculate the Jacobian for the transformation
#   jacobian_value <- abs(vars[2] - vars[3]) / 4
# 
#   return(LG_value * jacobian_value)
# }
# 
# # Define the limits of integration
# lower_limit <- c(0, 0, 0)
# upper_limit <- c(2 * pi, Inf, Inf)
# 
# # Perform the numerical integration using the cubature package
# result <- adaptIntegrate(integrand, lowerLimit = lower_limit, upperLimit = upper_limit, tol = tol1)
# 
# # End the timer
# end_time <- Sys.time()
# 
# # Calculate the elapsed time in minutes
# elapsed_time_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
# 
# # Print the result of the integral
# print(result$integral)
# 
# # Print the computation time in minutes
# cat("Computation time:", round(elapsed_time_minutes, 2), "minutes\n")

################
## Estimators ##
################

# Density estimator

hat_f <- function(XX, S, b, method = "WK") {
  estimator <- switch(method,
                      "WK" = mean(sapply(XX, function(X) LaplacesDemon::dwishart(X, 1/b + d + 1, b * S))),
                      "LG" = mean(sapply(XX, function(X) LG(X, b, S))),
                      warning("Invalid method. Should be 'WK' or 'LG'.")
  )
  return(estimator)
}

# # Tests if the density estimator hat_f integrates to 1 over all j in JJ for both WK and LG methods
# 
# # Parameters for testing
# n_WK_test <- 10
# n_LG_test <- 2
# b_WK_test <- 0.1
# b_LG_test <- 1
# 
# # Generalized integrand function for both WK and LG methods
# integrand <- function(vars, XX_sample, bandwidth, method) {
#   # Construct the SPD matrix
#   S <- construct_X(vars[1], vars[2], vars[3])
#   
#   # Estimate the density using the chosen method (WK or LG)
#   density_value <- hat_f(XX_sample, S, bandwidth, method)
#   
#   # Calculate the Jacobian
#   jacobian_value <- abs(vars[2] - vars[3]) / 4
#   
#   return(density_value * jacobian_value)
# }
# 
# # Initialize parallel cluster and load necessary libraries and variables
# cl <- setup_parallel_cluster()
# 
# # Start the timer
# start_time <- Sys.time()
# 
# # Run the integral computation in parallel for all 6 target densities
# results <- parLapply(cl, JJ, function(j) {
#   cat("Calculating integrals for j =", j, "\n")
#   
#   # Generate a sample of SPD matrices for the given j
#   XX_sample_WK <- XX(j, n_WK_test)
#   XX_sample_LG <- XX(j, n_LG_test)
#   
#   # Perform the numerical integration for the WK method
#   lower_limit <- c(0, 0, 0)
#   upper_limit <- c(2 * pi, Inf, Inf)
#   result_WK <- adaptIntegrate(function(vars) integrand(vars, XX_sample_WK, b_WK_test, "WK"),
#                               lowerLimit = lower_limit, upperLimit = upper_limit, tol = tol1)
#   
#   # Perform the numerical integration for the LG method
#   result_LG <- adaptIntegrate(function(vars) integrand(vars, XX_sample_LG, b_LG_test, "LG"),
#                               lowerLimit = lower_limit, upperLimit = upper_limit, tol = tol1)
#   
#   cat("Result for j =", j, ":", "WK =", result_WK$integral, ", LG =", result_LG$integral, "\n")
#   return(c(j, result_WK$integral, result_LG$integral))
# })
# 
# # Stop the cluster after computation
# stopCluster(cl)
# 
# # End the timer
# end_time <- Sys.time()
# 
# # Convert the results to a 6x3 matrix and set column names
# results_matrix <- do.call(rbind, results)
# colnames(results_matrix) <- c("Density", "WK Integral", "LG Integral")
# 
# # Calculate the elapsed time in minutes
# elapsed_time_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
# 
# # Print the results matrix
# print(results_matrix)
# 
# # Print the elapsed time in minutes
# cat("Computation time:", round(elapsed_time_minutes, 2), "minutes\n")

# # Analytical verification that hat_f integrates to 1 over all j in JJ for the WK method
# 
# # Parameters for testing
# n_test <- 1
# b_WK_test <- 0.01
# 
# # Start the timer
# start_time <- Sys.time()
# 
# # Run the integral computation for all 6 target densities sequentially
# results <- lapply(JJ, function(j) {
#   cat("Calculating WK integral for j =", j, "\n")
#   
#   # Generate a sample of SPD matrices for the given j
#   XX_sample <- XX(j, n_test)
#   
#   # Analytical calculation for the WK method
#   integral_WK <- {
#     numerator <- CholWishart::mvgamma(1 / (2 * b_WK_test), d)
#     denominator <- (2 * b_WK_test) ^ (d * (d + 1) / 2) * CholWishart::mvgamma(1 / (2 * b_WK_test) + (d + 1) / 2, d)
#     numerator / denominator
#   }
#   
#   cat("Result for j =", j, "WK =", integral_WK, "\n")
#   return(c(j, integral_WK))
# })
# 
# # End the timer
# end_time <- Sys.time()
# 
# # Convert the results to a matrix and set column names
# results_matrix <- do.call(rbind, results)
# colnames(results_matrix) <- c("Density", "WK Integral")
# 
# # Calculate the elapsed time in seconds
# elapsed_time_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
# 
# # Print the results matrix
# print(results_matrix)
# 
# # Print the elapsed time in seconds
# cat("Computation time:", round(elapsed_time_seconds, 2), "seconds\n")

###########################################
## Criterion to optimize (exact version) ##
###########################################

# Least Squares Cross Validation (LSCV) criterion (exact version)

LSCV <- function(XX, b, i, j, K, method, tolerance = tol1) {
  integrand <- function(vars) {
    # Construct SPD matrix S
    S <- construct_X(vars[1], vars[2], vars[3])
    
    # Compute the squared difference between the estimator and the target density
    diff_squared <- (hat_f(XX, S, b, method) - f(i, j, K, S)) ^ 2
    jacobian_value <- abs(vars[2] - vars[3]) / 4
    
    # Return the product of diff_squared and the Jacobian factor
    return(diff_squared * jacobian_value)
  }
  
  # Define the limits of integration
  lower_limit <- c(0, delta, delta)
  upper_limit <- c(2 * pi, 1 / delta, 1 / delta)
  
  # Perform the numerical integration using the cubature package
  result <- adaptIntegrate(integrand, lowerLimit = lower_limit, upperLimit = upper_limit, tol = tolerance)
  
  # Return the result of the integral
  return(result$integral)
}

# Parameters for testing
II_test <- 1:2
JJ_test <- 1:3
K_test <- 2
n_test <- 10
b_test <- 0.1

# Initialize a matrix to store results for both methods
results_matrix <- matrix(NA, nrow = length(II_test) * length(JJ_test), ncol = 4)
results_matrix[, 1:2] <- as.matrix(expand.grid(II_test, JJ_test))
colnames(results_matrix) <- c("i", "j", "WK LSCV", "LG LSCV")

# Create a list of pairs in the same order as expand.grid(II_test, JJ_test)
combinations <- lapply(1:nrow(expand.grid(II_test, JJ_test)), function(k) as.numeric(expand.grid(II_test, JJ_test)[k, ]))

# Start the timer
start_time <- Sys.time()

# Loop over both methods
for (method in c("WK", "LG")) {
  
  # Initialize parallel cluster and load necessary libraries and variables
  cl <- setup_parallel_cluster()
  
  # Parallelize computation for all i and j combinations for the current method
  results <- parLapply(cl, combinations, function(params) {
    i <- params[[1]]
    j <- params[[2]]
    
    cat("Testing LSCV for i =", i, ", j =", j, ", K =", K_test, "and method =", method, "\n")
    
    # Call the LSCV function for the combination of i, j, and K_test
    result <- LSCV(XX = XX(i, j, K_test, n_test), b = b_test, i = i, j = j, K = K_test, method = method)
    
    cat("Result for i =", i, ", j =", j, ", K =", K_test, "and method =", method, ":\n", result, "\n\n")
    return(result)
  })
  
  # Store the results in the results_matrix
  if (method == "WK") {
    results_matrix[, 3] <- unlist(results)
  } else if (method == "LG") {
    results_matrix[, 4] <- unlist(results)
  }
  
  # Stop the cluster after computation
  stopCluster(cl)
}

# End the timer
end_time <- Sys.time()

# Calculate the elapsed time in minutes
elapsed_time_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))

# Print the results
print(results_matrix)

# Print the elapsed time in minutes
cat("Elapsed time: ", round(elapsed_time_minutes, 2), "minutes\n")

#################################################
## Criterion to optimize (Monte Carlo version) ##
#################################################

# Least Squares Cross Validation (LSCV) criterion (Monte Carlo version)

LSCV_MC <- function(XX, b, ii, jj, K, method) {
  n <- length(XX)
  d <- nrow(XX[[1]])
  r_d <- d * (d + 1) / 2
  
  sum1 <- 0
  sum2 <- 0
  
  if (method == "WK") {
    second_term <- 0 # initialization of second term using h-block cross-validation
    
    for (i in 1:n) {
      n_h <- 0 # number of indices j for a given i in h-block cross-validation
      sum2 <- 0 # for h-block cross-validation
      
      for (j in 1:n) {
        Xi <- XX[[i]]
        Xj <- XX[[j]]
        Sij <- Xi + Xj
        
        # Avoid numerical errors due to Stirling's formula
        term2 <- exp((-d / b - r_d) * log(2) - (r_d / 2) * log(b) +
                       lmvgamma(1 / b + (d + 1) / 2, d) -
                       2 * lmvgamma(1 / (2 * b) + (d + 1) / 2, d))
        
        # Combined determinants to avoid numerical errors
        term3 <- det(solve(Sij %*% Sij, Xi %*% Xj)) ^ (1 / (2 * b))
        term4 <- det(Sij) ^ (-(d + 1) / 2)
        
        sum1 <- sum1 + term2 * term3 * term4
        
        if (abs(j - i) > h(n)) { # h-block cross-validation
          sum2 <- sum2 + LaplacesDemon::dwishart(Xi, nu = 1/b + d + 1, S = b * Xj)
          n_h <- n_h + 1 # update the number of indices j for a given i in h-block cross-validation
        }
      }
      second_term <- second_term + 2 * sum2 / (n * n_h) # second term update for h-block cross-validation
    }
    first_term <- (2 ^ (d / b) * b ^ (-r_d / 2)) * sum1 / (n ^ 2)
    
  } else if (method == "LG") {
    second_term <- 0 # initialization of second term using h-block cross-validation
    
    for (i in 1:n) {
      n_h <- 0 # number of indices j for a given i in h-block cross-validation
      sum2 <- 0 # for h-block cross-validation
      
      for (j in 1:n) {
        Xi <- XX[[i]]
        Xj <- XX[[j]]
        Xi_log <- logm(Xi)
        Xj_log <- logm(Xj)
        
        etr_term <- exp(tr((-Xi_log ^ 2 - Xj_log ^ 2 + (Xi_log + Xj_log) ^ 2 / 2) / (2 * b)))
        sum1 <- sum1 + etr_term
        
        if (abs(j - i) > h(n)) { # h-block cross-validation
          sum2 <- sum2 + G(Xj_log - Xi_log, b)
          n_h <- n_h + 1 # update the number of indices j for a given i in h-block cross-validation
        }
      }
      second_term <- second_term + 2 * sum2 / (n * n_h) # second term update for h-block cross-validation
    }
    first_term <- sum1 / ((2 * pi * b) ^ (r_d / 2) * 2 ^ (d / 2) * n ^ 2)
    
  } else {
    stop("Invalid method. Should be 'WK' or 'LG'.")
  }
  
  return(first_term - second_term)
}

# # Parameters for testing
# II_test <- 1:2
# JJ_test <- 1:3
# K_test <- 2
# n_test <- 100
# b_test <- 0.1
# 
# # Initialize a matrix to store results for both methods
# results_matrix <- matrix(NA, nrow = length(II_test) * length(JJ_test), ncol = 4)
# results_matrix[, 1:2] <- as.matrix(expand.grid(II_test, JJ_test))
# colnames(results_matrix) <- c("i", "j", "WK LSCV_MC", "LG LSCV_MC")
# 
# # Create a list of pairs in the same order as expand.grid(II_test, JJ_test)
# combinations <- lapply(1:nrow(expand.grid(II_test, JJ_test)), function(k) as.numeric(expand.grid(II_test, JJ_test)[k, ]))
# 
# # Start the timer
# start_time <- Sys.time()
# 
# # Loop over both methods
# for (method in c("WK", "LG")) {
#   
#   # Initialize parallel cluster and load necessary libraries and variables
#   cl <- setup_parallel_cluster()
#   
#   # Parallelize computation for all i and j combinations for the current method
#   results <- parLapply(cl, combinations, function(params) {
#     i <- params[[1]]
#     j <- params[[2]]
#     
#     cat("Testing LSCV_MC for i =", i, ", j =", j, ", K =", K_test, "and method =", method, "\n")
#     
#     # Call the LSCV_MC function for the combination of i, j, and K
#     result <- LSCV_MC(XX = XX(i, j, K_test, n_test), b = b_test, ii = i, jj = j, K = K_test, method = method)
#     
#     cat("Result for i =", i, ", j =", j, ", K =", K_test, "and method =", method, ":\n", result, "\n\n")
#     return(result)
#   })
#   
#   # Store the results in the results_matrix
#   if (method == "WK") {
#     results_matrix[, 3] <- unlist(results)
#   } else if (method == "LG") {
#     results_matrix[, 4] <- unlist(results)
#   }
#   
#   # Stop the cluster after computation
#   stopCluster(cl)
# }
# 
# # End the timer
# end_time <- Sys.time()
# 
# # Calculate the elapsed time in minutes
# elapsed_time_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
# 
# # Print the results
# print(results_matrix)
# 
# # Print the elapsed time in minutes
# cat("Elapsed time: ", round(elapsed_time_minutes, 2), "minutes\n")

###############################################
## Optimal bandwidth (parallel grid version) ##
###############################################

b_opt_MC_grid <- function(XX, i, j, K, method, return_LSCV_MC = FALSE) {
  
  # Generate grid points based on the method's bandwidth sequence
  b_grid <- BB(method)  # Use BB(method) to generate the grid of bandwidths
  
  # Initialize parallel cluster
  cl <- setup_parallel_cluster()
  
  # Parallelize computation for all b values on the grid
  LSCV_MC_values <- parSapply(cl, b_grid, function(b) {
    LSCV_MC(XX, b, i, j, K, method)
  })
  
  # Stop the cluster after computation
  stopCluster(cl)
  
  # Find the index of the b that minimizes LSCV_MC_values
  min_index <- which.min(LSCV_MC_values)
  
  # Get the optimal b value and the corresponding LSCV_MC value
  b_opt_value <- b_grid[min_index]
  min_LSCV_MC_value <- LSCV_MC_values[min_index]
  
  # Return the desired value(s) based on the return_LSCV_MC argument
  if (return_LSCV_MC) {
    return(min_LSCV_MC_value)
  } else {
    return(b_opt_value)
  }
}

# # Sequentially tests the b_opt_MC_grid function over all i in II_test and j in JJ_test for both WK and LG methods
# 
# # Parameters for testing
# II_test <- 1:2
# JJ_test <- 1:3
# K_test <- 4
# n_test <- 100
# 
# # Initialize a matrix to store results for both methods
# results_matrix <- matrix(NA, nrow = length(II_test) * length(JJ_test), ncol = 4)
# results_matrix[, 1:2] <- as.matrix(expand.grid(II_test, JJ_test))
# colnames(results_matrix) <- c("i", "j", "WK b_opt_MC_grid", "LG b_opt_MC_grid")
# 
# # Start the timer
# start_time <- Sys.time()
# 
# # Loop over both methods (WK and LG) sequentially
# for (method in c("WK", "LG")) {
# 
#   # Loop sequentially over all i and j combinations
#   for (i in II_test) {
#     for (j in JJ_test) {
#       cat("Calculating b_opt_MC_grid for i =", i, ", j =", j, ", K =", K_test, "and method =", method, "\n")
# 
#       # Generate sample matrices for the given i, j
#       XX_sample <- XX(i, j, K_test, n_test)
# 
#       # Calculate the optimal bandwidth using the grid search
#       result <- b_opt_MC_grid(XX_sample, i, j, K_test, method)
# 
#       # Store the result in the appropriate column
#       if (method == "WK") {
#         results_matrix[which(results_matrix[, 1] == i & results_matrix[, 2] == j), 3] <- result
#       } else if (method == "LG") {
#         results_matrix[which(results_matrix[, 1] == i & results_matrix[, 2] == j), 4] <- result
#       }
# 
#       cat("Result for i =", i, ", j =", j, "and method =", method, ":\n", result, "\n\n")
#     }
#   }
# }
# 
# # End the timer
# end_time <- Sys.time()
# 
# # Calculate the elapsed time in minutes
# elapsed_time_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
# 
# # Print the results
# print(results_matrix)
# 
# # Print the elapsed time in minutes
# cat("Elapsed time: ", round(elapsed_time_minutes, 2), "minutes\n")

#######################################
## Plotting LSCV_MC vs b_opt_MC_grid ##
#######################################

# # Set parameters for the test
# II_test <- 1:2
# JJ_test <- 1:3
# K_test <- 4
# n_test <- 100
# results_grid <- list()
# 
# # Start the timer for the entire process
# start_time_total <- Sys.time()
# 
# # Loop over all i, j values, and methods
# for (i in II_test) {
#   for (j in JJ_test) {
#     cat(paste("\nPlotting LSCV_MC vs b_opt_MC_grid for i =", i, ", j =", j, "\n"))
#     
#     results_grid[[paste0(i, "_", j)]] <- list()
#     
#     for (method in MM) {
#       cat(paste("Method:", method, "\n"))
#       
#       # Start the timer for the current combination of i, j, and method
#       start_time_individual <- Sys.time()
#       
#       # Generate sample matrices for Monte Carlo integration
#       XX_sample <- XX(i, j, K_test, n_test)
#       
#       # Initialize parallel cluster for this combination
#       cl <- setup_parallel_cluster()
#       
#       # Generate grid points based on the method's bandwidth sequence
#       b_values <- BB(method)
#       
#       # Calculate LSCV_MC for different values of b in parallel
#       LSCV_MC_values <- unlist(parLapply(cl, b_values, function(b) {
#         LSCV_MC(XX_sample, b, i, j, K_test, method)
#       }))
#       
#       # Stop the cluster after the parallel computation for LSCV_MC
#       stopCluster(cl)
#       
#       # Use the b_opt_MC_grid function to find the optimal b value (with internal parallelization)
#       b_opt_MC_grid_value <- b_opt_MC_grid(XX_sample, i, j, K_test, method)
#       
#       # Create a data frame for plotting
#       results_df <- data.frame(b = b_values, LSCV_MC = LSCV_MC_values)
#       
#       # Plot LSCV_MC values and indicate the optimal b
#       plot <- ggplot(results_df, aes(x = b, y = LSCV_MC)) +
#         geom_line(color = "blue") +
#         geom_vline(xintercept = b_opt_MC_grid_value, color = "red", linetype = "dashed") +
#         labs(
#           title = paste("LSCV_MC for i =", i, ", j =", j, "and method =", method),
#           x = "b",
#           y = "LSCV_MC"
#         ) +
#         theme_minimal()
#       
#       # Save the plot as a PDF file
#       pdf(file = file.path(path, paste0("LSCV_MC_vs_b_opt_MC_grid_i", i, "_j", j, "_", method, ".pdf")))
#       print(plot)
#       dev.off()
#       
#       # End the timer for the current combination of i, j, and method
#       end_time_individual <- Sys.time()
#       
#       # Calculate elapsed time in minutes for the current combination
#       elapsed_time_individual <- as.numeric(difftime(end_time_individual, start_time_individual, units = "mins"))
#       
#       # Store the result and time taken
#       results_grid[[paste0(i, "_", j)]][[method]] <- list(
#         b_opt_value = b_opt_MC_grid_value,
#         time_taken = elapsed_time_individual
#       )
#       
#       # Print the result and the time taken
#       cat(paste("b_opt_MC_grid_value for method", method, "with i =", i, ", j =", j, ":", b_opt_MC_grid_value, "\n"))
#       cat(paste("Time taken for method", method, "with i =", i, ", j =", j, ":", round(elapsed_time_individual, 2), "minutes\n"))
#     }
#   }
# }
# 
# # End the timer for the entire process
# end_time_total <- Sys.time()
# 
# # Calculate total elapsed time in minutes
# elapsed_time_total <- as.numeric(difftime(end_time_total, start_time_total, units = "mins"))
# 
# # Print the total elapsed time in minutes
# cat("Total elapsed time: ", round(elapsed_time_total, 2), "minutes\n")

#####################################
## Integrated Squared Errors (ISE) ##
#####################################

ISE <- function(XX, i, j, K, method, tolerance = tol1) {
  b_opt_MC_grid_value <- b_opt_MC_grid(XX, i, j, K, method)
  
  return(LSCV(XX, b_opt_MC_grid_value, i, j, K, method, tolerance))
}

# # Parameters for testing
# II_test <- 1:2
# JJ_test <- 1:3
# K_test <- 4
# n_test <- 500
# tol_test <- 0.01
# 
# # Initialize a matrix to store results for both methods
# results_matrix <- matrix(NA, nrow = length(II_test) * length(JJ_test), ncol = 4)
# results_matrix[, 1:2] <- as.matrix(expand.grid(II_test, JJ_test))
# colnames(results_matrix) <- c("i", "j", "WK ISE", "LG ISE")
# 
# # Start the timer
# start_time <- Sys.time()
# 
# # Loop over both methods (WK and LG)
# for (method in c("WK", "LG")) {
#   
#   # Loop over all i and j combinations sequentially
#   for (combo in 1:nrow(expand.grid(II_test, JJ_test))) {
#     i <- as.numeric(expand.grid(II_test, JJ_test)[combo, 1])
#     j <- as.numeric(expand.grid(II_test, JJ_test)[combo, 2])
#     
#     cat("Calculating ISE for i =", i, ", j =", j, "and method =", method, "\n")
#     
#     # Generate sample matrices for the given combination of i, j, and K
#     XX_sample <- XX(i, j, K_test, n_test)
#     
#     # Calculate the ISE using the optimal bandwidth from b_opt_MC_grid
#     result <- ISE(XX_sample, i, j, K_test, method, tol_test)
#     
#     # Store the result in the appropriate column
#     if (method == "WK") {
#       results_matrix[combo, 3] <- result
#     } else if (method == "LG") {
#       results_matrix[combo, 4] <- result
#     }
#     
#     cat("ISE for i =", i, ", j =", j, "and method =", method, ":\n", result, "\n\n")
#   }
# }
# 
# # End the timer
# end_time <- Sys.time()
# 
# # Calculate the elapsed time in minutes
# elapsed_time_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
# 
# # Print the results
# print(results_matrix)
# 
# # Print the elapsed time in minutes
# cat("Elapsed time: ", round(elapsed_time_minutes, 2), "minutes\n")

###############################
## Main code (exact version) ##
###############################

.libPaths("~/R/library")

# Disable the check for random number generation misuse in doFuture
options(doFuture.rng.onMisuse = "ignore")

# Register the doFuture parallel backend
registerDoFuture()

# Tweak the batchtools_slurm with the custom template and resources
myslurm <- tweak(
  batchtools_slurm,
  template = "batchtools.slurm.dependent.tmpl",
  resources = resources_list
)

# Set the plan for future
plan(list(myslurm, multisession))

# Create empty data frames to store the results
raw_results <- data.frame(
  n = integer(),
  i = integer(),
  j = integer(),
  method = character(),
  ISE = numeric(),
  stringsAsFactors = FALSE
)

# Capture the start time
start_time <- Sys.time()

# Parallel loop over the replications (RR), each node processes one set of RR values
res <- foreach(r = RR, .combine = "rbind", 
               .export = vars_to_export,
               .packages = libraries_to_load) %dopar% {
                 # Set a unique seed for each node (replication)
                 set.seed(r)
                 
                 # Set library paths within each worker node
                 .libPaths("~/R/library")
                 
                 local_raw_results <- data.frame(
                   n = integer(),
                   i = integer(),
                   j = integer(),
                   method = character(),
                   ISE = numeric(),
                   stringsAsFactors = FALSE
                 )
                 
                 # Loop over combinations of i, j, n, and method within each worker
                 for (i in II) {
                   for (j in JJ) {
                     for (n in NN) {
                       for (method in MM) {
                         XX_data <- XX(i, j, K, n)
                         ISE_value <- ISE(XX_data, i, j, K, method)
                         
                         # Store the result for this specific replication
                         local_raw_results <- rbind(
                           local_raw_results,
                           data.frame(
                             n = n,
                             i = i,
                             j = j,
                             method = method,
                             ISE = ISE_value,
                             stringsAsFactors = FALSE
                           )
                         )
                       }
                     }
                   }
                 }
                 
                 # Return the raw results for this replication
                 return(local_raw_results)
               }

# Combine results from all nodes
raw_results <- res

# Stop parallel execution
plan(sequential)

# Calculate the duration in minutes
elapsed_time_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
print(paste("Elapsed time:", round(elapsed_time_minutes, 2), "minutes"))

# Save the raw results to a CSV file
raw_output_file <- file.path(path, "raw_ISE_dependent.csv")
write.csv(raw_results, raw_output_file, row.names = FALSE)

print("Raw results saved to raw_ISE_dependent.csv")

#########################
## Process the results ##
#########################

# Create a data frame to store the summary results
summary_results <- data.frame(
  n = integer(),
  i = integer(),
  j = integer(),
  method = character(),
  mean_ISE = numeric(),
  sd_ISE = numeric(),
  median_ISE = numeric(),
  IQR_ISE = numeric(),
  stringsAsFactors = FALSE
)

# Loop through the results to compute the summary statistics
for (i in II) {
  for (j in JJ) {
    for (n in NN) {
      for (method in MM) {
        filtered_results <- raw_results %>%
          filter(i == !!i, j == !!j, n == !!n, method == !!method)
        
        ISE_values <- filtered_results$ISE
        mean_ISE <- mean(ISE_values)
        sd_ISE <- sd(ISE_values)
        median_ISE <- median(ISE_values)
        IQR_ISE <- IQR(ISE_values)
        
        # Store the summary results
        summary_results <- rbind(
          summary_results,
          data.frame(
            n = n,
            i = i,
            j = j,
            method = method,
            mean_ISE = mean_ISE,
            sd_ISE = sd_ISE,
            median_ISE = median_ISE,
            IQR_ISE = IQR_ISE,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}

# Save the summary results to a CSV file in the specified path
summary_output_file <- file.path(path, "ISE_results_dependent.csv")
write.csv(summary_results, summary_output_file, row.names = FALSE)

print("Summary results saved to ISE_results_dependent.csv")
