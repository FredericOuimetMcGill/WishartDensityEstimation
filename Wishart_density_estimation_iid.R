####################################################################################################
## Density estimation on the space of positive definite matrices using Wishart asymmetric kernels ##
####################################################################################################

## Written by Frederic Ouimet (March 2025)

require("CholWishart")        # for the multivariate gamma function
require("cubature")           # for integrals
require("doFuture")           # for parallel execution with foreach
require("expm")               # for matrix logarithms and exponentials
require("fs")                 # for filesystem operations (dependency)
require("future.batchtools")  # for batchtools integration with future
require("ggplot2")            # for plotting
require("LaplacesDemon")      # for the Wishart and inverse Wishart distributions
require("matrixcalc")         # for functions like is.positive.definite
require("matrixsampling")     # for the matrix-type-II Beta distribution
require("optimx")             # for improved bandwidth optimization
require("parallel")           # for parallelization of calculations
require("Rcpp")               # for C++ functions
require("RcppEigen")          # for C++ array/matrix algebra
require("tidyverse")          # for data manipulation and visualization
require("writexl")            # to write output to Excel files

##############################
## Parallelization on cores ##
##############################

# Define the list of libraries to load on each cluster node

libraries_to_load <- c(
  "CholWishart",
  "cubature",
  "doFuture",
  "expm",
  "future.batchtools",
  "fs",
  "ggplot2",
  "LaplacesDemon",
  "matrixcalc",
  "matrixsampling",
  "optimx",
  "parallel",
  "Rcpp",
  "RcppEigen",
  "tidyverse",
  "writexl"
)

# Define the list of variables/functions to export to the worker nodes

vars_to_export <- c(
  "BB",
  "CholWishart",
  "G",
  "ISE",
  "ISE_MC",
  "ISE_value",
  "JJ",
  "LG",
  "LSCV",
  "LSCV_MC",
  "LaplacesDemon",
  "MM",
  "NN",
  "RR",
  "S1",
  "S2",
  "S3",
  "S4",
  "XX",
  "XX_data",
  "XX_sample",
  "b_LG_test",
  "b_WK_test",
  "b_opt_MC",
  "b_opt_MC_grid",
  "b_test",
  "compute_ISE",
  "construct_X",
  "cores_per_node",
  "cubature",
  "d",
  "delta",
  "dmatrixbeta_typeII",
  "elapsed_time_seconds", # (***)
  "expm",
  "f",
  "hat_f",
  "integrand",
  "j",
  "local_raw_results",
  "lmvgamma",
  "logm",
  "matrixcalc",
  "matrixsampling",
  "method",
  "n_LG_test",
  "n_WK_test",
  "n_test",
  "optimx",
  "path",
  "raw_results",
  "resources_list",
  "rotation_matrix",
  "setup_parallel_cluster",
  "symmetrize",
  "test_estimator_integral",
  "time_seconds", # (***)
  "tol1",
  "tol_test",
  "vars_to_export",
  "writexl"
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
  
  # Each parallel worker must compile & load these C++ functions
  invisible(clusterEvalQ(cl, {
    sourceCpp("3d_array_inversion.cpp")
    sourceCpp("3d_array_ratio_conjugation.cpp")
    sourceCpp("3d_array_symmetrization.cpp")
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
#   "C://Users//fred1//Dropbox//Ouimet_Genest_projects//GMOR_2024_Wishart_asym_kernels//_simulations", 
#   fsep = .Platform$file.sep
# )
path <- getwd()
# setwd(path)

sourceCpp("3d_array_inversion.cpp")
sourceCpp("3d_array_ratio_conjugation.cpp")
sourceCpp("3d_array_symmetrization.cpp")

################
## Parameters ##
################

d <- 2 # width of the square matrices
delta <- 0.1 # lower bound on the eigenvalues of SPD matrices in LSCV

MM <- list("WK", "LG") # list of density estimation methods
NN <- c(10, 20) # sample sizes
JJ <- 1:6 # target density function indices
RR <- 1:1 # replication indices

cores_per_node <- 31 # number of cores for each node in the super-computer

tol1 <- 1e-1

##############################
## Parallelization on nodes ##
##############################

resources_list <- list(
  cpus_per_task = cores_per_node,
  mem = "240G",
  walltime = "4:30:00",
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

# Construct X = rotation_matrix(theta) . diag(lambda1, lambda2) . rotation_matrix(theta) ^ T

construct_X <- function(theta, lambda1, lambda2) {
  R_theta <- rotation_matrix(theta)
  diag_lambda <- diag(c(lambda1, lambda2))
  X <- R_theta %*% diag_lambda %*% t(R_theta)
  return(symmetrize(X))
}

# Bandwidth sequence based on the method

BB <- function(method) {
  if (method == "WK") {
    return(seq(0.02, 2, length.out = cores_per_node))  # BB_WK sequence
  } else if (method == "LG") {
    return(seq(0.02, 2, length.out = cores_per_node))  # BB_LG sequence
  } else {
    stop("Invalid method. Should be 'WK' or 'LG'.")
  }
}

##################
## Observations ##
##################

# Scale matrices

S1 <- matrix(c(1, 0.1, 0.1, 1), nrow = 2)
S2 <- matrix(c(1, -0.9, -0.9, 1), nrow = 2)
S3 <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)
S4 <- matrix(c(1, -0.5, -0.5, 1), nrow = 2)

# Random generation of observations

XX <- function(j, n) {
  if (j == 1) {
    bern <- as.numeric(runif(n) < 0.5)
    arr1 <- stats::rWishart(n, 4, S1)
    arr2 <- stats::rWishart(n, 5, S2)
    out <- arr1 * array(rep(bern, each = d * d), dim = c(d, d, n)) +
      arr2 * array(rep(1 - bern, each = d * d), dim = c(d, d, n))
    symmetrize_3d_array(out)
    
  } else if (j == 2) {
    bern <- as.numeric(runif(n) < 0.5)
    arr1 <- stats::rWishart(n, 5, S3)
    arr2 <- stats::rWishart(n, 6, S4)
    out <- arr1 * array(rep(bern, each = d * d), dim = c(d, d, n)) +
      arr2 * array(rep(1 - bern, each = d * d), dim = c(d, d, n))
    symmetrize_3d_array(out)
    
  } else if (j == 3) {
    arr <- stats::rWishart(n, 5, S2)
    symmetrize_3d_array(invert_3d_array(arr))
    
  } else if (j == 4) {
    arr <- stats::rWishart(n, 6, S3)
    symmetrize_3d_array(invert_3d_array(arr))
    
  } else if (j == 5) {
    arr1 <- stats::rWishart(n, 2 * 2, diag(d))
    arr2 <- stats::rWishart(n, 2 * 4, diag(d))
    symmetrize_3d_array(conjugated_ratio_3d_array(arr1, arr2))
    
  } else if (j == 6) {
    arr1 <- stats::rWishart(n, 2 * 3, diag(d))
    arr2 <- stats::rWishart(n, 2 * 5, diag(d))
    symmetrize_3d_array(conjugated_ratio_3d_array(arr1, arr2))
    
  } else {
    warning("Invalid value of j. Should be 1, 2, 3, 4, 5, or 6.")
    NULL
  }
}

##############################
## Target density functions ##
##############################

# Density function for matrix-type-II Beta distribution

dmatrixbeta_typeII <- function(X, a, b) { # X is an SPD matrix of size d x d, and a,b > (d - 1)/2
  d <- nrow(X)
  numerator <- mvgamma(a + b, d)
  denominator <- mvgamma(a, d) * mvgamma(b, d)
  det_X <- det(X)
  det_I_plus_X <- det(diag(d) + X)
  
  density <- (numerator / denominator) * det_X ^ (a - (d + 1)/2) * det_I_plus_X ^ (-(a + b))
  return(density)
}

# # Tests if dmatrixbeta_typeII(X, a, b) integrates to 1
# 
# # Parameters for testing
# a_test <- 2
# b_test <- 4
# tol_test <- 1e-3
# 
# # Start the timer
# start_time <- Sys.time()
# 
# # Define the integrand for cubature
# integrand <- function(vars) {
#   # Construct SPD matrix X from the parameters theta, lambda1, lambda2
#   X <- construct_X(vars[1], vars[2], vars[3])
# 
#   # Ensure X is positive definite
#   if (!is.positive.definite(X)) {
#     return(0)
#   }
# 
#   # Compute the density using the dmatrixbeta_typeII function
#   density_value <- dmatrixbeta_typeII(X, a_test, b_test)
# 
#   # Return the density multiplied by the Jacobian adjustment for polar coordinates
#   jacobian_value <- abs(vars[2] - vars[3]) / 4
#   return(density_value * jacobian_value)
# }
# 
# # Set integration limits for theta (0 to 2*pi), lambda1, and lambda2 (positive values)
# lower_limit_test <- c(0, 0, 0)
# upper_limit_test <- c(2 * pi, Inf, Inf)
# 
# # Perform the numerical integration using adaptIntegrate from the cubature package
# result_test <- adaptIntegrate(integrand, lowerLimit = lower_limit_test, upperLimit = upper_limit_test, tol = tol_test)
# 
# # End the timer
# end_time <- Sys.time()
# 
# # Calculate the elapsed time in seconds
# elapsed_time_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
# 
# # Print the result of the integral
# cat("The integral of dmatrixbeta_typeII with test parameters is approximately:", result_test$integral, "\n")
# 
# # Print the computation time in seconds
# cat("Computation time:", elapsed_time_seconds, "seconds\n")

# Target densities

f <- function(j, X) { # X is an SPD matrix of size d x d
  if (j == 1) {
    # Case when j = 1
    res <- 0.5 * LaplacesDemon::dwishart(X, 4, S1) + 0.5 * LaplacesDemon::dwishart(X, 5, S2)
  } else if (j == 2) {
    # Case when j = 2
    res <- 0.5 * LaplacesDemon::dwishart(X, 5, S3) + 0.5 * LaplacesDemon::dwishart(X, 6, S4)
  } else if (j == 3) {
    # Case when j = 3
    res <- LaplacesDemon::dinvwishart(X, 5, S2)
  } else if (j == 4) {
    # Case when j = 4
    res <- LaplacesDemon::dinvwishart(X, 6, S3)
  } else if (j == 5) {
    # Case when j = 5
    res <- dmatrixbeta_typeII(X, 2, 4)
  } else if (j == 6) {
    # Case when j = 6
    res <- dmatrixbeta_typeII(X, 3, 5)
  } else {
    warning("Invalid value of j. Should be 1, 2, 3, 4, 5, or 6.")
    res <- NULL
  }
  return(res)
}

# # Tests if the target density f(j, X) integrates to 1 over all j in JJ
# 
# tol_test <- 1e-3 # Test value for the tolerance
# 
# # Initialize parallel cluster and load necessary libraries and variables
# cl <- setup_parallel_cluster()
# 
# # Start the timer
# start_time <- Sys.time()
# 
# # Run the integral computation in parallel for all 6 target densities
# results <- parLapply(cl, JJ, function(j) {
#   cat("Calculating integral for j =", j, "\n")
# 
#   # Integrand function for the triple integral
#   integrand <- function(vars) {
#     density_value <- f(j, construct_X(vars[1], vars[2], vars[3]))
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
#   result <- adaptIntegrate(integrand, lowerLimit = lower_limit, upperLimit = upper_limit, tol = tol_test)
# 
#   cat("Result for j =", j, ":", result$integral, "\n")
#   return(result$integral)
# })
# 
# # Stop the cluster after computation
# stopCluster(cl)
# 
# # End the timer
# end_time <- Sys.time()
# 
# # Create a 6x2 matrix with columns "Density" and "Integral"
# results_matrix <- matrix(NA, nrow = 6, ncol = 2)
# results_matrix[, 1] <- JJ
# results_matrix[, 2] <- unlist(results)
# colnames(results_matrix) <- c("Density", "Integral")
# 
# # Calculate the elapsed time in minutes
# elapsed_time_minutes <- as.numeric(difftime(end_time, start_time, units = "secs"))
# 
# # Print the results as a 6x2 matrix
# print(results_matrix)
# 
# # Print the elapsed time in minutes
# cat("Elapsed time: ", round(elapsed_time_minutes, 2), "seconds\n")

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
# tol_test <- 1e-3 # Test value for the tolerance
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
# result <- adaptIntegrate(integrand, lowerLimit = lower_limit, upperLimit = upper_limit, tol = tol_test)
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
# # Parameters for testing
# X_test <- matrix(c(2, 0.5, 0.5, 1), nrow = 2)
# b_test <- 1
# tol_test <- 1e-3
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
# result <- adaptIntegrate(integrand, lowerLimit = lower_limit, upperLimit = upper_limit, tol = tol_test)
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
                      "WK" = mean(apply(XX, 3, function(X) LaplacesDemon::dwishart(X, 1/b + d + 1, b * S))),
                      "LG" = mean(apply(XX, 3, function(X) LG(X, b, S))),
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
# tol_test <- 1e-3
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
#                               lowerLimit = lower_limit, upperLimit = upper_limit, tol = tol_test)
# 
#   # Perform the numerical integration for the LG method
#   result_LG <- adaptIntegrate(function(vars) integrand(vars, XX_sample_LG, b_LG_test, "LG"),
#                               lowerLimit = lower_limit, upperLimit = upper_limit, tol = tol_test)
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

LSCV <- function(XX, b, j, method, tolerance = tol1) {
  integrand <- function(vars) {
    # Construct SPD matrix S
    S <- construct_X(vars[1], vars[2], vars[3])
    
    diff_squared <- (hat_f(XX, S, b, method) - f(j, S)) ^ 2
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

# # Tests the LSCV function (exact version) over all j in JJ for both WK and LG methods
# 
# # Parameters for testing
# n_test <- 10
# b_test <- 0.1
# 
# # Initialize a matrix to store results for both methods
# results_matrix <- matrix(NA, nrow = 6, ncol = 3)
# results_matrix[, 1] <- JJ
# colnames(results_matrix) <- c("Density", "WK LSCV", "LG LSCV")
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
#   # Parallelize computation for all j's for the current method
#   results <- parLapply(cl, JJ, function(j) {
#     cat("Testing LSCV for j =", j, "and method =", method, "\n")
#     result <- LSCV(XX = XX(j, n_test), b = b_test, j = j, method = method)
#     cat("Result for j =", j, "and method =", method, ":\n", result, "\n\n")
#     return(result)
#   })
# 
#   # Store the results in the results_matrix
#   if (method == "WK") {
#     results_matrix[, 2] <- unlist(results)
#   } else if (method == "LG") {
#     results_matrix[, 3] <- unlist(results)
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

#################################################
## Criterion to optimize (Monte Carlo version) ##
#################################################

# Least Squares Cross Validation (LSCV) criterion (Monte Carlo version)

LSCV_MC <- function(XX, b, jj, method) {
  dims <- dim(XX)
  n <- dims[3]
  d <- dims[1]
  r_d <- d * (d + 1) / 2
  
  sum1 <- 0
  sum2 <- 0
  
  if (method == "WK") {
    for (i in 1:n) {
      for (j in 1:n) {
        Xi <- XX[,,i]
        Xj <- XX[,,j]
        Sij <- Xi + Xj
        
        # Avoid numerical errors due to Stirling's formula
        term2 <- exp((-d / b - r_d) * log(2) - (r_d / 2) * log(b) +
                       lmvgamma(1 / b + (d + 1) / 2, d) -
                       2 * lmvgamma(1 / (2 * b) + (d + 1) / 2, d))
        
        # Combined determinants to avoid numerical errors
        term3 <- det(solve(Sij %*% Sij, Xi %*% Xj)) ^ (1 / (2 * b))
        term4 <- det(Sij) ^ (-(d + 1) / 2)
        
        sum1 <- sum1 + term2 * term3 * term4
        
        if (i != j) {
          sum2 <- sum2 + LaplacesDemon::dwishart(Xi, nu = 1/b + d + 1, S = b * Xj)
        }
      }
    }
    first_term <- (2 ^ (d / b) * b ^ (-r_d / 2)) * sum1 / (n ^ 2)
    second_term <- 2 * sum2 / (n * (n - 1))
    
  } else if (method == "LG") {
    for (i in 1:n) {
      for (j in 1:n) {
        Xi <- XX[,,i]
        Xj <- XX[,,j]
        Xi_log <- logm(Xi)
        Xj_log <- logm(Xj)
        
        etr_term <- exp(tr((-Xi_log ^ 2 - Xj_log ^ 2 + (Xi_log + Xj_log) ^ 2 / 2) / (2 * b)))
        sum1 <- sum1 + etr_term
        
        if (i != j) {
          sum2 <- sum2 + G(Xj_log - Xi_log, b)
        }
      }
    }
    first_term <- sum1 / ((2 * pi * b) ^ (r_d / 2) * 2 ^ (d / 2) * n ^ 2)
    second_term <- 2 * sum2 / (n * (n - 1))
    
  } else {
    stop("Invalid method. Should be 'WK' or 'LG'.")
  }
  
  return(first_term - second_term)
}

# # Tests the LSCV_MC function (Monte Carlo version) over all j in JJ for both WK and LG methods
# 
# # Parameters for testing
# n_test <- 100
# b_test <- 0.1
# 
# # Initialize a matrix to store results for both methods
# results_matrix <- matrix(NA, nrow = 6, ncol = 3)
# results_matrix[, 1] <- JJ
# colnames(results_matrix) <- c("Density", "WK LSCV_MC", "LG LSCV_MC")
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
#   # Parallelize computation for all j's for the current method
#   results <- parLapply(cl, JJ, function(j) {
#     cat("Testing LSCV_MC for j =", j, "and method =", method, "\n")
#     result <- LSCV_MC(XX = XX(j, n_test), b = b_test, j = j, method = method)
#     cat("Result for j =", j, "and method =", method, ":\n", result, "\n\n")
#     return(result)
#   })
# 
#   # Store the results in the results_matrix
#   if (method == "WK") {
#     results_matrix[, 2] <- unlist(results)
#   } else if (method == "LG") {
#     results_matrix[, 3] <- unlist(results)
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

#######################
## Optimal bandwidth ##
#######################

b_opt_MC <- function(XX, j, method) {
  optimize(
    function(b) LSCV_MC(XX, b, j, method), 
    interval = c(min(BB(method)), max(BB(method)))
  )$minimum
}

# # Tests b_opt_MC over all j in JJ for both WK and LG methods
# 
# # Parameters for testing
# n_test <- 100
# 
# # Initialize a matrix to store results for both methods
# results_matrix <- matrix(NA, nrow = 6, ncol = 3)
# results_matrix[, 1] <- JJ
# colnames(results_matrix) <- c("Density", "WK b_opt_MC", "LG b_opt_MC")
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
#   # Parallelize computation for all j's for the current method
#   results <- parLapply(cl, JJ, function(j) {
#     cat("Calculating b_opt_MC for j =", j, "and method =", method, "\n")
# 
#     XX_sample <- XX(j, n_test)
#     result <- b_opt_MC(XX_sample, j, method)
# 
#     cat("Result for j =", j, "and method =", method, ":\n", result, "\n\n")
#     return(result)
#   })
# 
#   # Store the results in the results_matrix
#   if (method == "WK") {
#     results_matrix[, 2] <- unlist(results)
#   } else if (method == "LG") {
#     results_matrix[, 3] <- unlist(results)
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

b_opt_MC_grid <- function(XX, j, method, return_LSCV_MC = FALSE) {
  
  # Generate grid points based on the method's bandwidth sequence
  b_grid <- BB(method)  # Use BB(method) to generate the grid of bandwidths
  
  # Initialize parallel cluster
  cl <- setup_parallel_cluster()
  
  # Parallelize computation for all b values on the grid
  LSCV_MC_values <- parSapply(cl, b_grid, function(b) {
    LSCV_MC(XX, b, j, method)
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
# # Tests the b_opt_MC_grid function over all j in JJ for both WK and LG methods

# # Parameters for testing
# n_test <- 100
# 
# # Initialize a matrix to store results for both methods
# results_matrix <- matrix(NA, nrow = length(JJ), ncol = 3)
# results_matrix[, 1] <- JJ
# colnames(results_matrix) <- c("Density", "WK b_opt_MC_grid", "LG b_opt_MC_grid")
# 
# # Start the timer
# start_time <- Sys.time()
# 
# # Loop over both methods (WK and LG)
# for (method in c("WK", "LG")) {
# 
#   # Loop over all j's sequentially
#   for (j in JJ) {
#     cat("Calculating b_opt_MC_grid for j =", j, "and method =", method, "\n")
# 
#     # Generate sample matrices for the given j
#     XX_sample <- XX(j, n_test)
# 
#     # Calculate the optimal bandwidth using the grid search
#     result <- b_opt_MC_grid(XX_sample, j, method)
# 
#     # Store the result in the appropriate column
#     if (method == "WK") {
#       results_matrix[JJ == j, 2] <- result
#     } else if (method == "LG") {
#       results_matrix[JJ == j, 3] <- result
#     }
# 
#     cat("Result for j =", j, "and method =", method, ":\n", result, "\n\n")
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
# n_test <- 100  # Set sample size
# results_grid <- list()  # To store the results
# 
# # Start the timer for the entire process
# start_time_total <- Sys.time()
# 
# # Loop over all methods and j values
# for (j in JJ) {
#   cat(paste("\nPlotting LSCV_MC vs b_opt_MC_grid for j =", j, "\n"))
#   
#   results_grid[[j]] <- list()
#   
#   for (method in MM) {
#     cat(paste("Method:", method, "\n"))
#     
#     # Start the timer for the current combination of j and method
#     start_time_individual <- Sys.time()
#     
#     # Generate sample matrices for Monte Carlo integration
#     XX_sample <- XX(j, n_test)
#     
#     # Initialize parallel cluster for this combination
#     cl <- setup_parallel_cluster()
#     
#     # Generate grid points based on the method's bandwidth sequence
#     b_values <- BB(method)
#     
#     # Calculate LSCV_MC for different values of b in parallel
#     LSCV_MC_values <- unlist(parLapply(cl, b_values, function(b) {
#       LSCV_MC(XX_sample, b, j, method)
#     }))
#     
#     # Stop the cluster after the parallel computation for LSCV_MC
#     stopCluster(cl)
#     
#     # Use the b_opt_MC_grid function to find the optimal b value (with internal parallelization)
#     b_opt_MC_grid_value <- b_opt_MC_grid(XX_sample, j, method)
#     
#     # Create a data frame for plotting
#     results_df <- data.frame(b = b_values, LSCV_MC = LSCV_MC_values)
#     
#     # Plot LSCV_MC values and indicate the optimal b
#     plot <- ggplot(results_df, aes(x = b, y = LSCV_MC)) +
#       geom_line(color = "blue") +
#       geom_vline(xintercept = b_opt_MC_grid_value, color = "red", linetype = "dashed") +
#       labs(
#         title = paste("LSCV_MC for j =", j, "and method =", method),
#         x = "b",
#         y = "LSCV_MC"
#       ) +
#       theme_minimal()
#     
#     # Save the plot as a PDF file
#     pdf(file = file.path(path, paste0("LSCV_MC_vs_b_opt_MC_grid_j", j, "_", method, ".pdf")))
#     print(plot)
#     dev.off()
#     
#     # End the timer for the current combination of j and method
#     end_time_individual <- Sys.time()
#     
#     # Calculate elapsed time in minutes for the current combination
#     elapsed_time_individual <- as.numeric(difftime(end_time_individual, start_time_individual, units = "mins"))
#     
#     # Store the result and time taken
#     results_grid[[j]][[method]] <- list(
#       b_opt_value = b_opt_MC_grid_value,
#       time_taken = elapsed_time_individual
#     )
#     
#     # Print the result and the time taken
#     cat(paste("b_opt_MC_grid_value for method", method, "with j =", j, ":", b_opt_MC_grid_value, "\n"))
#     cat(paste("Time taken for method", method, "with j =", j, ":", round(elapsed_time_individual, 2), "minutes\n"))
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

ISE <- function(XX, j, method, tolerance = tol1) {
  b_opt_MC_grid_value <- b_opt_MC_grid(XX, j, method)
    
  return(LSCV(XX, b_opt_MC_grid_value, j, method, tolerance))
}

# # Tests the ISE function over all j in JJ for both WK and LG methods
# 
# # Parameters for testing
# n_test <- 500
# 
# # Initialize a matrix to store results for both methods
# results_matrix <- matrix(NA, nrow = length(JJ), ncol = 3)
# results_matrix[, 1] <- JJ
# colnames(results_matrix) <- c("Density", "WK ISE", "LG ISE")
# 
# # Start the timer
# start_time <- Sys.time()
# 
# # Loop over both methods (WK and LG)
# for (method in c("WK", "LG")) {
#   
#   # Loop over all j's sequentially
#   for (j in JJ) {
#     cat("Calculating ISE for j =", j, "and method =", method, "\n")
#     
#     # Generate sample matrices for the given j
#     XX_sample <- XX(j, n_test)
#     
#     # Calculate the ISE using the optimal bandwidth from b_opt_MC_grid
#     result <- ISE(XX_sample, j, method)
#     
#     # Store the result in the appropriate column
#     if (method == "WK") {
#       results_matrix[JJ == j, 2] <- result
#     } else if (method == "LG") {
#       results_matrix[JJ == j, 3] <- result
#     }
#     
#     cat("ISE for j =", j, "and method =", method, ":\n", result, "\n\n")
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
  template = "batchtools.slurm.iid.tmpl",
  resources = resources_list
)

# Set the plan for future
plan(list(myslurm, multisession))

# Create empty data frames to store the results
raw_results <- data.frame(
  j = integer(),
  n = integer(),
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
    j = integer(),
    n = integer(),
    method = character(),
    ISE = numeric(),
    time_seconds = numeric(), # New column for elapsed time (***)
    stringsAsFactors = FALSE
  )

  # Loop over combinations of j, n, and method within each worker
  for (j in JJ) {
    for (n in NN) {
      XX_data <- XX(j, n) # Generate the observations

      for (method in MM) {
        start_time <- Sys.time() # Start the timer (***)
        ISE_value <- ISE(XX_data, j, method) # Calculate ISE for the current replication
        end_time <- Sys.time() # End the timer (***)
        elapsed_time_seconds <- as.numeric(difftime(end_time, start_time, units = "secs")) # (***)
        
        # Store the result for this specific replication
        local_raw_results <- rbind(
          local_raw_results,
          data.frame(
            j = j,
            n = n,
            method = method,
            ISE = ISE_value,
            time_seconds = elapsed_time_seconds, # Save elapsed time (***)
            stringsAsFactors = FALSE
          )
        )
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

# Save the raw results to a CSV file in the specified path
raw_output_file <- file.path(path, "raw_ISE_results_iid.csv")
write.csv(raw_results, raw_output_file, row.names = FALSE)

print("Raw results saved to raw_ISE_results_iid.csv")

#########################
## Process the results ##
#########################

# Create a data frame to store the summary results
summary_results <- data.frame(
  j = integer(),
  n = integer(),
  method = character(),
  mean_ISE = numeric(),
  sd_ISE = numeric(),
  median_ISE = numeric(),
  IQR_ISE = numeric(),
  mean_time_seconds = numeric(), # New column for mean elapsed time (***)
  stringsAsFactors = FALSE
)

# Loop through the results to compute the summary statistics
for (j in JJ) {
  for (n in NN) {
    for (method in MM) {
      # Filter the raw results by j, n, and method
      filtered_results <- raw_results %>%
        filter(j == !!j, n == !!n, method == !!method)
      
      ISE_values <- filtered_results$ISE
      elapsed_times <- filtered_results$time_seconds # Extract elapsed times (***)
      
      mean_ISE <- mean(ISE_values)
      sd_ISE <- sd(ISE_values)
      median_ISE <- median(ISE_values)
      IQR_ISE <- IQR(ISE_values)
      mean_time_seconds <- mean(elapsed_times) # Calculate mean elapsed time (***)
      
      # Store the summary results
      summary_results <- rbind(
        summary_results,
        data.frame(
          j = j,
          n = n,
          method = method,
          mean_ISE = mean_ISE,
          sd_ISE = sd_ISE,
          median_ISE = median_ISE,
          IQR_ISE = IQR_ISE,
          mean_time_seconds = mean_time_seconds, # Save mean elapsed time (***)
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

# Save the summary results to a CSV file in the specified path
summary_output_file <- file.path(path, "ISE_results_iid.csv")
write.csv(summary_results, summary_output_file, row.names = FALSE)

print("Summary results saved to ISE_results_iid.csv")
