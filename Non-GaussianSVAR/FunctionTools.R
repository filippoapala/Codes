#THIS SCRIPT CONTAINS A SET OF FUNCTIONS THAT ARE USED THROUGH THE ANALYSIS


#COMPUTE CHI STATISTICS FOR TESTING THE JOINT HYPOTHESIS OF BOOTSTRAP MATRICES
chi_stat = function(boot_B, B_hat){

  vectorized_matrices = sapply(FAboot_B, as.vector)#vectorize matrices
  vectorized_matrices = t(vectorized_matrices)
  vectorized_matrices =  vectorized_matrices[,c(4,7,8)] #set the desired positions
  cov_B_star_vec = cov(vectorized_matrices) #computes covariance of the bootstrap results
  R = matrix(0, nrow = 3, ncol = 3^2) #set the desired dimension

  # Set the appropriate positions to 1 to impose restrictions
  R[1,4] =  1  # Restriction on b_12
  R[2,7] = 1  # Restriction on b_13
  R[3,8] = 1  # Restriction on b_23

  r = numeric(3)
  r = t(matrix(r, nrow  = 1))

  chi = t((R %*% as.vector(B_hat) - r)) %*% solve(cov_B_star_vec) %*% (R %*% as.vector(B_hat) - r)


}


#COMPUTE THE ILMONEN INDEX FOR MONTE CARLO SIMULATIONS, it computes a mesure of similarity between the original matrix and the proposed one,
computeIlmonenIndex = function(B_true, B_est) {
 
  # Ensure matrices are square and have the same dimensions
  k_true = nrow(B_true)
  n_true = ncol(B_true)
  k_est = nrow(B_est)
  n_est = ncol(B_est)
  
  k = k_true  # Size of the square matrix
  
  # Define the objective function for minimization: Frobenius norm calculation
  objectiveFn = function(C) {
    # Convert the vectorized form back to matrix form
    C_matrix = matrix(C, nrow = k, ncol = k)
    
    # Frobenius norm || C * inv(B_est) * B_true - I_k ||
    distance_matrix = C_matrix %*% solve(B_est) %*% B_true - diag(k)
    d = norm(distance_matrix, type = "F")  # Frobenius norm
    return(d)
  }
  
  # Construct initial guess for C as an identity matrix
  C_init = diag(k)
  
  # Use optimization to find the best C
  C_opt = optim(par = as.vector(C_init), fn = objectiveFn, method = "Nelder-Mead")$par
  
  # Reshape optimized vector back to matrix form
  C_opt_matrix = matrix(C_opt, nrow = k, ncol = k)
  
  # Compute the final Ilmonen index
  distance = norm(C_opt_matrix %*% solve(B_est) %*% B_true - diag(k), type = "F")
  ilmonen_index = distance / sqrt(k - 1)
  
  return(ilmonen_index)
}


#EXTRACT THE VAR COEFFICIENTS FROM A VAR MODEL OBJECT
var_coefficients = function(var_model, ts_data) {
  
  # Number of variables (M) and lags (p) used in the VAR model
  M = ncol(ts_data)  # Assuming M = number of variables in the VAR model
  p = var_model$p   # Number of lags
  
  # Initialize matrix to store VAR coefficients in the required format
  A_full = matrix(0, nrow = M, ncol = (M * p) + 1)  # +1 for the constant term
  
  
  # Extract the coefficients of the VAR model
  coefficients = as.matrix(coef(var_model))
  
  # Place the lagged coefficients in the appropriate blocks
  for (m in 1:M) {
    # Each block of lagged coefficients is extracted from the model
    coef_matrix = coefficients[[m]]
    equation = t(coef_matrix[,1])
    A_full[m,] = equation
    
  }
  return((A_full))
}


library(combinat) # for permutations function

# Function to align bootstrap replication to original estimate using Frobenius norm
align_bootstrap_sample = function(B_hat, B_boot) {
  n = ncol(B_hat)
  perms = permn(1:n)  # Generate all permutations of columns
  sign_combinations = expand.grid(rep(list(c(-1, 1)), n))  # All combinations of sign flips
  
  # Initialize minimum norm
  min_norm = Inf
  best_B_boot_aligned = B_boot
  
  # Loop over each permutation and sign combination
  for (perm in perms) {
    for (signs in 1:nrow(sign_combinations)) {
      # Apply permutation and signs
      B_boot_perm = B_boot[, perm] * as.numeric(sign_combinations[signs, ])
      
      # Compute Frobenius norm
      frob_norm = norm(B_hat - B_boot_perm, type = "F")
      
      # Update if this configuration has the smallest norm
      if (frob_norm < min_norm) {
        min_norm = frob_norm
        best_B_boot_aligned = B_boot_perm
      }
    }
  }
  
  return(best_B_boot_aligned)
}



library(gtools)  # for generating permutations
#SELECT THE BEST PERMUTATION AND SIGN FOR A MATRIX WITH RESPECT TO A WELL SHAPED OBJECT

select_best_permutation_and_sign <- function(D_p, D_r, Omega) {
  # Ensure D_p, D_r, and Omega are square matrices with the same dimensions
  #if (!all(dim(D_p) == dim(D_r) && dim(D_p) == dim(Omega))) {
   # stop("D_p, D_r, and Omega must have the same dimensions")
  #}
  
  # Get the dimension of the matrices
  K = nrow(D_p)
  
  # Extract diagonal elements of Omega for rescaling
  omega_diag = diag(Omega)
  
  # Generate all column permutations
  permutations = permutations(K, K)
  
  # Generate all combinations of sign flips for K columns
  sign_combinations = expand.grid(rep(list(c(1, -1)), K))
  sign_combinations = as.matrix(sign_combinations)  # Ensure sign_combinations is numeric
  
  # Initialize variables to track the best matrix and minimum criterion value
  min_criterion_value = Inf
  best_matrix = NULL
  
  # Loop through each permutation of columns
  for (perm in 1:nrow(permutations)) {
    # Apply the permutation to the columns of D_r
    permuted_D_r = D_r[, permutations[perm, ]]
    
    # Loop through each combination of sign flips
    for (signs in 1:nrow(sign_combinations)) {
      # Apply sign flips to the permuted columns
      signed_D_r = sweep(permuted_D_r, 2, sign_combinations[signs, ], `*`)
      
      # Calculate the criterion value for the signed and permuted D_r
      criterion_value = 0
      for (i in 1:K) {
        for (j in 1:K) {
          # Calculate squared difference and apply rescaling
          diff = D_p[i, j] - signed_D_r[i, j]
          scaled_diff = (diff^2) / omega_diag[i]
          
          # Apply the indicator function to check for opposite signs
          if (D_p[i, j] * signed_D_r[i, j] < 0) {
            criterion_value = criterion_value + scaled_diff
          }
        }
      }
      
      # Update best_matrix if the current matrix has a lower criterion value
      if (criterion_value < min_criterion_value) {
        min_criterion_value = criterion_value
        best_matrix = signed_D_r
      }
    }
  }
  
  # Return the best matrix that minimizes the criterion
  return(best_matrix)
}




