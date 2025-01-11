#THIS SCRIPT CONTAINS A SET OF FUNCTIONS THAT ARE USED THROUGH THE ANALYSIS
# Load required library
library(gtools)
library(combinat) 


# -------------------------COMPUTE CHI STATISTICS FOR TESTING THE JOINT HYPOTHESIS OF ML ESTIMATES---------------------------

chi_statML = function(restricted_log, unrestricted_log){
  
  chi = 2(unrestricted_log - restricted_log)
  
  return(chi)
}


# -------------------------COMPUTE CHI STATISTICS FOR TESTING THE JOINT HYPOTHESIS OF BOOTSTRAP MATRICES---------------------------

chi_statBoot = function(boot_B, B_hat){
  
  vectorized_matrices = sapply(boot_B, as.vector) #vectorize matrices
  vectorized_matrices = t(vectorized_matrices)
  vectorized_matrices =  vectorized_matrices[,c(4,7,8)] #set the desired positions
  cov_B_star_vec = cov(vectorized_matrices) #computes covariance of the bootstrap results
  R = matrix(0, nrow = 3, ncol = 3^2) #set the desired dimension
  
  # Set the appropriate positions to 1 to impose restrictions
  R[1,4] =  1  # Restriction
  R[2,7] = 1  # Restriction
  R[3,8] = 1  # Restriction
  
  r = numeric(3)
  r = t(matrix(r, nrow  = 1))
  chi = t((R %*% as.vector(B_hat) - r)) %*% solve(cov_B_star_vec) %*% (R %*% as.vector(B_hat) - r)
  
  return(chi)
  
}


# -------------------------EXTRACT THE VAR COEFFICIENTS FROM A VAR MODEL OBJECT---------------------------

var_coefficients = function(var_model, ts_data) {
  

  M = ncol(ts_data)  # Number of variables in the VAR model
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



# ------------------------- USW THE FROBENIOUS NORM TO GET THE CLOSEST VERSION TO B_HAT---------------------------
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




# ------------------------- USE A SIMILARITY CRITERION TO GET THE CLOSEST VERSION TO B_HAT---------------------------

select_best_permutation_and_sign = function(D_p, D_r, Omega) {

  
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





# ------------------------- GET THE COLUMN PERMUTATION THAT MAXIMIZES THE ABSOLUTE VALUE OF THE DIAGONAL---------------------------

maximize_absolute_product_columns = function(B) {
  
  # Number of columns (assumes B is square)
  n = ncol(B)
  
  # Generate all permutations of column indices
  column_perms = permutations(n = n, r = n, v = 1:n)
  
  
  # Initialize variables to store the best result
  max_prod = -Inf
  best_B = NULL
  
  # Iterate over all permutations of columns
  for (i in 1:nrow(column_perms)) {
    
    # Permute the columns
    permuted_B = B[, column_perms[i, ]]
    
    
    # Extract the diagonal elements (row i, column i after reordering)
    diagonal = sapply(1:n, function(k) permuted_B[k, k])
    
    # Calculate the product
    current_prod = prod(diag(permuted_B))
    
    # Update if the absolute product is greater
    if (abs(current_prod) > max_prod) {
      max_prod = abs(current_prod)
      best_B = permuted_B
      
    }
  }
  
  
  # Ensure all diagonal elements are positive
  for (i in 1:ncol(best_B)) {
    if (best_B[i, i] < 0) {
      best_B[, i] = -best_B[, i]
    }
  }
  
  return(best_B)
}




# ------------------------- CREATE THE UNIT VECTOR FOR ML ESTIMATION---------------------------

Unit_vector_gen = function(M){
  
  # Initialize an empty list to store unit vectors
  l = list()
  
  
  # Loop to create each unit vector
  for (i in 1:M) {
    # Initialize a zero vector of length M
    unit_vector = t(rep(0, M))
    
    # Set the i-th element to 1
    unit_vector[i] = 1
    
    # Store the unit vector in the list
    l[[i]] = matrix(unit_vector, ncol = 1)
    
    
  }
  
  return(l)
  
}

