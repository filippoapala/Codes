#THIS CODE AIMS AT GENERATING I.I.D. BOOTSTRAP SAMPLES AND THEN ESTIMATING THE MODEL OVER EACH SAMPLE IN ORDER TO BOOTSTRAP THE IMPULSE RESPONSE FUNCTIONS
#IT CONTAINS BOTH THE CODE FOR BOTH A RESIDUAL I.I.D. AND MOVING BLOCK BOOTSTRAP


##GENERATE BOOTSTRAP SAMPLE WITH MOVING BLOCK BOOTSTRAP  
generate_moving_block_sample = function(var_model,ts_data, n_lags){
  
 Residuals = residuals(var_model)
 n_lags = var_model$p
 n_vars = ncol(Residuals)
 coef_matrix = var_coefficients(var_model,ts_data)
 
  # Step 3: Define block length and create bootstrap residuals
  block_length = 20  # Adjust block length (ℓ) as needed
  T = nrow(Residuals)  # Number of time points
  N = ceiling(T / block_length)  # Number of blocks
  bootstrap_residuals = vector("list", N)
  
  # Create blocks and perform block sampling
  for (i in 1:N) {
    block_start = sample(1:(T - block_length), 1)
    bootstrap_residuals[[i]] = Residuals[block_start:(block_start + block_length - 1), ]
  }
  
  # Combine blocks to form the bootstrap residuals
  bootstrap_residuals_combined = do.call(rbind, bootstrap_residuals)
  if (length(bootstrap_residuals_combined) > T) {
    bootstrap_residuals_combined = bootstrap_residuals_combined[1:T, ]
  }
  
  print(nrow(bootstrap_residuals_combined))
  
  # Step 4: Center the bootstrap residuals
  bootstrap_residuals_centered = bootstrap_residuals_combined - colMeans(bootstrap_residuals_combined)
  
  # Step 5: Generate the bootstrap sample
  y_star <- matrix(0, nrow = T, ncol = ncol(Residuals))  # Placeholder for bootstrap sample
  
  # Initialize the bootstrap sample with sample means for the first p observations
  y_star[1 : n_lags, ] =  colMeans(Residuals) # Set initial values as data sample mean
  
  
  for (t in (n_lags + 1):(T)) {
    # Start with the constant term (intercept) from the last column of A_full
    # Extracts the intercept for each variabl
    const = coef_matrix[,n_lags * n_vars +1]
    
    #simulated_data[t, ] = A_full[, n_vars] 
    new_value =  numeric(n_vars)
    
    # Add contributions from lagged values
    for (lag in 1:n_lags) {
      # Get the lagged data as a column vector
      lagged_values = as.numeric(y_star[t - lag, ])
      
      # Extract the corresponding block of coefficients for the current lag
      coef_block = coef_matrix[, ((lag - 1) * n_vars + 1):(lag * n_vars)]
      
      
      # Update `new_value` by adding the matrix product
      new_value = new_value + coef_block %*% lagged_values
    }
    # Add the resampled residual
    y_star[t, ] = t(new_value) + bootstrap_residuals_centered[t-n_lags,] + t(const)
  }
  
  
 return(y_star)
  
}




#GENERATE BOOTSTRAP SAMPLE THROUGH I.I.D. RESIDUAL BOOTSTRAP 
generate_bootstrap_sample = function(var_model,residuals,ts_data, n_lags) {
  # Get the number of observations and variables
  n_obs = nrow(ts_data)
  n_vars = ncol(ts_data)
  
  # Initialize the new bootstrap sample using the original data for initial conditions
  simulated_data = matrix(NA, nrow = n_obs, ncol = n_vars)
  column_means = colMeans(ts_data)
  
  for (i in seq_along(ts_data)) {
    simulated_data[1:n_lags, i] = column_means[i]
  }
  
  # init = as.matrix(ts_data[1:n_lags, ])
  # simulated_data[1:n_lags, ] = init
  
  # Resample the residuals with replacement
  boot_residuals = residuals[sample(1:nrow(residuals), size = nrow(residuals), replace = TRUE), ]
  boot_residuals = boot_residuals - mean(boot_residuals)
  #boot_residuals = boot_residuals - mean(boot_residuals)
  # Coefficients from the VAR model
  coef_matrix = var_coefficients(var_model,ts_data)
  
  
  for (t in (n_lags + 1):(n_obs)) {
    # Start with the constant term (intercept) from the last column of A_full
    # Extracts the intercept for each variabl
    const = coef_matrix[,n_lags * n_vars +1]
    
    #simulated_data[t, ] = A_full[, n_vars] 
    new_value =  numeric(n_vars)
    
    # Add contributions from lagged values
    for (lag in 1:n_lags) {
      # Get the lagged data as a column vector
      lagged_values = as.numeric(simulated_data[t - lag, ])
      
      # Extract the corresponding block of coefficients for the current lag
      coef_block = coef_matrix[, ((lag - 1) * M + 1):(lag * M)]
      
     
      # Update `new_value` by adding the matrix product
      new_value = new_value + coef_block %*% lagged_values
    }
    # Add the resampled residual
    simulated_data[t, ] =  t(new_value) + boot_residuals[t-n_lags,] + t(const)
  }
  
  # Convert to time series format
  
  return(simulated_data)
  
}




#GENERATE 1000 BOOTSTRAP IRF'S FOR A CHOLEKSY SVAR 
generateBootIRFs = function(n_bootstrap, H, pi){
  
  # Assuming you have the necessary function to estimate the SVAR model, e.g., SVAR()
  n_bootstrap = 1000 #SET NUMBER OF BOOTSTRAP REPETITIONS
  H = 15 #SET IRF'S HORIZON 
  pi = 2 #SET THE NUMBER OF LAGS
  
  # Storage for VAR and SVAR bootstrap results
  var_boot_results = vector("list", n_bootstrap)
  svar_boot_IRF = vector("list", n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    
    # Generate bootstrap sample with moving block bootstrap
    boot_sample = generate_moving_block_sample(var_model2T, data_neriCompleto, pi) #you can change with i.i.d. residual bootstrap by changing name and arguments of the function
    
    # Refit the VAR model on the bootstrap sample
    boot_var_model = vars::VAR(boot_sample, p = pi, type = "const")
    
    B_boot = t(chol(cov(residuals(boot_var_model)))) #extract the cholesky decomposition of the residiual's covariance matrix
    Boot_shocks = residuals(boot_var_model) %*% solve(B_boot) #invert the structural matrix to get the shocks
    B_boot[,1] = B_boot[,1] *  sd(Boot_shocks[,1]) #normalizes the structural matrix to preserve the unit variance of the shocks
    B_boot[,2] = B_boot[,2] *  sd(Boot_shocks[,2])
    B_boot[,3] = B_boot[,3] *  sd(Boot_shocks[,3])

    #extract the coefficients matrix of the VAR model estimated on the boostrap sample
    A_matrix = var_coefficients(boot_var_model,data_neriCompleto)
    
    # Estimate the SVAR model on the refitted VAR model (using your SVAR methodology)
    svar_boot_IRF[[i]] = compute_irf_H(A_matrix,B_boot,H)
  }
  
  
  
  #GENERATE BOOTSTRAP IRF'S USING FAST ICA FOR STATISTICAL IDENTIFICATION OF THE STRUCTURAL MATRIX, USES SVARS PACKAGE FOR ICA PROCEDURE
  #YOU CAN USE DIFFERENT METHODS IN ORDER TO GET THE STRUCUTRAL MATRIX, BOTH FROM THE SVARS PACKAGES AND THE STUDENT-T ML ALGORITHM IN THE OTHER FILE
  
  #ICA bootstrap
  n_bootstrap = 1000
  H = 15
  pi = 2
  #DEFINE OBEJECTS IN ORDER TO STORE RESULTS
  var_boot_results = vector("list", n_bootstrap)
  ICAsvar_boot_IRF = vector("list", n_bootstrap)
  ICAboot_B = vector("list", n_bootstrap)
  B_sample = id.dc(var_model2T,PIT = FALSE) #replace the first input function with your VAR model object
  eigen_u = eigen(cov(residuals_var2T)) #eigendecomposition of the covariance matrix
  eigen_u_st = diag((eigen_u$values)^(1/2))
 
  
  for (i in 1:n_bootstrap) {
    
    # Generate bootstrap sample
    boot_sample = generate_moving_block_sample(var_model2T,data_neriCompleto, pi)  # Generate bootstrap sample, you can change it
    
    # Refit the VAR model on the bootstrap sample
    boot_var_model = VAR(boot_sample, p = pi, type = "const")
    boot_var_residuals = residuals(boot_var_model)
    eigen_u_star = eigen(cov(boot_var_residuals))
    eigen_u_star_st = (eigen_u_star$values)^(-1/2)
    eigen_u_star_st = diag(eigen_u_star_st)
    
    boot_ICA = id.dc(boot_var_model,PIT = FALSE) #replace with your ICA estimator, here I used distance covariance from the svars package
    B_boot_st = boot_ICA$B #extract structural matrix 
   
    B_bootShocks = boot_var_residuals %*%  solve(B_boot_st) #extract the shocks
    B_boot_st[,1] = B_boot_st[,1] * sd(B_bootShocks[,1]) 
    B_boot_st[,2] = B_boot_st[,2] * sd(B_bootShocks[,2])
    B_boot_st[,3] = B_boot_st[,3] * sd(B_bootShocks[,3])
  
    
    #SELECT BEST PERMUTATION AND SIGN GIVEN THE WELL SHAPED MATRIX AND THE BOOTSTRAP ESTIMATES
    B_boot = select_best_permutation_and_sign(B_sample$B, B_boot_st, eigen_u_st) #replace the first input function with your sample B matrix estimates
    
    #extract VAR coefficients matrix from the VAR boot model, replace the second input element with your dataset
    A_matrix = var_coefficients(boot_var_model,data_neriCompleto)
        
    ICAsvar_boot_IRF[[i]] = compute_irf_H(A_matrix,B_boot,H) #store IRF's and B^** matrix
    ICAboot_B[[i]] = B_boot
  
  }
 

}













