#MAXIMUM LIKELIHOOD ESTIMATOR OF THE STRUCTRAL MATRIX USING STUDENT-T DENSITIES
#1)density of a STUDENT-T WITH UNIT SCALE
#2)sum of the likelihoods of the univariate sudent-t
#3) optimize over lambda and beta to get structural parameters and degrees of freedom


#1)
#DENSITY OF AN UNIVARIATE STUDENT-T DISTRIBUTION WITH SCALE 1
tdensity_i = function(lambda, sigma, l, U, t, i, B) {
  
  # Small constant to avoid zero or negative values
  epsilon = 1e-10
  
  # Ensure sigma and lambda are positive
  sigma = sigma + epsilon
  lambda = lambda + epsilon
  
  # Compute the log-likelihood components
  term1 = lgamma((lambda[i] + 1) / 2)
  term2 = lgamma(lambda[i] / 2)
  
  # Constant term
  term3 = -0.5 * log(lambda[i] * pi)
  
  
  # Log density computation
  term4 = -((lambda[i] + 1) / 2) * log(1 + ((l[[i]] %*% t(U[t, ] %*% solve(B)))^2) / lambda[i])
  
  # Total log-density
  log_f = term1 - term2 + term3 + term4
  
  return(log_f)
}


#2)
#SUM OF THE LOG-LIKELIHOOD STUDENT-T
log_like = function(theta, U, l, B4) {
  
  # Reshape B from theta based on the dimensions of B4
  B = matrix(theta[1:9], nrow = nrow(B4), ncol = ncol(B4))

  
  # Extract lambda from theta
  lambda = theta[10:12]
  
  # Fixed sigma values
  sigma = c(1, 1, 1)
  
  # Initialize log-likelihood
  lt = 0
  
  # Loop over data points
  for (t in 1:nrow(U)) {
    
    # Initialize sum of log f_i for each observation
    sum_log_fi = 0
    
    for (i in 1:length(lambda)) {
      # Calculate density for each observation
      sum_log_fi = sum_log_fi + tdensity_i(lambda, sigma, l, U, t, i, B)
    }
    
    # Adjust with determinant of B and sigma terms
    sum_log_fi = sum_log_fi - log(abs(det(B))) - sum(log(sigma))
    
    # Add the sum to the log-likelihood
    lt = lt + sum_log_fi
  }
  
  # Average the log-likelihood over the number of observations
  lt = -lt / nrow(U)
  
  print(lt)
  return(lt)
}





#3
#OPTIMIZATION OF THE LOG-LIKELIHOOD WITH RESPCT TO BETA AND LAMBDA
mle_theta = function(U, l, B0, random_lambda, random_sigma) {
  
  # Initial beta as a vector from B0 matrix
  initial_beta = as.vector(B0)
  
  # Initial sigma and lambda
  initial_sigma = random_sigma
  initial_lambda = random_lambda
  
  # Combine initial parameters into theta vector
  initial_theta = c(initial_beta, initial_lambda)
  
  # Optimization settings
  options = list(maxit = 20000, reltol = 10^-10)
  
  # Define lower and upper bounds
  lb = c(rep(-Inf, length(initial_beta)), rep(0, length(initial_sigma)), rep(0, length(initial_lambda)))
  ub = c(rep(Inf, length(initial_beta)), rep(Inf, length(initial_sigma)), rep(Inf, length(initial_lambda)))
  
  # Define the function for optimization
  my_fun4 = function(theta) log_like(theta, U, l, B0)
  
  # Run the optimization using optim
  result = optim(par = initial_theta, fn = my_fun4, method = "L-BFGS", lower = lb, upper = ub, control = options, hessian = TRUE)
  
  # Extract optimized values
  optimal_theta = result$par
  fval = result$value
  exitflag = result$convergence
  output = result$message
  grad = result$gradient
  hessian = result$hessian
  
  # Reshape optimal_beta
  optimal_beta = matrix(optimal_theta[1:9], nrow = nrow(B0), ncol = ncol(B0))
  
  # Extract lambda and sigma
  optimal_lambda = optimal_theta[10:12]
  optimal_sigma = c(1, 1, 1)
  
  return(list(optimal_beta = optimal_beta, optimal_sigma = optimal_sigma, optimal_lambda = optimal_lambda,
              fval = fval, exitflag = exitflag, output = output, grad = grad, hessian = hessian))
}





# Initialize lists to store results
Bthetas4 = list()
best_sigmas = list()
best_lambdas = list()
gradients = list()
hessians = list()
SE_Hessian_matrix = list()


# Initialize an empty list to store unit vectors
l = list()

# Define the number of unit vectors to create
M = 3  # Replace with desired value of M

# Loop to create each unit vector
for (i in 1:M) {
  # Initialize a zero vector of length M
  unit_vector = rep(0, M)
  
  # Set the i-th element to 1
  unit_vector[i] = 1
  
  # Store the unit vector in the list
  l[[i]] = unit_vector
}





lower_bound = -1
upper_bound = 1
# TRY DIFFERENT STARTING POINTS FOR GLOBAL OPTIMIZER
for (i in 1:15) {
  random_lambda = 1 + (10 - 1) * runif(3)
  random_sigma = 1 + (2 - 1) * runif(3)
  random_B = lower_bound + (upper_bound - lower_bound) * matrix(runif(9), nrow = 3, ncol = 3)
  
  # Run MLE optimization
  mle_result = mle_theta(residuals_var3, l, random_B, random_lambda, random_sigma)
  
  # Store results
  Bthetas4[[i]] = mle_result$optimal_beta
  best_sigmas[[i]] = mle_result$optimal_sigma
  best_lambdas[[i]] = mle_result$optimal_lambda
  gradients[[i]] = mle_result$grad
  hessians[[i]] = mle_result$hessian
  
  # Calculate standard errors from Hessian
  SE_Hessian_matrix[[i]] = sqrt(diag(solve(hessians[[i]])))
}

