


# ------------------------- DENSITY OF AN UNIVARIATE STUDENT-T DISTRIBUTION WITH SCALE 1 ---------------------------

tdensity_i = function(lambda, sigma, l, U, t, i, B) {
  
  
  epsilon = 1e-10 #Small constant to avoid zero or negative values
  
  #sigma = exp(sigma)
  lambda = exp(lambda)
  sigma = sigma + epsilon  #Ensure sigma and lambda are positive
  lambda = lambda + epsilon
  
  #Compute the log-likelihood components
  term1 = lgamma((lambda[i] + 1) / 2)
  term2 = lgamma(lambda[i] / 2)
  
  #Constant term
  term3 = -0.5 * log(lambda[i] * pi)
  
  #Log density computation
  term4 = -((lambda[i] + 1) / 2) * log(1 + ((( t(l[[i]])  %*% solve(B)) %*% (U[t, ]))^2) /  lambda[i])
  
  #Total log-density
  log_f = term1 - term2 + term3 + term4
  
  return(log_f)
}

# ------------------------- SUM OF THE LOG-LIKELIHOOD STUDENT-T ---------------------------

log_like = function(theta, U, l, B4) {
  
  #Reshape B from theta based on the dimensions of B4
  B = matrix(theta[1:9], nrow = nrow(B4), ncol = ncol(B4))
  
  #Extract lambda from theta
  lambda = theta[10:12]
  
  # Fixed sigma values
  sigma = c(1, 1, 1)
  
  #Initialize log-likelihood
  lt = 0
  #Loop over data points
  for (t in 1:nrow(U)) {
    
    #Initialize sum of log f_i for each observation
    sum_log_fi = 0
    
    for (i in 1:length(lambda)) {
      sum_log_fi = sum_log_fi + tdensity_i(lambda, sigma, l, U, t, i, B) #Calculate density for each observation
    }
    
    sum_log_fi = sum_log_fi - log(abs(det(B))) - sum(log(sigma))  #Adjust with determinant of B and sigma terms
    
    lt = lt + sum_log_fi #Add the sum to the log-likelihood
  }
  
  lt = -lt / nrow(U) #Average the log-likelihood over the number of observations
  
  print(lt)
  return(lt)
}





# ------------------------- OPTIMIZATION OF THE LOG-LIKELIHOOD WITH RESPECT TO BETA AND LAMBDA ---------------------------

mle_theta = function(U, l, B0, random_lambda, random_sigma) {
  
  #Initial beta as a vector from B0 matrix
  initial_beta = as.vector(B0)
  
  #Initial sigma and lambda
  initial_lambda = log(random_lambda)
  initial_sigma = random_sigma
  
  initial_theta = c(initial_beta, initial_lambda) #Combine initial parameters into theta vector
  
  options = list(maxit = 20000, reltol = 10^-15) #Optimization settings
  
  # Define lower and upper bounds
  #lb = c(rep(-Inf, length(initial_beta)), rep(0, length(initial_sigma)), rep(0, length(initial_lambda)))
  #ub = c(rep(Inf, length(initial_beta)), rep(Inf, length(initial_sigma)), rep(Inf, length(initial_lambda)))
  
  # Define the function for optimization
  my_fun4 = function(theta) log_like(theta, U, l, B0)
  
  # Run the optimization using optim
  #result = optim(par = initial_theta, fn = my_fun4, method = "L-BFGS", lower = lb, upper = ub, control = options, hessian = TRUE) #Constrained optimization
  result = optim(par = initial_theta, fn = my_fun4, method = "L-BFGS", control = options, hessian = TRUE) #Unconstrained optimization
  
  #Extract optimized values
  optimal_theta = result$par
  fval = result$value
  exitflag = result$convergence
  output = result$message
  grad = result$gradient
  hessian = result$hessian
  
  optimal_beta = matrix(optimal_theta[1:9], nrow = nrow(B0), ncol = ncol(B0), byrow = FALSE) #Reshape optimal_beta
  optimal_beta = maximize_absolute_product_columns(optimal_beta) #Get the permutations of the columns that maximizes the absolute product of the diagonal elements
  
  # Extract lambda and sigma
  optimal_lambda = optimal_theta[10:12]
  optimal_lambda = exp(optimal_lambda)
  optimal_sigma = c(1,1,1)
  
  return(list(optimal_beta = optimal_beta, optimal_sigma = optimal_sigma, optimal_lambda = optimal_lambda,
              fval = fval, exitflag = exitflag, output = output, grad = grad, hessian = hessian))
}





