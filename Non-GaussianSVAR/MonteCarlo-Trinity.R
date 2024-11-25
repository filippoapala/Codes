#THIS SCRIPT SIMULATES DATA FROM A 3 VARIABLES MONETARY POLICY SVAR(2) MODEL WITH STUDENT-T SHOCKS
#HERE THE PARAMETERS ARE THOUGH TO MIMIC A CALIBRATED SVAR MODEL IN THE LITERATURE, BUT JUST LIKE THE STUDENT-T SHOCKS CAN BE EASILY CHANGED


simulate_data = function(n) {
  
  # Structural matrix B
  D = matrix(c(2.32, -0.48, -0.41,
                      0.72, 2.32, -0.22,
                      0.98, 1.57, 0.76), nrow = 3, byrow = TRUE)
  
  
  phi =  matrix(c(0.74, -0.09, -0.16,
                        0.13, 0.44, -0.06,
                        0.24, 0.30, 0.53), nrow = 3, byrow = TRUE)
  
  Q =  matrix(c(0.5, 0, 0,
                  0, 0.5, 0,
                  0, 0, 0.5), nrow = 3, byrow = TRUE)
  
  
  # Coefficients for Y(t-1)
  A_1 = phi + D %*%Q%*%solve(D)
  
  # Coefficients for Y(t-2)
  A_2 = -(D %*%Q%*%solve(D))%*%phi
  
  
  
  # Simulation parameters
  numObservations = n  # Set the number of observations to 'n'
  numVars = 3          # Number of variables in the VAR model
  numLags = 2          # Number of lags
  
  # Degrees of freedom and scale parameters for Student's t distribution
  nu = c(12.3,8, 5)      # Degrees of freedom for each variable.  #first simulation (2,5,8) scales always 1
  scales = c(1, 1, 1)  # Scale parameters for each variable
  constant = c(6, 15, 11)  # Constant term for each variable
  
  # Load necessary package for Student's t-distribution
  if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
  library(MASS)
  
  # Generate Student's t-distributed shocks
 # set.seed(123)  # Set seed for reproducibility
  e = matrix(0, nrow = numObservations, ncol = numVars)
  
  for (i in 1:numVars) {
    e[, i] = scales[i] * rt(numObservations, df = nu[i])
  }
  
  # Calculate innovations by multiplying the shocks with the transpose of B matrix
  innovations = e %*% t(D)
  
  # Pre-allocate space for simulated data
  Y = matrix(0, nrow = numObservations, ncol = numVars)
  
  # Set initial values for the first numLags observations
  Y[1:numLags, ] = innovations[1:numLags, ] + matrix(rep(constant, each = numLags), nrow = numLags, byrow = TRUE)
  
  # Generate the time series data
  for (t in (numLags + 1):numObservations) {
    # Calculate lagged terms
    laggedTerms = A_1 %*% Y[t - 1, ] + A_2 %*% Y[t - 2, ]
    
    # Simulate current period values
    Y[t, ] = laggedTerms + innovations[t, ] + constant
  }
  
  return(Y)
}



