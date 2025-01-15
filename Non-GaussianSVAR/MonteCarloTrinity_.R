
library(MASS) #package for Student-t distribution

# ------------------------- SIMULATE DATA FROM A TRIVARIATE SVAR MODEL WITH STUDENT-SHOCKS  ---------------------------

simulate_data = function(n) {
  
  # Structural matrix B
  D = matrix(c(0.932, 0.1, -0.2,
               -0.1, 0.47, 0.26,
               0.18, 0.057, 0.76), nrow = 3, byrow = TRUE)
  
  
  phi =  matrix(c(0.74, -0.09, -0.16,
                  0.13, 0.44, -0.06,
                  0.24, 0.30, 0.53), nrow = 3, byrow = TRUE)
  
  Q =  matrix(c(0.5, 0, 0,
                0, 0.5, 0,
                0, 0, 0.5), nrow = 3, byrow = TRUE)
  
  
  
  A_1 = phi + D %*%Q%*%solve(D)  #Coefficients for Y(t-1)
  
  #1.24 -0.09 -0.16
  #0.13  0.94 -0.06
  #0.24  0.30  1.03
  

  A_2 = -(D %*%Q%*%solve(D))%*%phi #Coefficients for Y(t-2)
  
  #-0.370  0.045  0.080
  #-0.065 -0.220  0.030
  #-0.120 -0.150 -0.265
  
  #Simulation parameters
  numObservations = n  #Set the number of observations to 'n'
  numVars = 3          #Number of variables in the VAR model
  numLags = 2          #Number of lags
  
  #Degrees of freedom and scale parameters for Student's t distribution
  nu = c(4,5, 3)      #Degrees of freedom for each variable
  scales = c(1, 1, 1)  #Scale parameters for each variable
  constant = c(6, 15, 11)  #Constant term for each variable
  
  
  
  e = matrix(0, nrow = numObservations, ncol = numVars) #Generate Student's t-distributed shocks
  
  for (i in 1:numVars) {
    e[, i] = scales[i] * rt(numObservations, df = nu[i])
  }
  
  #Calculate innovations by multiplying the shocks with the transpose of B matrix
  innovations = e %*% t(D)
  
  #Pre-allocate space for simulated data
  Y = matrix(0, nrow = numObservations, ncol = numVars)
  
  #Set initial values for the first numLags observations
  Y[1:numLags, ] = innovations[1:numLags, ] + matrix(rep(constant, each = numLags), nrow = numLags, byrow = TRUE)
  
  #Generate the time series data
  for (t in (numLags + 1): numObservations) {
    #Calculate lagged terms
    laggedTerms = A_1 %*% Y[t - 1, ] + A_2 %*% Y[t - 2, ]
    
    #Simulate current period values
    Y[t, ] = laggedTerms + innovations[t, ] + constant
  }
  
  return(Y)
}




# ------------------------- TEST THE MODELS ON THE SIMULATED DATA  ---------------------------


#Simulate 1000 datasets of 100 observations each
YMC = vector("list", 1000)
for (i  in 1:1000){ 
  YMC[[i]] = simulate_data(100)
  
}  

STvar = StVAR(YMC[[400]],Trend=1,lag=2,v=5,maxiter=1000,meth="BFGS",hes="TRUE",init="na") # Student-t VAR model estimated with ML

varY = VAR(YMC[[35]], p = 2, type = "const") #Select a sample 
ia = id.dc(varY, PIT = FALSE).       #ICA by distance covariance
b = ia$B

B_bootShocks = residuals(varY)%*%  solve(b) #Extract the shocks
b[,1] = b[,1] * sd(B_bootShocks[,1])
b[,2] = b[,2] * sd(B_bootShocks[,2])
b[,3] = b[,3] * sd(B_bootShocks[,3])


resid_monte = residuals(varY) #Extract the residuals from the VAR model

#Set initial parameters of B for minimization
B0 = matrix(c(0.32, -0.4, 0.91,
              -0.32, 1.1, 2.6,
              0.28, 0.17, 0.76), nrow = 3, byrow = TRUE)

random_lambda = c(5,5,6) #Choose initial values of lambda for minimization
random_sigma = c(1,1,1)  #Choose initial values of sigma for minimization

ML_monte = mle_theta(resid_monte,l, B0,random_lambda, random_sigma) #Perform Maximum likelihood
maximize_absolute_product_columns(ML_monte$optimal_beta) #Permute the columms 
hessian = ML_monte$hessian #Extract the Hessian 
inv_hessian = - solve(hessian) #Invert the Hessian
SE = sqrt(abs(diag(inv_hessian))) #Obtain the SE from the Hessian






