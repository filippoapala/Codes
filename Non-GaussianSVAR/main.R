install.packages("ggplot2")
install.packages("gridExtra")
install.packages("scales")
library(svars)  #or steady ICA
library(vars)
library(StReg)



# -------------------------LOAD THE DATA   ---------------------------

data_neriCompleto = readxl::read_xlsx("Daticompleti.xlsx")
data_neriCompleto = data_neriCompleto[, c("Quarterly inflation", "Output gap CL", "Euribor 3 months")] #Order the variables

# -------------------------LOAD THE FUNCTIONS   ---------------------------
source("SVAR_Bootstrap_.R") 
source("MaximumLikelihood_.R")
source("IRF_.R")
source("FunctionsTools_.R")


# -------------------------ESTIMATE THE VAR MODEL  ---------------------------

M = 3
H = 15


var_neri = VARselect(data_neriCompleto, lag.max = 10, type = "const") #lag selection
var_model2T = vars::VAR(data_neriCompleto, p = 2, type = "const") #estimate the VAR with 2 lags
A_full = var_coefficients(var_model2T, data_neriCompleto) #extract the VAR coefficients
residuals_var2T = residuals(var_model2T) #extract the residuals
STvar_model2T = StVAR(as.matrix(data_neriCompleto),Trend=1,lag=2,v=8,maxiter=1000,meth="BFGS",hes="TRUE",init="na") # Student-t VAR model estimated with ML, don't work well with simulated data


# -------------------------CHOLEKSY SVAR WITH BOOTSTRAP ---------------------------

B_chol = t(chol(cov(residuals_var2T))) #cholesky factor on the sample
Chol_IRF = compute_irf_H(A_full, B_chol, 15) 
svar_boot_IRF = generateBootIRFs_Chol(3000, 15, 3, data_neriCompleto, var_model2T)
plot_IRF(Chol_IRF, svar_boot_IRF, 3, 15, 3000)

# -------------------------ICA SVAR WITH BOOTSTRAP ---------------------------

B_hat = id.dc(var_model2T, PIT = FALSE)$B #ICA by distance covariance on the sample
DC_IRF = compute_irf_H(A_full, B_hat, 15) 
svar_boot_IRF = generateBootIRFs_ICA(3000, 15, 3, data_neriCompleto, var_model2T, 1)
plot_IRF(DC_IRF, svar_boot_IRF$ICAsvar_boot_IRF, 3, 15, 3000)

# -------------------------ICA SVAR STUDENT-T MAXIMUM LIKELIHOOD ---------------------------

random_B = B_chol       #set initial conditions for the parameters
random_lambda = c(6,4,3)
random_sigma = c(1,1,1)
l = Unit_vector_gen(3)

ML = mle_theta(residuals_var2T, l, random_B, random_lambda, random_sigma) #ICA by Student-t maximum likelihood
B_ML = ML$optimal_beta
ML_IRF = compute_irf_H(A_full, B_ML, 15) 
svar_boot_IRF = generateBootIRFs_ICA(10, 15, 3, data_neriCompleto, var_model2T, 0)
plot_IRF(ML_IRF, svar_boot_IRF$ICAsvar_boot_IRF, 3, 15, 10)
