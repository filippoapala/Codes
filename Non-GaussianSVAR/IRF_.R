# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(scales)

Populate_4D = function(boot_IRF, M, H, n_bootstrap){

irf_array = array(NA, dim = c(M, M, H, n_bootstrap))
# Loop through each bootstrap iteration and populate the 4D array

for (b in 1: n_bootstrap) {
  #Assuming each element in svar_bootstrap is a matrix of dimensions (M, H)
  #irf_array[, , , b] = ICAsvar_boot_IRF[[b]] #here you can place either the ICAsvar_boot_IRF object or the svar_boot_IRF obeject from the bootstrap 
  irf_array[, , , b] = boot_IRF[[b]]
}

  return(irf_array)
}


# ------------------------- COMPUTE THE IRF'S  ---------------------------

compute_irf_H = function (A_full, B, H) {
  
  
  M = nrow(A_full)  # Number of variables
  p = (ncol(A_full) - 1) / M  # Number of lags
  
  # Separate the constant term (not needed for IRFs)
  constants = A_full[, ncol(A_full)]  # The last column contains constants
  A_matrix = A_full[, 1:(ncol(A_full) - 1)]  # The rest contains VAR coefficients
  
  # Reshape the large matrix A_full into a list A for each lag
  A = vector("list", p)
  for (lag in 1:p) {
    A[[lag]] = A_matrix[, ((lag - 1) * M + 1):(lag * M)]
  }
  
  # Construct companion matrix C (size Mp x Mp)
  C = matrix(0, nrow = M * p, ncol = M * p)
  for (i in 1:p) {
    C[1:M, ((i - 1) * M + 1):(i * M)] = A[[i]]  # Place A matrices in the first row block
  }
  C[(M + 1):(M * p), 1:(M * (p - 1))] = diag(M * (p - 1))  # Set identity matrices for lag structure
  
  # Define selection matrix R (size M x Mp)
  R = matrix(0, nrow = M, ncol = M * p)
  diag(R) = 1
  
  # Compute SVAR IRFs up to horizon H
  IRFs = array(0, dim = c(M, M, H))  # Store IRFs in 3D array (M x M x H)
  
  # Compute the IRF at each horizon
  for (h in 0:(H - 1)) {
    Ch = matrix_power(C,h)  # Compute C to the power of h
    IRFs[, , h + 1] = (R %*% Ch %*% t(R)) %*% B  # Compute the IRF at horizon h
  }
  
  return(IRFs)
}



#Helper function to compute matrix power (R does not have this by default)
matrix_power = function(A, n) {
  if (n == 0) return(diag(nrow(A)))  # Identity matrix
  result = A
  for (i in 1:n) {
    result = result %*% A
  }
  
  return(result)
}




# -------------------------PLOT THE IRF'S ---------------------------

plot_IRF = function(sample_IRF, boot_IRF, M, H, n_bootstrap){
  
irf_array = Populate_4D(boot_IRF, M, H, n_bootstrap) #Populate the irf_array object
  

# Confidence levels for lower and upper bounds
conf_level = 0.68
conf_level2 = 0.9
conf_level3 = 0.95

# Initialize arrays to store results
irf_median = array(NA, dim = c(M, M, H))
irf_st = array(NA, dim = c(M, M, H))

lower_bound68 = array(NA, dim = c(M, M, H))
upper_bound68 = array(NA, dim = c(M, M, H))
lower_bound90 = array(NA, dim = c(M, M, H))
upper_bound90 = array(NA, dim = c(M, M, H))
lower_bound95 = array(NA, dim = c(M, M, H))
upper_bound95 = array(NA, dim = c(M, M, H))

# Prepare the data for plotting (median and confidence bands)
irf_data = data.frame()

for (i in 1:M) {
  
  for (j in 1:M) {
    
    for (h in 1:H) {
      # Extract all bootstrap IRFs for the (i,j) pair at horizon h
      boot_irfs = irf_array[i, j, h, ]  # 1D array of size N_bootstrap
      
      # Compute median and confidence bands 
      irf_median[i, j, h] = median(boot_irfs)
      irf_st[i,j,h] = sd(boot_irfs)
      lower_bound68[i, j, h] = quantile(boot_irfs, probs = (1 - conf_level) / 2) 
      upper_bound68[i, j, h] = quantile(boot_irfs, probs = (1 + conf_level) / 2) 
      lower_bound90[i, j, h] = quantile(boot_irfs, probs = (1 - conf_level2) / 2) 
      upper_bound90[i, j, h] = quantile(boot_irfs, probs = (1 + conf_level2) / 2) 
      lower_bound95[i, j, h] = quantile(boot_irfs, probs = (1 - conf_level3) / 2) 
      upper_bound95[i, j, h] = quantile(boot_irfs, probs = (1 + conf_level3) / 2) 
    }
    
  }
  
}




# Prepare the results for plotting in ggplot
for (i in 1:M) {
  
  for (j in 1:M) {
    
    for (h in 1:H) {
      
      irf_data = rbind(irf_data, data.frame(
        Response_Var = paste("Response Var", i),
        Impulse_Var = paste("Impulse Var", j),
        Horizon = h -1,  # Horizon starts from 0
        Median = irf_median[i, j, h],
        Lower = lower_bound68[i, j, h],
        Upper = upper_bound68[i, j, h],
        Lower2 = lower_bound90[i, j, h],
        Upper2 = upper_bound90[i, j, h],
        Lower3 = lower_bound95[i, j, h],
        Upper3 = upper_bound95[i, j, h],
        #sample = Dc_IRF[i, j, h] #Dc_IRF is obtained from compute_irf_H(A_full, B, H) where B is the estimated ICA structural matrix
        sample = sample_IRF[i,j,h] #Chol_IRF is obtained from compute_irf_H(A_full, B, H) where B is the estimated cholesky structural matrix
      ))
      
    }
  }
}



plots_list = list() #Create a list to store all plots


var_names = c("INF", "OG", "IR") #change with your own variable's names

# Loop through each combination of response and impulse variable
for (i in 1:M) {
  
  for (j in 1:M) {
    
    # Filter data for this specific response and impulse variable
    plot_subset = irf_data[irf_data$Response_Var == paste("Response Var", i) & 
                              irf_data$Impulse_Var == paste("Impulse Var", j), ]
    
    # Create a plot for each response variable to the impulse variable
      p = ggplot(plot_subset, aes(x = Horizon)) +
      geom_ribbon(aes(ymin = Lower3, ymax = Upper3), fill = "orange", alpha = 0.1) +  # 95% CI
      geom_ribbon(aes(ymin = Lower2, ymax = Upper2), fill = "orange", alpha = 0.23) +  # 90% CI
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "orange", alpha = 0.4) +  # 68% CI
      geom_line(aes(y = Median), size = 0.5, color = "orange") +  # Plot the median IRF line
      geom_line(aes(y = sample), size = 0.25 ,color = "black") +  # Sample IRF line
      geom_hline(yintercept = 0, linetype = "dashed",size = 0.1 ,color = "black") +  # Dashed line at zero
      labs(title = paste("Response of", var_names[i], "from shock in", var_names[j]),
           x = "Horizon",
           y = "Response") +
      scale_x_continuous(breaks = seq(0, max(plot_subset$Horizon), by = 1), 
                         labels = scales::comma_format()) +  # Finer x-axis units with clear labels
      scale_y_continuous(breaks = pretty_breaks(n = 6), 
                         labels = scales::comma_format()) +  # Finer y-axis units with clear labels
      
      theme_minimal() +
        
      theme(
        panel.border = element_rect(colour = "black", fill = NA, size = 0.1),  # Border around each plot
        panel.grid.major = element_line(color = "grey80", size = 0.01),  # Major grid lines
        plot.title = element_text(hjust = 0.5, size = 7, face = "bold"),  #Centered, smaller title
        axis.title = element_text(size = 7),   # Smaller axis titles
        axis.text = element_text(size = 6),    # Smaller axis labels
        axis.ticks = element_line(color = "black", size = 0.3),  # Black axis ticks
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")  # Adjusted margin around each plot
      )
    
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(col = "black", fill = NA, lwd = 1.5)),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
    
    
  
    plots_list[[length(plots_list) + 1]] = p #Add the plot to the list
  }
}


grid.arrange(grobs = plots_list, ncol = M, heights = unit(rep(2, M), "null")) #Arrange the plots in a grid (assuming 3x3 grid for M = 3)

}


