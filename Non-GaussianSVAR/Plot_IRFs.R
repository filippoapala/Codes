


#YOU NEED THIS TO PLOT THE IRF'S
irf_array <- array(NA, dim = c(M, M, H, n_bootstrap))

# Loop through each bootstrap iteration and populate the 4D array
for (b in 1: n_bootstrap) {
  #Assuming each element in svar_bootstrap is a matrix of dimensions (M, H)
  irf_array[, , , b] <- ICAsvar_boot_IRF[[b]] #here you put the object you obtained from the bootstrap code
}



#IT COMPUTES THE IRF'S
compute_irf_H <- function (A_full, B, H) {
  # A_full: Large matrix of VAR(p) coefficients including constant (Mx(M*p + 1))
  # B: Structural impact matrix (MxM)
  # H: the number of horizons for IRF
  # Returns IRFs for H horizons
  
  # Get dimensions
  M <- nrow(A_full)  # Number of variables
  p <- (ncol(A_full) - 1) / M  # Number of lags
  
  # Separate the constant term (not needed for IRFs)
  constants <- A_full[, ncol(A_full)]  # The last column contains constants
  A_matrix <- A_full[, 1:(ncol(A_full) - 1)]  # The rest contains VAR coefficients
  
  # Reshape the large matrix A_full into a list A for each lag
  A <- vector("list", p)
  for (lag in 1:p) {
    A[[lag]] <- A_matrix[, ((lag - 1) * M + 1):(lag * M)]
  }
  
  # Construct companion matrix C (size Mp x Mp)
  C <- matrix(0, nrow = M * p, ncol = M * p)
  for (i in 1:p) {
    C[1:M, ((i - 1) * M + 1):(i * M)] <- A[[i]]  # Place A matrices in the first row block
  }
  C[(M + 1):(M * p), 1:(M * (p - 1))] <- diag(M * (p - 1))  # Set identity matrices for lag structure
  
  # Define selection matrix R (size M x Mp)
  R <- matrix(0, nrow = M, ncol = M * p)
  diag(R) <- 1
  
  # Compute SVAR IRFs up to horizon H
  IRFs <- array(0, dim = c(M, M, H))  # Store IRFs in 3D array (M x M x H)
  
  # Compute the IRF at each horizon
  for (h in 0:(H - 1)) {
    Ch <- matrix_power(C, h)  # Compute C to the power of h
    IRFs[, , h + 1] <- (R %*% Ch %*% t(R)) %*% B  # Compute the IRF at horizon h
  }
  
  return(IRFs)
}


 #Helper function to compute matrix power (R does not have this by default)
  matrix_power <- function(A, n) {
    if (n == 0) return(diag(nrow(A)))  # Identity matrix
      result <- A
   for (i in 2:n) {
      result <- result %*% A
   }
      
  return(result)
}



# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Assuming svar_boot_IRF is a 3D array: (M, H, N_bootstrap)

# Number of variables and horizons
M <- 3 #Replace with number of variables
H <- 15  #Replace with Number of horizons


conf_level = 0.9

# Initialize arrays to store results
irf_median <- array(NA, dim = c(M, M, H))
irf_st = array(NA, dim = c(M,M, H))
lower_bound <- array(NA, dim = c(M, M, H))
upper_bound <- array(NA, dim = c(M, M, H))
# Prepare the data for plotting (median and confidence bands)
irf_data <- data.frame()

for (i in 1:M) {
  for (j in 1:M) {
    for (h in 1:H) {
      # Extract all bootstrap IRFs for the (i,j) pair at horizon h
      boot_irfs <- irf_array[i, j, h, ]  # 1D array of size N_bootstrap
      
      # Compute median and confidence bands (5th and 95th percentiles)
      irf_median[i, j, h] <- median(boot_irfs)
      irf_st[i,j,h] = sd(boot_irfs)
      lower_bound[i, j, h] <- quantile(boot_irfs, probs = (1 - conf_level) / 2) # e.g., 5th percentile
      upper_bound[i, j, h] <- quantile(boot_irfs, probs = (1 + conf_level) / 2) # e.g., 95th percentile
    }
  }
}



# Prepare the results for plotting in ggplot
for (i in 1:M) {
  for (j in 1:M) {
    for (h in 1:H) {
      irf_data <- rbind(irf_data, data.frame(
        Response_Var = paste("Response Var", i),
        Impulse_Var = paste("Impulse Var", j),
        Horizon = h -1,  # Horizon starts from 0
        Median = irf_median[i, j, h],
        Lower = lower_bound[i, j, h],
        Upper = upper_bound[i, j, h],  
        sample = CV_IRF[i,j,h] #Replace with your sample estimates
      ))
    }
  }
}



# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Assuming 'irf_data' is already prepared as per the previous example
# Now we will create the plots for each response and impulse variable

# Create a list to store all plots
plots_list <- list()


var_names <- c("INF", "OG", "IR") #Replace with your variable's name

# Loop through each combination of response and impulse variable
for (i in 1:M) {
  for (j in 1:M) {
    # Filter data for this specific response and impulse variable
    plot_subset <- irf_data[irf_data$Response_Var == paste("Response Var", i) & 
                              irf_data$Impulse_Var == paste("Impulse Var", j), ]
    
    # Create a plot for each response variable to the impulse variable
    p <- ggplot(plot_subset, aes(x = Horizon)) +
      geom_line(aes(y = Lower), linetype = "dashed", color = "black") +  # Dashed line for lower CI
      geom_line(aes(y = Upper), linetype = "dashed", color = "black") +  # Dashed line for upper CI
      geom_line(aes(y = Median), size = 0.5, color = "red") +  # Plot the median IRF line
      geom_line(aes(y = sample), size = 0.2,linetype = "dashed" ,color = "blue") +  # Sample IRF line
      geom_hline(yintercept = 0, linetype = "dashed",size = 0.5 ,color = "grey") +  # Dashed line at zero
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
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),  # Smaller title size
        
        axis.title = element_text(size = 7),   # Smaller axis title
        axis.text = element_text(size = 6),    # Smaller axis labels
        
        axis.ticks = element_line(color = "black", size = 0.3),  # Black axis ticks
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")  # Adjusted margin around each plot
      ) 
    
    
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(col = "black", fill = NA, lwd = 1.5)),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
    
    
    
    # Add the plot to the list
    plots_list[[length(plots_list) + 1]] <- p
  }
}

# Arrange the plots in a grid (assuming 3x3 grid for M = 3)
grid.arrange(grobs = plots_list, ncol = M, heights = unit(rep(2, M), "null"))



