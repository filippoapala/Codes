library(YieldCurve)

# =====================================
# IMPORT YIELD CURVE DATA
# =====================================


#import yield curve data
files <- list.files(path = "/Users/filippopalandri/Desktop/UNI/Pre-Doc/UW_Project/DATA/Risk_free", pattern = "\\.csv$", full.names = TRUE)
dfs <- lapply(files, read.csv)

interest_rates <- bind_rows(dfs)
interest_rates$Date <- gsub("/", "-", interest_rates$Date)



# =====================================
# OBTAIN THE YIELD CURVE FOR EACH BOND
# =====================================


#create new dataframe with unique cusip rows
dt_final_long_boa_rs = dt_final_long_boa %>%
  distinct(complete_cusip, .keep_all = TRUE)   # keep only one row per bond

#create a new dataframe with unique cusip dt
dt_one = dt %>%
  distinct(complete_cusip, .keep_all = TRUE)   # keep only one row per bond

#join columns from the dt_one dataframe
dt_final_long_boa_rs = dt_final_long_boa_rs %>% # obtain the maturity for each bond
  left_join(
    dt_one %>% select(complete_cusip, maturity, offering_date),
    by = "complete_cusip"
  )

#compute the maturity length as the difference between the maturity date and the offering date
dt_final_long_boa_rs <- dt_final_long_boa_rs %>%
  mutate(
    # Convert to Date if not already
    offering_date = as.Date(offering_date),
    maturity = as.Date(maturity),
    
    # Compute difference in years and months
    maturity_length_years = year(maturity) - year(offering_date),
    maturity_length_months = month(maturity) - month(offering_date),
    
    # Adjust if months negative
    maturity_length = maturity_length_years + maturity_length_months / 12
  )



#change the name of the column date and harmonize the date format between the bond dataframe and the interest rate one
interest_rates <- interest_rates %>%
  rename(
    offering_date = Date
  )
interest_rates <- interest_rates %>%
  mutate(offering_date = format(as.Date(offering_date, format = "%m-%d-%Y"), "%Y-%m-%d"))

interest_rates<- interest_rates %>%
  mutate(offering_date = as.Date(offering_date, format = "%Y-%m-%d"))


# Create new column order:
n <- ncol(interest_rates)
# 1:2 stay, n (last column) goes third, then columns 3:(n-1) follow
new_order <- c(1, 2, n, 3:(n-1))
# Reorder the data frame
interest_rates <- interest_rates[, new_order]


#join by offering_date
dt_final_long_boa_rs = left_join(dt_final_long_boa_rs, interest_rates, by = "offering_date")

#change the name of the columns
colnames(dt_final_long_boa_rs)[20] <- "0.083.year"
colnames(dt_final_long_boa_rs)[21] <- "0.16.year"
colnames(dt_final_long_boa_rs)[22] <- "0.25.year"
colnames(dt_final_long_boa_rs)[23] <- "0.5.year"
colnames(dt_final_long_boa_rs)[24] <- "1.year"
colnames(dt_final_long_boa_rs)[25] <- "2.year"
colnames(dt_final_long_boa_rs)[26] <- "3.year"
colnames(dt_final_long_boa_rs)[27] <- "5.year"
colnames(dt_final_long_boa_rs)[28] <- "7.year"
colnames(dt_final_long_boa_rs)[29] <- "10.year"
colnames(dt_final_long_boa_rs)[30] <- "20.year"
colnames(dt_final_long_boa_rs)[31] <- "30.year"


# ==================================================================================
# MATCH THE CLOSEST RISK FREE RATE TO EACH MATURITY, (should work but check better)
# ==================================================================================


# Extract numeric maturities from column names like "1.year", "5.year"
maturity_cols <- grep("\\.year$", names(dt_final_long_boa_rs), value = TRUE)
maturity_years <- as.numeric(sub("\\.year$", "", maturity_cols))

# For each row, pick the yield value with the closest maturity
dt_final_long_boa_rs$risk_free_rate <- apply(dt_final_long_boa_rs, 1, function(row) {
  target <- as.numeric(row["maturity_length"])
  yields <- as.numeric(row[maturity_cols])
  
  # Find index of the closest maturity
  closest_index <- which.min(abs(maturity_years - target))
  
  # Return yield for that maturity
  return(yields[closest_index])
})


#create the variable offering_spread as the differcence between the offering yield and the risk free rate
dt_final_long_boa_rs <- dt_final_long_boa_rs %>%
  mutate(offering_spread = offering_yield - risk_free_rate)




# =================================================================================================
# ??????? ESTIMATE THE NELSON-SIEGEL MODEL IN ORDER TO RECOVER THE EXACT TIME TO MATURITIES ???????
# =================================================================================================


obs_interest_rates = interest_rates[612, 2:ncol(interest_rates)]
x = c(1,2,3, 6, 12, 24, 36, 60, 81.6,84, 120, 240, 360) / 12
obs_interest_rates = obs_interest_rates[,ncol(obs_interest_rates)]

NSparameters = Nelson.Siegel(obs_interest_rates, time_years)




NS_yield <- function(params, tau) {
  # params: numeric vector of length 4: c(beta0, beta1, beta2, lambda)
  beta0 <- params[1]
  beta1 <- params[2]
  beta2 <- params[3]
  lambda <- params[4]
  
  # Nelson-Siegel formula
  y <- beta0 + 
    beta1 * (1 - exp(-lambda * tau)) / (lambda * tau) +
    beta2 * ((1 - exp(-lambda * tau)) / (lambda * tau) - exp(-lambda * tau))
  
  return(y)
}


risk_free = NS_yield(NSparameters, 6.8)
obs_interest_rates = cbind(obs_interest_rates, risk_free)


x = c(1,2,3, 6, 12, 24, 36, 60, 81.6,84, 120, 240, 360) / 12
y <- as.numeric(obs_interest_rates[1, ])  # first row as y (numeric)

# Basic plot
plot(x, y, type = "b", pch = 19, col = "blue",
     xlab = "Maturity", ylab = "Yield",
     main = "Yield Curve for First Row")
