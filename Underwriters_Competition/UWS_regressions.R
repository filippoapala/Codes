

library(broom)
library(tidyverse)
# ====================================================
#             REGRESSIONS WITH LAGGED BOND RATING
# ====================================================


# FIRST: CREATE A 2 QUARTERS AHEAD BOND RATING MEASURE
# SECOND PERFORM THE ANALYSIS (THE RISK FREE R SCRIPT NEEDS TO BE RUN BEFORE PERFORMING THE REGRESSIONS)


############## !!!!!! FIRST RUN THE risk_free R script ----> TO GET dt_final_long_boa_rs !!!!!! ########### 

# 1) FILTER DT ONLY FOR THE FINAL BONDS
# Keep only rows in dt where complete_cusip exists in dt_final_long_boa (the time series is still present)
dt_filtered <- dt %>%
  semi_join(dt_final_long_boa_rs %>% select(complete_cusip) %>% distinct(),
            by = "complete_cusip")


# 2) APPEND THE LAGGED BOND RATING COLUMN TO DT 
dt_filtered <- dt_filtered %>%
  mutate(qt = as.Date(as.yearqtr(qt)))

# 3) Create lagged bond rating dataset (2 quarters ahead)
dt_lagged <- dt_filtered %>%
  mutate(
    target_quarter = as.yearqtr(offering_date) + 0.5  # 2 quarters ahead
  ) %>%
  left_join(
    dt_filtered %>%
      mutate(qt_quarter = as.yearqtr(qt)) %>%
      select(complete_cusip, qt_quarter, bond.rat),
    by = c("complete_cusip", "target_quarter" = "qt_quarter")
  ) %>%
  rename(bond.rat_plus2q = bond.rat.y) %>%
  select(complete_cusip, qt, offering_date, bond.rat = bond.rat.x, bond.rat_plus2q) %>%
  mutate(
    # if bond.rat_plus2q is NA, fill with current bond.rat
    bond.rat_plus2q = coalesce(bond.rat_plus2q, bond.rat)
  )

# 4) Join back to the main dataset
dt_final_long_boa_rs <- dt_final_long_boa_rs %>%
  left_join(
    dt_lagged %>%
      select(complete_cusip, bond.rat_plus2q) %>%
      group_by(complete_cusip) %>%
      summarise(bond.rat_plus2q = first(bond.rat_plus2q), .groups = "drop"),
    by = "complete_cusip"
  ) %>%
  rename(bond.rat_lagged = bond.rat_plus2q)

#create new column that takes the value 1 high yields bonds
dt_final_long_boa_rs <- dt_final_long_boa_rs %>%
  mutate(is_high_yield_lagged = if_else(bond.rat_lagged > 10, 1, 0))

 


# ================================================================================
#  Run regression by quarter with fixed effect for HY/IG (within group variation)
# ================================================================================


#takes the rows that only contain both the offering_spread and the bond rating (also make the distinct but not necessary as dt_final_long_boa_rs is already cusip unique per row)
dt_clean <- dt_final_long_boa_rs %>%
  filter(!is.na(offering_spread), !is.na(bond.rat_lagged)) %>%
  distinct(complete_cusip, .keep_all = TRUE)

#remove the quarters that produce outliers in the estimates and distort the scale
dt_clean <- dt_clean %>%
  filter(!qt_off %in% c("2006 Q4", "2019 Q4"))


#performs the regression of the bond rating one the offering_spread with group specific fixed effects
reg_results <- dt_clean %>%
  group_by(qt_off) %>%
  group_modify(~tidy(lm(offering_spread ~ bond.rat_lagged + factor(is_high_yield_lagged), data = .x))) %>%
  filter(term == "bond.rat_lagged") %>%
  ungroup()

# ==================================================
# Plot the coefficients over time
# ==================================================
ggplot(reg_results, aes(x = qt_off, y = estimate, group = 1)) +
  geom_ribbon(aes(ymin = estimate - 1.96*std.error,
                  ymax = estimate + 1.96*std.error),
              fill = "lightcoral", alpha = 0.3, na.rm = TRUE) +
  geom_line(color = "coral", size = 1, na.rm = TRUE) +
  geom_point(color = "red", size = 2, na.rm = TRUE) +
  labs(
    title = "Bond rating on the offering spread with IG/HY fixed effects, coefficient estimates quarter by quarter",
    x = "Quarter",
    y = "Coefficient (slope)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# =================================================================
#  REGRESSION WITH ONLY THE DIFFERENCES IN MEAN BETWEEN HY AND IG
# =================================================================


reg_results_dummy <- dt_clean %>%
  group_by(qt_off) %>%
  group_modify(~ tidy(lm(offering_spread ~ factor(is_high_yield_lagged), data = .x))) %>%
  filter(term == "factor(is_high_yield_lagged)1") %>%  # Keep only the dummy coefficient
  ungroup() %>%
  mutate(
    ci_lower = estimate - 1.96 * std.error,
    ci_upper = estimate + 1.96 * std.error
  )



ggplot(reg_results_dummy, aes(x = qt_off, y = estimate, group = 1)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              fill = "skyblue", alpha = 0.3) +
  geom_line(color = "powderblue", size = 1) +
  geom_point(color = "deepskyblue1", size = 2) +
  labs(
    title = "High-Yield Dummy Coefficient on Offering Spread by Quarter",
    x = "Quarter",
    y = "Coefficient"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))










# ============================================================
#   PLOT THE AVERAGE BOND RATING PER QUARTER 
# ============================================================

# Calculate average bond rating per quarter
average_bond_rat_market <- dt_final_long_boa_rs %>%
  distinct(complete_cusip, .keep_all = TRUE) %>%  # keep unique bonds
  group_by(qt_off) %>%
  summarise(
    avg_bond_rat = mean(bond.rat, na.rm = TRUE),
    n_bonds = n(),
    .groups = "drop"
  ) %>%
  mutate(qt_fac = factor(qt_off, levels = unique(qt_off)))



ggplot(average_bond_rat_market, aes(x = qt_fac, y = avg_bond_rat, group = 1)) +
  geom_line(size = 1, color = "purple") +
  geom_point(size = 2, color = "purple") +
  labs(
    x = "Quarter",
    y = "Average bond rating",
    title = "Average bond rating per quarter — Market Level"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
  )

# ============================================================
#   PLOT THE AVERAGE BOND RATING PER QUARTER FOR HY AND IG
# ============================================================



# Step 1: Calcola la media per trimestre e categoria
average_bond_rat <- dt_final_long_boa_rs %>%
  distinct(complete_cusip, .keep_all = TRUE) %>%
  group_by(qt_off, is_high_yield_lagged) %>%
  summarise(
    avg_bond_rat = mean(bond.rat_lagged, na.rm = TRUE),
    n_bonds = n(),
    .groups = "drop"
  ) %>%
  mutate(
    category = if_else(is_high_yield_lagged == 1, "High Yield", "Investment Grade"),
    qt_fac = factor(qt_off, levels = unique(qt_off))
  )

# Step 2: Time series plot
ggplot(average_bond_rat, aes(x = qt_fac, y = avg_bond_rat, color = category, group = category)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "Quarter",
    y = "Average bond rating",
    color = NULL,
    title = "Average bond rating per quarter: High Yield vs Investment Grade"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "top"
  )
