# ===============================
# 1. Compute number of underwriters per quarter
# ===============================
uw_per_quarter <- dt_final_long_boa %>%
  group_by(qt) %>%
  summarise(
    n_underwriters = n_distinct(legal_name_canonical),
    .groups = "drop"
  )


# ============================================================
#  Weight the number of underwriters for the number of deals
# ============================================================

uw_weighted <- dt_final_long_boa %>%
  #  Make sure each row is a unique deal-underwriter pair
  distinct(qt, legal_name_canonical, complete_cusip) %>%
  
  #  Count how many deals each underwriter participated in per quarter
  group_by(qt, legal_name_canonical) %>%
  summarise(deals = n_distinct(complete_cusip), .groups = "drop_last") %>%
  
  #  Compute total deals per quarter and relative weights
  mutate(
    total_deals = sum(deals, na.rm = TRUE),
    weight = deals / total_deals
  ) %>%
  
  #  Compute the weighted effective number of underwriters (1 / sum(w^2))
  summarise(
    weighted_n_underwriters = 1 / sum(weight^2, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# Make the plot
# ============================================================

ggplot(uw_weighted, aes(x = as.Date(qt), y = weighted_n_underwriters)) +
  geom_line(color = "royalblue", size = 1.2) +
  geom_point(color = "darkblue", size = 2) +
  labs(
    title = "Effective Number of Underwriters per Quarter (Weighted by Deal Activity)",
    x = "Quarter",
    y = "Weighted Number of Underwriters"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey70", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey85", linewidth = 0.3),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



# ===============================
# 2. Compute total market value per quarter (each bond counted once)
# ===============================
market_value_per_quarter <- dt_final_long_boa %>%
  group_by(qt, complete_cusip) %>%
  summarise(bond_amt = max(off_amt, na.rm = TRUE), .groups = "drop") %>%
  group_by(qt) %>%
  summarise(market_value = sum(bond_amt, na.rm = TRUE), .groups = "drop")

# ===============================
# 3. Plot number of underwriters
# ===============================
plot_underwriters <- ggplot(uw_per_quarter, aes(x = as.Date(qt), y = n_underwriters)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 1.5) +
  labs(
    title = "Number of Underwriters per Quarter",
    x = "Quarter",
    y = "Number of Underwriters"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
print(plot_underwriters)
# ===============================
# 4. Plot market value
# ===============================
plot_market_value <- ggplot(market_value_per_quarter, aes(x = as.Date(qt), y = market_value / 1e9)) +
  geom_line(color = "darkorange", size = 1) +
  geom_point(color = "orange", size = 1.5) +
  labs(
    title = "Total Market Value of Bonds per Quarter",
    x = "Quarter",
    y = "Market Value (Billions)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
print(plot_market_value)
# ===============================
# 1. Compute average bond size per quarter
# ===============================
avg_bond_size_per_quarter <- dt_final_long_boa %>%
  group_by(qt, complete_cusip) %>%
  summarise(bond_amt = max(off_amt, na.rm = TRUE), .groups = "drop") %>%
  group_by(qt) %>%
  summarise(
    avg_bond_size = mean(bond_amt, na.rm = TRUE),
    .groups = "drop"
  )

# ===============================
# 2. Plot average bond size
# ===============================
plot_avg_bond_size <- ggplot(avg_bond_size_per_quarter, aes(x = as.Date(qt), y = avg_bond_size / 1e6)) +
  geom_line(color = "purple", size = 1) +
  geom_point(color = "violet", size = 1.5) +
  labs(
    title = "Average Bond Size per Quarter",
    x = "Quarter",
    y = "Average Size (Millions)"
  ) +
  theme_minimal() +
  theme(
    
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

print(plot_avg_bond_size)


# ===============================
# 1. Prepare data
# ===============================
dt_nuws <- dt_final_long_boa %>%
  mutate(qt = as.yearqtr(qt, format = "%Y Q%q")) %>%
  arrange(qt)

# ===============================
# 2. Compute number of underwriters per deal
# ===============================
uws_per_deal <- dt_nuws %>%
  distinct(qt, complete_cusip, legal_name_canonical) %>%  # one row per underwriter-deal
  group_by(qt, complete_cusip) %>%
  summarise(n_uw = n_distinct(legal_name_canonical), .groups = "drop") %>%
  group_by(qt) %>%
  summarise(avg_uw_per_deal = mean(n_uw, na.rm = TRUE),
            n_deals = n(),
            .groups = "drop") %>%
  arrange(qt)

# ===============================
# 3. Plot average number of underwriters per deal
# ===============================
ggplot(uws_per_deal, aes(x = as.Date(qt), y = avg_uw_per_deal)) +
  geom_line(color = "purple", size = 1) +
  geom_point(color = "purple", size = 1.5) +
  labs(
    title = "Average Number of Underwriters per Deal per Quarter",
    x = "Quarter",
    y = "Average Number of Underwriters"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )


#compute the correlation between the average bond size and the average number of underwriters per deal
cor(avg_bond_size_per_quarter$avg_bond_size, uws_per_deal$avg_uw_per_deal)



# ===============================================
#  SAME NUMBER OF UNDERWRITERS BUT WITH BOND TYPE
# ===============================================

uw_per_quarter_hy_ig <- dt_final_long_boa %>%
  group_by(qt, legal_name_canonical) %>%
  summarise(
    high_yield_only = all(is_high_yield == 1),
    inv_grade_only = all(is_high_yield == 0),
    .groups = "drop"
  ) %>%
  # Label the underwriter type
  mutate(type = case_when(
    high_yield_only ~ "High-Yield Only",
    inv_grade_only  ~ "Investment-Grade Only",
    TRUE            ~ "Mixed"
  )) %>%
  filter(type != "Mixed") %>%  # Remove underwriters issuing both
  group_by(qt, type) %>%
  summarise(n_underwriters = n_distinct(legal_name_canonical), .groups = "drop")

# 2. Plot
ggplot(uw_per_quarter_hy_ig, aes(x = as.Date(qt), y = n_underwriters, color = type)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  labs(
    title = "Number of Underwriters per Quarter by Bond Type",
    x = "Quarter",
    y = "Number of Underwriters",
    color = "Underwriter Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )



# ==================================================================
#  AVERAGE NUMBER OF BOND'S ISSUED BY ONE UNDERWRITER BY QUARTER
# ==================================================================


avg_bonds_per_uw <- dt_final_long_boa %>%
  group_by(qt) %>%
  summarise(
    n_bonds = n_distinct(complete_cusip),               # unique bonds in the quarter
    n_underwriters = n_distinct(legal_name_canonical), # unique underwriters in the quarter
    avg_bonds_per_uw = n_bonds / n_underwriters,       # average bonds per underwriter
    .groups = "drop"
  )

ggplot(avg_bonds_per_uw, aes(x = as.Date(qt), y = avg_bonds_per_uw)) +
  geom_line(color = "purple", size = 1) +
  geom_point(color = "violet", size = 2) +
  labs(
    title = "Average Number of Bonds per Underwriter per Quarter",
    x = "Quarter",
    y = "Average Bonds per Underwriter"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# =====================================================================
#  PLOT THE RATIO OF FIRMS THAT ISSUE ONLY HIGH YIELD BONDS PER QUARTER
# =====================================================================


IG_HY_firms <- dt_final_long_boa_rs %>%
  distinct(qt_off, gvkey, is_high_yield) %>%   # unique firm-quarter-rating combos
  group_by(qt_off, gvkey) %>%                 # groups all the bonds issued by that firm in the same quarter
  summarise(
    only_HY = all(is_high_yield == 1),         # TRUE if each firm-quarter-rating combination has only high yield 
    .groups = "drop_last"
  ) %>%
  summarise(
    n_firms_total = n(),
    n_firms_HY_only = sum(only_HY, na.rm = TRUE),
    share_HY_only = n_firms_HY_only / n_firms_total,
    .groups = "drop"
  )


ggplot(IG_HY_firms, aes(x = qt_off, y = share_HY_only, group = 1)) +
  geom_line(color = "#007FFF") +
  geom_point(color = "#007FFF") +
  labs(
    title = "Share of firms issuing only HY bonds per quarter",
    x = "Quarter",
    y = "Share of HY-only issuers"
  ) +
  theme_minimal()

# =====================================================================
#  TRY SMOOTHING VIA AVERAGING CLOSE POINTS
# =====================================================================

IG_HY_firms$share_HY_only_smooth <- rollmean(IG_HY_firms$share_HY_only, k = 3, fill = NA, align = "center")

ggplot(IG_HY_firms, aes(x = qt_off)) +
  geom_line(aes(y = share_HY_only_smooth), color = "#007FFF") +
  geom_point(aes(y = share_HY_only), color = "#007FFF") +
  labs(title = "Share of firms issuing only HY bonds per quarter (smoothed)",
       x = "Quarter",
       y = "Share of HY-only issuers") +
  theme_minimal()
