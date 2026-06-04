# ------------------------------------------------------------------------------
# CONCENTRATION ANALYSIS, HHI and BC
# ------------------------------------------------------------------------------


# check for outliers in offering amount, (there are two)
top10_2014Q4 <- dt_final_long_boa %>%
  filter(qt_off == "2014 Q4") %>%     # keep only 2014 Q4
  arrange(desc(off_amt)) %>%          # sort by off_amt (descending)
  slice_head(n = 10)                  # keep the top 10


# ===============================================
# HHI, assuming that each UW parteciates evenly
# ===============================================


# ===============================
# 1. Prepare data
# ===============================
dt_concentration <- dt_final_long_boa %>%
  mutate(qt = as.yearqtr(qt, format = "%Y Q%q")) %>%
  arrange(qt)

# All quarters in the sample
quarters <- dt_concentration %>%
  distinct(qt) %>%
  arrange(qt) %>%
  pull(qt)

# ===============================
# 2. Helper: compute market-level distribution
# ===============================
compute_market_dist <- function(qt, data) {
  window <- data %>%
    filter(qt == !!qt, off_amt >= 0)
  
  if (nrow(window) == 0) return(NULL)
  
  # Total market issuance = sum of unique bonds
  total_market <- window %>%
    distinct(complete_cusip, .keep_all = TRUE) %>%
    summarize(total_off = sum(off_amt, na.rm = TRUE)) %>%
    pull(total_off)
  
  # Adjust each bond's amount evenly among underwriters
  dist <- window %>%
    group_by(complete_cusip) %>%
    mutate(n_underwriters = n(),
           adj_off_amt = off_amt / n_underwriters) %>%
    ungroup() %>%
    distinct(legal_name_canonical, complete_cusip, .keep_all = TRUE) %>%
    group_by(legal_name_canonical) %>%
    summarize(uw_amt = sum(adj_off_amt, na.rm = TRUE), .groups = "drop") %>%
    mutate(share = uw_amt / total_market) %>%
    ungroup()
  
  return(dist)
}

# ===============================
# 3. Compute market-level HHI per quarter (scaled by 10,000)
# ===============================
hhi_market <- map_dfr(quarters, function(qt) {
  dist <- compute_market_dist(qt, dt_concentration)
  
  if (is.null(dist)) return(tibble(qt = qt, HHI_market = NA_real_))
  
  tibble(
    qt = qt,
    HHI_market = 10000 * sum(dist$share^2, na.rm = TRUE)
  )
})

# ===============================
# 4. Plot market-level HHI
# ===============================
ggplot(hhi_market, aes(x = as.Date(qt), y = HHI_market)) +
  geom_line(color = "lightcoral", size = 1) +
  geom_point(color = "lightblue", size = 1.5) +
  labs(
    title = "Market-Level Herfindahl–Hirschman Index (HHI × 10,000)",
    x = "Quarter",
    y = "HHI (scaled by 10,000)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )



# ===============================================================================================================================
# HHI, scaling each underwriter share by the proportion of times it partecipates deals of that quarter and then rescaling
# ===============================================================================================================================


# ===============================
# 1. Prepare data
# ===============================
dt_concentration <- dt_final_long_boa %>%
  mutate(qt = as.yearqtr(qt, format = "%Y Q%q")) %>%
  arrange(qt)

# All quarters in the sample
quarters <- dt_concentration %>%
  distinct(qt) %>%
  arrange(qt) %>%
  pull(qt)

# ===============================
# 2. Helper: compute market-level distribution (weighted by bond participation share)
# ===============================
compute_market_dist <- function(qt, data) {
  window <- data %>%
    filter(qt == !!qt, off_amt >= 0)
  
  if (nrow(window) == 0) return(NULL)
  
  # Total market issuance = sum of unique bonds' amounts
  total_market <- window %>%
    distinct(complete_cusip, .keep_all = TRUE) %>%
    summarise(total_off = sum(off_amt, na.rm = TRUE)) %>%
    pull(total_off)
  
  total_bonds <- window %>%
    distinct(complete_cusip) %>%
    nrow()
  
  # Adjust each bond's amount evenly among underwriters
  dist <- window %>%
    group_by(complete_cusip) %>%
    mutate(n_underwriters = n(),
           adj_off_amt = off_amt / n_underwriters) %>%
    ungroup() %>%
    distinct(legal_name_canonical, complete_cusip, .keep_all = TRUE) %>%
    group_by(legal_name_canonical) %>%
    summarise(
      uw_amt = sum(adj_off_amt, na.rm = TRUE),
      n_bonds = n_distinct(complete_cusip),
      .groups = "drop"
    ) %>%
    mutate(
      share = uw_amt / total_market,
      bond_particip_share = n_bonds / total_bonds,
      scaled = share * bond_particip_share,
      final_share = scaled / sum(scaled, na.rm = TRUE)
    ) %>%
    ungroup()
  
  return(dist)
}

# ===============================
# 3. Compute market-level HHI per quarter (scaled by 10,000)
# ===============================
hhi_market_weighted <- map_dfr(quarters, function(qt) {
  dist <- compute_market_dist(qt, dt_concentration)
  
  if (is.null(dist)) return(tibble(qt = qt, HHI_market_weighted = NA_real_))
  
  tibble(
    qt = qt,
    HHI_market_weighted = 10000 * sum(dist$final_share^2, na.rm = TRUE)
  )
})

# ===============================
# 4. Plot market-level HHI (frequency-weighted)
# ===============================
ggplot(hhi_market_weighted, aes(x = as.Date(qt), y = HHI_market_weighted)) +
  geom_line(color = "purple", size = 1) +
  geom_point(color = "orange", size = 1.5) +
  labs(
    title = "Market-Level Herfindahl–Hirschman Index (HHI × 10,000)\nWeighted by Bond Participation Frequency",
    x = "Quarter",
    y = "Weighted HHI (scaled by 10,000)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

