# =====================================================================
#  MUTUAL FUNDS ANALYSIS AND SHARE TRACKING
# =====================================================================

holdings <- holdings %>%
  mutate(
    #identify the type of investor
    inv_group = case_when(
      fund_type %in% c("BAL", "MMM", "MUT", "END", "QUI", "FOF") ~ "Mutual Funds",
      fund_type %in% c("INS", "LIN", "PIN")                      ~ "Insurance",
      fund_type %in% c("CPF", "GPE", "UPE")                      ~ "Pension Funds",
      TRUE ~ "Other"
    ),
    # Portiamo la data a formato "quarter"
    qt = as.yearqtr(date)
  )

#filter all the rows that contains BlackRock or similar in parent_name
blackrock = holdings %>%
  filter(str_detect(parent_name, regex("BlackRock", ignore_case = TRUE)))

#filter all the rows that contains Fidelity or similar in parent_name
fidelity = holdings %>%
  filter(str_detect(parent_name, regex("Fidelity", ignore_case = TRUE)))

#filter all the rows that contains BlackRock or similar in parent_name
state_street = holdings %>%
  filter(str_detect(parent_name, regex("State Street", ignore_case = TRUE)))

#filter all the rows that contains BlackRock or similar in parent_name
vanguard = holdings %>%
  filter(str_detect(parent_name, regex("Vanguard", ignore_case = TRUE)))


# Check uniqueness
is_unique_per_qt <- holdings %>%
  count(qt, cusip) %>%
  summarise(all_unique = all(n == 1)) %>%
  pull(all_unique)

is_unique_per_qt


dt_final_long_boa_rs <- dt_final_long_boa_rs %>%
  left_join(
    dt_unique %>% select(complete_cusip, cusip_8),
    by = "complete_cusip"
  )



holdings_HY = holdings %>%
  left_join(dt_final_long_boa_rs %>% select(cusip_8, is_high_yield),
            by = c("cusip" = "cusip_8")
  )


# =====================================================================
#  COMPUTE AND PLOT THE SHARE OF THE MARKET HOLD BY BLACKROCK OVER TIME
# =====================================================================

blackrock_shares <- holdings %>%
  group_by(qt) %>%
  summarise(
    total_holding = sum(dol_holding, na.rm = TRUE),
    black_holding = sum(
      dol_holding[str_detect(parent_name, regex("BlackRock", ignore_case = TRUE))],
      na.rm = TRUE
    ),
    black_share = black_holding / total_holding * 100,
    .groups = "drop"
  )

blackrock_shares_HY <- holdings_HY %>%
  group_by(qt,is_high_yield) %>%
  summarise(
    total_holding = sum(dol_holding, na.rm = TRUE),
    black_holding = sum(
      dol_holding[str_detect(parent_name, regex("BlackRock", ignore_case = TRUE))],
      na.rm = TRUE
    ),
    black_share = black_holding / total_holding * 100,
    .groups = "drop"
  )


ggplot(blackrock_shares_HY, aes(x = as.Date(qt), y = black_share, color = factor(is_high_yield))) +
  geom_line(size = 1) +
  geom_point(size = 1.8) +
  scale_color_manual(
    name = "Market type",
    values = c("0" = "#1E90FF", "1" = "#FF4500"),
    labels = c("Investment Grade", "High Yield")
  ) +
  labs(
    title = "BlackRock's Market Share in IG vs HY Bond Holdings",
    x = "Quarter",
    y = "Market share (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.major = element_line(color = "grey60", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey70", linewidth = 0.2),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



#counting cusip only once per parent
blackrock_shares <- holdings %>%
  distinct(qt, cusip, parent_name, dol_holding) %>%   # ensure one row per bond per quarter per parent
  group_by(qt) %>%
  summarise(
    total_holding = sum(dol_holding, na.rm = TRUE),            # total across all unique bonds
    black_holding = sum(
      dol_holding[str_detect(parent_name, regex("BlackRock", ignore_case = TRUE))],
      na.rm = TRUE
    ),
    black_share = black_holding / total_holding * 100,
    .groups = "drop"
  )


ggplot(blackrock_shares, aes(x = qt, y = black_share, group = 1)) +
  geom_line(color = "#007FFF", size = 1.2) +
  geom_point(color = "#007FFF", size = 2) +
  labs(
    title = "BlackRock's Share of Total Holdings Over Time",
    x = "Quarter",
    y = "Share of Total Holdings (%)"
  ) +
  theme_minimal() +
  theme(

    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey85")
  )



# =====================================================================
#  COMPUTE AND PLOT THE SHARE OF THE MARKET HOLD BY FIDELITY OVER TIME
# =====================================================================

#some bonds in 2002 Q3 are counted more than once
fidelity_shares <- holdings %>%
  group_by(qt) %>%
  summarise(
    total_holding = sum(dol_holding, na.rm = TRUE),
    fidelity_holding = sum(
      dol_holding[str_detect(parent_name, regex("Fidelity", ignore_case = TRUE))],
      na.rm = TRUE
    ),
    fidelity_shares = fidelity_holding / total_holding * 100,
    .groups = "drop"
  )

fidelity_shares_HY <- holdings_HY %>%
  group_by(qt, is_high_yield) %>%
  summarise(
    total_holding = sum(dol_holding, na.rm = TRUE),
    fidelity_holding = sum(
      dol_holding[str_detect(parent_name, regex("Fidelity", ignore_case = TRUE))],
      na.rm = TRUE
    ),
    fidelity_shares = fidelity_holding / total_holding * 100,
    .groups = "drop"
  )


ggplot(fidelity_shares_HY, aes(x = as.Date(qt), y = fidelity_shares, color = factor(is_high_yield))) +
  geom_line(size = 1) +
  geom_point(size = 1.8) +
  scale_color_manual(
    name = "Market type",
    values = c("0" = "purple", "1" = "orange"),
    labels = c("Investment Grade", "High Yield")
  ) +
  labs(
    title = "Fidelity's Market Share in IG vs HY Bond Holdings",
    x = "Quarter",
    y = "Market share (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
    axis.text.x = element_text(angle = 45, hjust = 1)

  )

#counting cusip only once per parent
fidelity_shares <- holdings %>%
  distinct(qt, cusip, parent_name, dol_holding) %>%   # ensure one row per bond per quarter per parent
  group_by(qt) %>%
  summarise(
    total_holding = sum(dol_holding, na.rm = TRUE),            # total across all unique bonds
    fidelity_holding = sum(
      dol_holding[str_detect(parent_name, regex("Fidelity", ignore_case = TRUE))],
      na.rm = TRUE
    ),
    fidelity_shares = fidelity_holding / total_holding * 100,
    .groups = "drop"
  )


ggplot(fidelity_shares, aes(x = qt, y = fidelity_shares, group = 1)) +
  geom_line(color = "orange", size = 1.2) +
  geom_point(color = "orange", size = 2) +
  labs(
    title = "Fidelity's Share of Total Holdings Over Time",
    x = "Quarter",
    y = "Share of Total Holdings (%)"
  ) +
  theme_minimal() +
  theme(
    ,
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey85")
  )



# ========================================================================  
#  COMPUTE AND PLOT THE SHARE OF THE MARKET HOLD BY STATE STREET OVER TIME
# ========================================================================

state_street_shares <- holdings %>%
  group_by(qt) %>%
  summarise(
    total_holding = sum(dol_holding, na.rm = TRUE),
    state_street_holding = sum(
      dol_holding[str_detect(parent_name, regex("State Street", ignore_case = TRUE))],
      na.rm = TRUE
    ),
    state_street_shares = state_street_holding / total_holding * 100,
    .groups = "drop"
  )


state_street_shares_HY <- holdings_HY %>%
  group_by(qt,is_high_yield) %>%
  summarise(
    total_holding = sum(dol_holding, na.rm = TRUE),
    state_street_holding = sum(
      dol_holding[str_detect(parent_name, regex("State Street", ignore_case = TRUE))],
      na.rm = TRUE
    ),
    state_street_shares = state_street_holding / total_holding * 100,
    .groups = "drop"
  )


ggplot(state_street_shares_HY, aes(x = as.Date(qt), y = state_street_shares, color = factor(is_high_yield))) +
  geom_line(size = 1) +
  geom_point(size = 1.8) +
  scale_color_manual(
    name = "Market type",
    values = c("0" = "#007FFF", "1" = "pink"),
    labels = c("Investment Grade", "High Yield")
  ) +
  labs(
    title = "State Street's Market Share in IG vs HY Bond Holdings",
    x = "Quarter",
    y = "Market share (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



#counting cusip only once per parent
state_street_shares <- holdings %>%
  distinct(qt, cusip, parent_name, dol_holding) %>%   # ensure one row per bond per quarter per parent
  group_by(qt) %>%
  summarise(
    total_holding = sum(dol_holding, na.rm = TRUE),            # total across all unique bonds
    state_street_holding = sum(
      dol_holding[str_detect(parent_name, regex("State Street", ignore_case = TRUE))],
      na.rm = TRUE
    ),
    state_street_shares = state_street_holding / total_holding * 100,
    .groups = "drop"
  )



ggplot(state_street_shares, aes(x = qt, y = state_street_shares, group = 1)) +
  geom_line(color = "purple", size = 1.2) +
  geom_point(color = "purple", size = 2) +
  labs(
    title = "State Street's Share of Total Holdings Over Time",
    x = "Quarter",
    y = "Share of Total Holdings (%)"
  ) +
  theme_minimal() +
  theme(
    
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey85")
  )




# ========================================================================  
#  COMPUTE AND PLOT THE SHARE OF THE MARKET HOLD BY VANGUARD OVER TIME
# ========================================================================

Vanguard_shares <- holdings %>%
  group_by(qt) %>%
  summarise(
    total_holding = sum(dol_holding, na.rm = TRUE),
    Vanguard_holding = sum(
      dol_holding[str_detect(parent_name, regex("Vanguard", ignore_case = TRUE))],
      na.rm = TRUE
    ),
    Vanguard_shares = Vanguard_holding / total_holding * 100,
    .groups = "drop"
  )

Vanguard_shares = Vanguard_shares %>%
  filter(qt >= "2002 Q1") 


Vanguard_shares_HY <- holdings_HY %>%
  group_by(qt, is_high_yield) %>%
  summarise(
    total_holding = sum(dol_holding, na.rm = TRUE),
    Vanguard_holding = sum(
      dol_holding[str_detect(parent_name, regex("Vanguard", ignore_case = TRUE))],
      na.rm = TRUE
    ),
    Vanguard_shares = Vanguard_holding / total_holding * 100,
    .groups = "drop"
  )

Vanguard_shares_HY = Vanguard_shares_HY %>%
  filter(qt >= "2002 Q1") 


ggplot(Vanguard_shares_HY, aes(x = as.Date(qt), y = Vanguard_shares, color = factor(is_high_yield))) +
  geom_line(size = 1) +
  geom_point(size = 1.8) +
  scale_color_manual(
    name = "Market type",
    values = c("0" = "grey", "1" = "black"),
    labels = c("Investment Grade", "High Yield")
  ) +
  labs(
    title = "Vanguard's Market Share in IG vs HY Bond Holdings",
    x = "Quarter",
    y = "Market share (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



#counting cusip only once per parent
Vanguard_shares <- holdings %>%
  distinct(qt, cusip, parent_name, dol_holding) %>%   # ensure one row per bond per quarter per parent
  group_by(qt) %>%
  summarise(
    total_holding = sum(dol_holding, na.rm = TRUE),            # total across all unique bonds
    Vanguard_holding = sum(
      dol_holding[str_detect(parent_name, regex("Vanguard", ignore_case = TRUE))],
      na.rm = TRUE
    ),
    Vanguard_shares = Vanguard_holding / total_holding * 100,
    .groups = "drop"
  )

Vanguard_shares = Vanguard_shares %>%
  filter(qt >= "2002 Q1") 




ggplot(Vanguard_shares, aes(x = qt, y = Vanguard_shares, group = 1)) +
  geom_line(color = "black", size = 1.2) +
  geom_point(color = "black", size = 2) +
  labs(
    title = "Vanguard's Share of Total Holdings Over Time",
    x = "Quarter",
    y = "Share of Total Holdings (%)"
  ) +
  theme_minimal() +
  theme(
    
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey85")
  )




# ------------------------------------------------------------------------------
# SHARE OF THE MARKET WITH RESPECT TO THE BOND SIZE
# ------------------------------------------------------------------------------


holdings_with_size <- holdings %>% filter(str_detect(parent_name, regex("BlackRock", ignore_case = TRUE)) | str_detect(parent_name, regex("Fidelity", ignore_case = TRUE)) | str_detect(parent_name, regex("State Street", ignore_case = TRUE)) | str_detect(parent_name, regex("Vanguard", ignore_case = TRUE)) )


dt_final_long_boa_rs <- dt_final_long_boa_rs %>%
  left_join(
    dt_unique %>% select(complete_cusip, cusip_8),
    by = "complete_cusip"
  )


only_issuance_holdings = holdings %>%
  semi_join(
    dt_final_long_boa_rs,
    by = c("qt" = "qt_off", "cusip" = "cusip_8")
  )


only_issuance_holdings <- only_issuance_holdings%>%
  left_join(dt_final_long_boa_rs %>% select(cusip_8, complete_cusip, off_amt),
            by = c("cusip" = "cusip_8")) %>%
  mutate(
    
    holding_share = dol_holding / off_amt   # quota detenuta rispetto alla size
  )


blackrock_holdings <- only_issuance_holdings %>%
  filter(str_detect(parent_name, regex("BlackRock", ignore_case = TRUE))) %>%  # only BlackRock
  group_by(qt) %>%
  summarise(
    total_black_holding = sum(holding_share, na.rm = TRUE),   # total fraction of market held by BlackRock
    avg_black_holding   = mean(holding_share, na.rm = TRUE),  # average share per bond held
    .groups = "drop"
  )



ggplot(blackrock_holdings, aes(x = qt, y = avg_black_holding)) +
  geom_line(color = "#0072B2", size = 1.2) +
  geom_point(color = "#0072B2", size = 2) +
  labs(
    title = "BlackRock’s Share of the Corporate Bond Market",
    x = "Quarter",
    y = "Market Share (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey70"),  # lighter grid
    panel.grid.minor = element_line(color = "grey85")
  )




#compute and plot the share of the top 4 investors


top_4_holdings_HY <- holdings_HY %>%
  group_by(qt, is_high_yield) %>%
  summarise(
    total_holding = sum(dol_holding, na.rm = TRUE),
    top_4_holdings = sum(
      dol_holding[
        str_detect(parent_name, regex("BlackRock", ignore_case = TRUE)) |
          str_detect(parent_name, regex("Fidelity", ignore_case = TRUE)) |
          str_detect(parent_name, regex("State Street", ignore_case = TRUE)) |
          str_detect(parent_name, regex("Vanguard", ignore_case = TRUE))
      ],
      na.rm = TRUE
    ),
    top_4_share_HY = top_4_holdings / total_holding * 100,
    .groups = "drop"
  ) %>%
  filter(qt >= as.yearqtr("2002 Q1"))




ggplot(top_4_holdings_HY, 
       aes(x = as.Date(qt), 
           y = top_4_share_HY, 
           color = factor(is_high_yield, labels = c("Investment Grade", "High Yield")))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Top 4 Asset Managers' Market Share by Segment",
    subtitle = "Share of total holdings (%) for Investment Grade vs High Yield bonds",
    x = "Quarter",
    y = "Top 4 Holdings Share (%)",
    color = "Market Segment"
  ) +
  scale_color_manual(values = c(
    "Investment Grade" = "#1E90FF",   # electric blue
    "High Yield" = "#FF4500"          # vivid orange-red
  )) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
