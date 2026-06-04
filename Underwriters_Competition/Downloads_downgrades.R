
# ==============================================================================
# SET WORKING DIRECTORY
# ==============================================================================
path <- "C:/Users/david/Desktop/PreDoc/OTC-Markets/shared_filippo_davide/"
setwd(path)


# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

# Data Handling
library(tidyverse)
library(lubridate)
library(abind)
library(entropy)
library(forcats)
library(stringdist)

# Time Series Analysis
library(tseries)
library(zoo)

# I/O
library(writexl)
library(readxl)
library(xtable)

# Plotting
library(ggplot2)
library(forcats)
library(scales)

# conflict
library(conflicted)
conflict_prefer("select", "dplyr", quiet = TRUE)
conflict_prefer("filter", "dplyr", quiet = TRUE)

# ==============================================================================
# LOAD DATA AND FUNCTIONS
# ==============================================================================

# Load the sample of bonds (quarterly data)
load("data/def_final.RData")

load("data/dng_final.RData")

holdings_total_filled <- readRDS("data/holdings_total_filled.rds")
ls(holdings_total_filled)
# ==============================================================================
# Let’s check the holdings of institutions in the quarter of issuance (as a proxy for placement).
# Do we see an increase in mutual fund compared to insurance holdings?
# ==============================================================================

# Crea classificazione
holdings_class <- holdings_total_filled %>%
  mutate(
    inv_group = case_when(
      fund_type %in% c("BAL", "MMM", "MUT", "END", "QUI", "FOF") ~ "Mutual Funds",
      fund_type %in% c("INS", "LIN", "PIN")                      ~ "Insurance",
      fund_type %in% c("CPF", "GPE", "UPE")                      ~ "Pension Funds",
      TRUE ~ "Other"
    ),
    qt = as.yearqtr(date)  # portiamo la data a formato quarter
  )

# Aggrega per quarter e gruppo
holdings_qt <- holdings_class %>%
  group_by(qt, inv_group) %>%
  summarise(total_holding = sum(dol_holding, na.rm = TRUE), .groups = "drop")

# Percentuali (per comparare MF vs Insurance)
holdings_qt <- holdings_qt %>%
  group_by(qt) %>%
  mutate(perc = total_holding / sum(total_holding) * 100) %>%
  ungroup()

head(holdings_qt)

# Visualization 

ggplot(
  holdings_qt %>% filter(inv_group %in% c("Mutual Funds", "Insurance")),
  aes(x = qt, y = perc, color = inv_group)
) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    x = "Quarter",
    y = "Share of Holdings (%)",
    color = "Investor group",
    title = "Mutual Funds vs Insurance holdings over time"
  ) +
  theme_minimal()


# Relative to the size of the bond 

holdings_with_size <- holdings_total_filled %>%
  left_join(dt_issuance %>% select(complete_cusip, issue_id, off_amt, qt_off_yq),
            by = c("cusip" = "complete_cusip"))

holdings_with_size <- holdings_with_size %>%
  mutate(
    qt = as.yearqtr(date),
    holding_share = dol_holding / off_amt   # quota rispetto alla size bond
  )

holdings_class_rel <- holdings_with_size %>%
  mutate(
    inv_group = case_when(
      fund_type %in% c("BAL", "MMM", "MUT", "END", "QUI", "FOF") ~ "Mutual Funds",
      fund_type %in% c("INS", "LIN", "PIN")                      ~ "Insurance",
      fund_type %in% c("CPF", "GPE", "UPE")                      ~ "Pension Funds",
      TRUE ~ "Other"
    ),
    qt = as.yearqtr(date)
  ) %>%
  group_by(qt, inv_group) %>%
  summarise(
    avg_share = mean(holding_share, na.rm = TRUE),   # media quota detenuta
    total_share = sum(holding_share, na.rm = TRUE),  # somma delle quote
    .groups = "drop"
  )

ggplot(
  holdings_class_rel %>% filter(inv_group %in% c("Mutual Funds", "Insurance")),
  aes(x = qt, y = avg_share * 100, color = inv_group)
) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    x = "Quarter",
    y = "Average share of bond size held (%)",
    color = "Investor group",
    title = "Mutual Funds vs Insurance: Holdings relative to bond size"
  ) +
  theme_minimal()

# Costruisci la dummy per default

# Porta date a yearqtr
def.final <- def.final %>%
  mutate(default_qt = as.yearqtr(default_date))

dt_issuance_qt <- dt_issuance %>%
  mutate(issue_qt = as.yearqtr(qt_off_yq))

# Merge
issuance_defaults <- dt_issuance_qt %>%
  left_join(def.final %>% select(issue_id, default_qt), by = "issue_id") %>%
  mutate(default_dummy = if_else(!is.na(default_qt) & default_qt > issue_qt, 1, 0))


# Costruisci la dummy per downgrade

dng <- dng %>%
  mutate(rating_qt = as.yearqtr(rating_date))

issuance_dng <- dt_issuance_qt %>%
  left_join(dng %>% select(issue_id, rating_qt, rating), by = "issue_id") %>%
  mutate(dng_dummy = if_else(!is.na(rating_qt) & rating_qt > issue_qt, 1, 0))

# Default share per quarter
default_share <- issuance_defaults %>%
  group_by(issue_qt) %>%
  summarise(share_default = mean(default_dummy, na.rm = TRUE),
            .groups = "drop")

# Downgrade share per quarter
dng_share <- issuance_dng %>%
  group_by(issue_qt) %>%
  summarise(share_dng = mean(dng_dummy, na.rm = TRUE),
            .groups = "drop")
library(ggplot2)

ggplot(default_share, aes(x = issue_qt, y = share_default)) +
  geom_line(color = "red") +
  geom_point() +
  labs(x = "Quarter of issuance", y = "Share of defaults",
       title = "Share of bonds defaulting after issuance") +
  theme_minimal()

ggplot(dng_share, aes(x = issue_qt, y = share_dng)) +
  geom_line(color = "blue") +
  geom_point() +
  labs(x = "Quarter of issuance", y = "Share of downgrades",
       title = "Share of bonds downgraded after issuance") +
  theme_minimal()

# 1. Cumulative incidence curve (default entro N anni dall’emissione)

# Prepariamo emissioni + default
issuance_def <- dt_issuance %>%
  mutate(issue_qt = as.yearqtr(qt_off_yq)) %>%
  left_join(def.final %>% mutate(default_qt = as.yearqtr(default_date)),
            by = "issue_id") %>%
  mutate(time_to_default = as.numeric((default_qt - issue_qt) * 4))  # trimestri

# Dummy: default entro N trimestri
max_horizon <- 12  # es: 5 anni = 20 trimestri
surv_data <- issuance_def %>%
  mutate(default_within = if_else(!is.na(time_to_default) &
                                    time_to_default <= max_horizon &
                                    time_to_default >= 0, 1, 0)) %>%
  group_by(time_to_default) %>%
  summarise(cum_prob = mean(default_within, na.rm = TRUE), .groups = "drop")

# Plot cumulative incidence
ggplot(surv_data, aes(x = time_to_default/4, y = cum_prob)) +
  geom_line(color = "red", linewidth = 1.2) +
  labs(x = "Years since issuance",
       y = "Cumulative share of defaults",
       title = "Cumulative incidence of defaults after issuance") +
  theme_minimal()

# 2. HY vs IG – Share of defaults/downgrades over time
# Assumiamo che dt_issuance abbia rating_group (HY / IG)
issuance_defaults <- dt_issuance %>%
  mutate(issue_qt = as.yearqtr(qt_off_yq)) %>%
  left_join(def.final %>% mutate(default_qt = as.yearqtr(default_date)),
            by = "issue_id") %>%
  mutate(default_dummy = if_else(!is.na(default_qt) & default_qt > issue_qt, 1, 0))

default_share_rating <- issuance_defaults %>%
  group_by(issue_qt, rating_group) %>%
  summarise(share_default = mean(default_dummy, na.rm = TRUE), .groups = "drop")

ggplot(default_share_rating, aes(x = issue_qt, y = share_default, color = rating_group)) +
  geom_line(linewidth = 1) + geom_point() +
  labs(x = "Quarter of issuance", y = "Share defaults",
       title = "Share of defaults: HY vs IG") +
  theme_minimal()

# 3- New vs Old firms – Defaults
issuance_def_stage <- issuance_defaults %>%
  left_join(first_entry, by = "gvkey") %>%
  mutate(first_qt = as.yearqtr(first_qt),
         age_quarters = as.numeric((issue_qt - first_qt) * 4),
         firm_stage = if_else(age_quarters <= 20, "New", "Old")) %>%
  group_by(issue_qt, firm_stage) %>%
  summarise(share_default = mean(default_dummy, na.rm = TRUE), .groups = "drop")

ggplot(issuance_def_stage, aes(x = issue_qt, y = share_default, color = firm_stage)) +
  geom_line(linewidth = 1) + geom_point() +
  labs(x = "Quarter of issuance", y = "Share defaults",
       title = "Share of defaults: New vs Old firms") +
  theme_minimal()

# 4. Scatter Gross Spread vs Ex-post events
ggplot(issuance_defaults, aes(x = gross_spread, y = default_dummy, color = rating_group)) +
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
  labs(x = "Gross spread at issuance", y = "Default probability",
       title = "Relation between Gross Spread and ex-post defaults") +
  theme_minimal()




# Conta i default avvenuti in ciascun trimestre
default_events <- def.final %>%
  mutate(event_qt = as.yearqtr(default_date)) %>%
  group_by(event_qt) %>%
  summarise(n_default = n(), .groups = "drop")

# Conta i downgrade avvenuti in ciascun trimestre
dng_events <- dng %>%
  mutate(event_qt = as.yearqtr(rating_date)) %>%
  group_by(event_qt) %>%
  summarise(n_dng = n(), .groups = "drop")

# Merge default share con eventi
default_combined <- default_share %>%
  left_join(default_events, by = c("issue_qt" = "event_qt"))

# Merge downgrade share con eventi
dng_combined <- dng_share %>%
  left_join(dng_events, by = c("issue_qt" = "event_qt"))

library(ggplot2)

# Default: share + number of defaults
ggplot(default_combined, aes(x = issue_qt)) +
  geom_line(aes(y = share_default), color = "red", linewidth = 1) +
  geom_point(aes(y = share_default), color = "red") +
  geom_col(aes(y = n_default / max(n_default, na.rm=TRUE)), 
           fill = "grey", alpha = 0.4) +
  scale_y_continuous(
    name = "Share of defaults (emissions)",
    sec.axis = sec_axis(~ . * max(default_combined$n_default, na.rm=TRUE),
                        name = "Number of defaults (events)")
  ) +
  labs(x = "Quarter", title = "Defaults: share of emissions vs number of events") +
  theme_minimal()

# Downgrades: share + number of downgrades
ggplot(dng_combined, aes(x = issue_qt)) +
  geom_line(aes(y = share_dng), color = "blue", linewidth = 1) +
  geom_point(aes(y = share_dng), color = "blue") +
  geom_col(aes(y = n_dng / max(n_dng, na.rm=TRUE)), 
           fill = "grey", alpha = 0.4) +
  scale_y_continuous(
    name = "Share of downgrades (emissions)",
    sec.axis = sec_axis(~ . * max(dng_combined$n_dng, na.rm=TRUE),
                        name = "Number of downgrades (events)")
  ) +
  labs(x = "Quarter", title = "Downgrades: share of emissions vs number of events") +
  theme_minimal()

cap_value <- 300

default_combined <- default_combined %>%
  mutate(n_default_capped = pmin(n_default, cap_value))

dng_combined <- dng_combined %>%
  mutate(n_dng_capped = pmin(n_dng, cap_value))

ggplot(default_combined, aes(x = issue_qt)) +
  geom_line(aes(y = share_default), color = "red", linewidth = 1) +
  geom_point(aes(y = share_default), color = "red") +
  geom_col(aes(y = n_default_capped / max(n_default_capped, na.rm=TRUE)), 
           fill = "grey", alpha = 0.4) +
  scale_y_continuous(
    name = "Share of defaults (emissions)",
    sec.axis = sec_axis(~ . * max(default_combined$n_default_capped, na.rm=TRUE),
                        name = "Number of defaults (events, capped at 300)")
  ) +
  labs(x = "Quarter", title = "Defaults: share of emissions vs number of events") +
  theme_minimal()

ggplot(dng_combined, aes(x = issue_qt)) +
  geom_line(aes(y = share_dng), color = "blue", linewidth = 1) +
  geom_point(aes(y = share_dng), color = "blue") +
  geom_col(aes(y = n_dng_capped / max(n_dng_capped, na.rm=TRUE)), 
           fill = "grey", alpha = 0.4) +
  scale_y_continuous(
    name = "Share of downgrades (emissions)",
    sec.axis = sec_axis(~ . * max(dng_combined$n_dng_capped, na.rm=TRUE),
                        name = "Number of downgrades (events, capped at 300)")
  ) +
  labs(x = "Quarter", 
       title = "Downgrades: share of emissions vs number of events (capped)") +
  theme_minimal()

# Default share per rating group (ex-ante, emissioni)
default_share_rating <- issuance_defaults %>%
  group_by(issue_qt, rating_group) %>%
  summarise(share_default = mean(default_dummy, na.rm = TRUE),
            .groups = "drop")

# Conteggi di default osservati (ex-post) per rating group
default_events_rating <- def.final %>%
  left_join(dt_issuance %>% select(issue_id, rating_group), by = "issue_id") %>%
  mutate(event_qt = as.yearqtr(default_date)) %>%
  group_by(event_qt, rating_group) %>%
  summarise(n_default = n(), .groups = "drop")

# Merge
default_combined_rating <- default_share_rating %>%
  left_join(default_events_rating, by = c("issue_qt" = "event_qt", "rating_group"))
# Downgrade share per rating group (ex-ante)
dng_share_rating <- issuance_dng %>%
  group_by(issue_qt, rating_group) %>%
  summarise(share_dng = mean(dng_dummy, na.rm = TRUE),
            .groups = "drop")

# Conteggi downgrade osservati per rating group
dng_events_rating <- dng %>%
  left_join(dt_issuance %>% select(issue_id, rating_group), by = "issue_id") %>%
  mutate(event_qt = as.yearqtr(rating_date)) %>%
  group_by(event_qt, rating_group) %>%
  summarise(n_dng = n(), .groups = "drop")

# Merge
dng_combined_rating <- dng_share_rating %>%
  left_join(dng_events_rating, by = c("issue_qt" = "event_qt", "rating_group"))

# Senza capping
ggplot(default_combined_rating, aes(x = issue_qt)) +
  geom_line(aes(y = share_default, color = rating_group), linewidth = 1) +
  geom_point(aes(y = share_default, color = rating_group)) +
  geom_col(aes(y = n_default / max(n_default, na.rm=TRUE), fill = rating_group),
           alpha = 0.3, position = "dodge") +
  scale_y_continuous(
    name = "Share of defaults (emissions)",
    sec.axis = sec_axis(~ . * max(default_combined_rating$n_default, na.rm=TRUE),
                        name = "Number of defaults (events)")
  ) +
  labs(x = "Quarter", title = "Defaults: share vs number of events by Rating") +
  theme_minimal()

# Con capping
ggplot(default_combined_rating, aes(x = issue_qt)) +
  geom_line(aes(y = share_default, color = rating_group), linewidth = 1) +
  geom_point(aes(y = share_default, color = rating_group)) +
  geom_col(aes(y = n_default_capped / max(n_default_capped, na.rm=TRUE), fill = rating_group),
           alpha = 0.3, position = "dodge") +
  scale_y_continuous(
    name = "Share of defaults (emissions)",
    sec.axis = sec_axis(~ . * max(default_combined_rating$n_default_capped, na.rm=TRUE),
                        name = "Number of defaults (events, capped)")
  ) +
  labs(x = "Quarter", title = "Defaults (capped at 300): share vs number of events by Rating") +
  theme_minimal()

# Senza capping
ggplot(dng_combined_rating, aes(x = issue_qt)) +
  geom_line(aes(y = share_dng, color = rating_group), linewidth = 1) +
  geom_point(aes(y = share_dng, color = rating_group)) +
  geom_col(aes(y = n_dng / max(n_dng, na.rm=TRUE), fill = rating_group),
           alpha = 0.3, position = "dodge") +
  scale_y_continuous(
    name = "Share of downgrades (emissions)",
    sec.axis = sec_axis(~ . * max(dng_combined_rating$n_dng, na.rm=TRUE),
                        name = "Number of downgrades (events)")
  ) +
  labs(x = "Quarter", title = "Downgrades: share vs number of events by Rating") +
  theme_minimal()

# Con capping
ggplot(dng_combined_rating, aes(x = issue_qt)) +
  geom_line(aes(y = share_dng, color = rating_group), linewidth = 1) +
  geom_point(aes(y = share_dng, color = rating_group)) +
  geom_col(aes(y = n_dng_capped / max(n_dng_capped, na.rm=TRUE), fill = rating_group),
           alpha = 0.3, position = "dodge") +
  scale_y_continuous(
    name = "Share of downgrades (emissions)",
    sec.axis = sec_axis(~ . * max(dng_combined_rating$n_dng_capped, na.rm=TRUE),
                        name = "Number of downgrades (events, capped)")
  ) +
  labs(x = "Quarter", title = "Downgrades (capped at 300): share vs number of events by Rating") +
  theme_minimal()


cap_value <- 300

default_combined_rating <- default_combined_rating %>%
  mutate(n_default_capped = pmin(n_default, cap_value))

dng_combined_rating <- dng_combined_rating %>%
  mutate(n_dng_capped = pmin(n_dng, cap_value))

ggplot(default_combined_rating, aes(x = issue_qt)) +
  geom_line(aes(y = share_default, color = rating_group), linewidth = 1) +
  geom_point(aes(y = share_default, color = rating_group)) +
  geom_col(aes(y = n_default_capped / max(n_default_capped, na.rm=TRUE),
               fill = rating_group),
           alpha = 0.3, position = "dodge") +
  scale_y_continuous(
    name = "Share of defaults (emissions)",
    sec.axis = sec_axis(~ . * max(default_combined_rating$n_default_capped, na.rm=TRUE),
                        name = "Number of defaults (events, capped)")
  ) +
  labs(x = "Quarter", title = "Defaults (capped at 300): share vs number of events by Rating") +
  theme_minimal()

ggplot(dng_combined_rating, aes(x = issue_qt)) +
  geom_line(aes(y = share_dng, color = rating_group), linewidth = 1) +
  geom_point(aes(y = share_dng, color = rating_group)) +
  geom_col(aes(y = n_dng_capped / max(n_dng_capped, na.rm=TRUE),
               fill = rating_group),
           alpha = 0.3, position = "dodge") +
  scale_y_continuous(
    name = "Share of downgrades (emissions)",
    sec.axis = sec_axis(~ . * max(dng_combined_rating$n_dng_capped, na.rm=TRUE),
                        name = "Number of downgrades (events, capped)")
  ) +
  labs(x = "Quarter", title = "Downgrades (capped at 300): share vs number of events by Rating") +
  theme_minimal()


# ==============================================================================
# ANALISI DEL RUOLO DEI MUTUAL FUNDS NEL MERCATO OBBLIGAZIONARIO
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Classificazione per tipo di investitore e conversione date in trimestri
# ------------------------------------------------------------------------------

holdings_class <- holdings_total_filled %>%
  mutate(
    # Creiamo una variabile "inv_group" che identifica il tipo di investitore
    inv_group = case_when(
      fund_type %in% c("BAL", "MMM", "MUT", "END", "QUI", "FOF") ~ "Mutual Funds",
      fund_type %in% c("INS", "LIN", "PIN")                      ~ "Insurance",
      fund_type %in% c("CPF", "GPE", "UPE")                      ~ "Pension Funds",
      TRUE ~ "Other"
    ),
    # Portiamo la data a formato "quarter"
    qt = as.yearqtr(date)
  )

# ------------------------------------------------------------------------------
# 2. Quota complessiva detenuta dai Mutual Funds nel tempo
# ------------------------------------------------------------------------------

# Aggrega per trimestre e gruppo
holdings_qt <- holdings_class %>%
  group_by(qt, inv_group) %>%
  summarise(total_holding = sum(dol_holding, na.rm = TRUE), .groups = "drop")

# Calcola la quota percentuale per ogni gruppo all’interno del trimestre
holdings_qt <- holdings_qt %>%
  group_by(qt) %>%
  mutate(perc = total_holding / sum(total_holding) * 100) %>%
  ungroup()

# Filtra solo i Mutual Funds
mf_share_qt <- holdings_qt %>% filter(inv_group == "Mutual Funds")

# Grafico: evoluzione temporale della quota di mercato dei Mutual Funds
ggplot(mf_share_qt, aes(x = qt, y = perc)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "darkred") +
  labs(
    title = "Quota di mercato dei Mutual Funds nel tempo",
    x = "Quarter",
    y = "Share of Holdings (%)"
  ) +
  theme_minimal()

# ------------------------------------------------------------------------------
# 3. Detenzioni dei Mutual Funds relative alla dimensione dei bond
# ------------------------------------------------------------------------------

# Aggiunge la size dei bond (off_amt) alle holdings
holdings_with_size <- holdings_total_filled %>%
  left_join(dt_issuance %>% select(complete_cusip, cusip_8, issue_id, off_amt, qt_off_yq),
            by = c("cusip" = "cusip_8")) %>%
  mutate(
    qt = as.yearqtr(date),
    holding_share = dol_holding / off_amt   # quota detenuta rispetto alla size
  )

# Calcola media e somma per trimestre solo per i Mutual Funds
mf_holdings_rel <- holdings_with_size %>%
  mutate(
    inv_group = case_when(
      fund_type %in% c("BAL", "MMM", "MUT", "END", "QUI", "FOF") ~ "Mutual Funds",
      TRUE ~ "Other"
    ),
    qt = as.yearqtr(date)
  ) %>%
  filter(inv_group == "Mutual Funds") %>%
  group_by(qt) %>%
  summarise(
    avg_share = mean(holding_share, na.rm = TRUE),   # quota media detenuta
    total_share = sum(holding_share, na.rm = TRUE),  # quota totale detenuta
    .groups = "drop"
  )

# Grafico: quota media detenuta rispetto alla dimensione del bond
ggplot(mf_holdings_rel, aes(x = qt, y = avg_share * 100)) +
  geom_line(color = "firebrick", linewidth = 1.2) +
  geom_point(color = "firebrick") +
  labs(
    title = "Quota media dei bond detenuta dai Mutual Funds",
    x = "Quarter",
    y = "Average share of bond size held (%)"
  ) +
  theme_minimal()

# ------------------------------------------------------------------------------
# 4. Sintesi tabellare
# ------------------------------------------------------------------------------

# Crea tabella riassuntiva per esportazione o consultazione
mf_summary <- mf_share_qt %>%
  left_join(mf_holdings_rel, by = "qt") %>%
  select(qt, perc, avg_share, total_share)

# Stampa a schermo
print(head(mf_summary))

# (Opzionale) Esporta in Excel per analisi successive
# write_xlsx(mf_summary, "output/mutual_funds_summary.xlsx")

# ------------------------------------------------------------------------------
# OUTPUT
# ------------------------------------------------------------------------------
#  - Grafico 1: quota di mercato dei Mutual Funds nel tempo
#  - Grafico 2: quota media detenuta rispetto alla size dei bond
#  - Tabella mf_summary: contiene per ogni trimestre:
#       • perc = quota % sul totale delle holdings
#       • avg_share = quota media di ciascun bond
#       • total_share = quota aggregata rispetto alla size complessiva
