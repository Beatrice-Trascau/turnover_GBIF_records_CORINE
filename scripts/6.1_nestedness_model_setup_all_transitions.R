##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 6.1_nestedness_model_setup_all_transitions
# This script contains code which sets up GLS models exploring the impact of 
# land cover transitions on beta_jne for:
#   - 3 land-cover change types: Forest→TWS, TWS→Forest, All→Urban
#   - 3 taxonomic groups: All occurrences, Vascular plants, Birds
##----------------------------------------------------------------------------##

library(here)
source(here("scripts", "0_setup.R"))

# 1. LOAD DATA -----------------------------------------------------------------

# Turnover data for all occurrences
load(here("data", "derived_data", 
          "all_periods_turnover_all_land_cover_climate_15km.rda"))

# Turnover data for vascular plants
load(here("data", "derived_data", 
          "vascular_plants_turnover_all_land_cover_climate_15km.rda"))

# Turnover data for birds
load(here("data", "derived_data", 
          "bird_turnover_all_land_cover_climate_15km.rda"))


# 2. FOREST -> TWS NESTEDNESS (BETA_JNE) ---------------------------------------

## 2.1. All occurrences --------------------------------------------------------

# Select only Forest -> TWS columns
turnover_forest_tws_15km_all <- all_periods_turnover_with_climate |>
  select(-c('2000-2006_TWS no change', '2000-2006_TWS to Forest',
            '2000-2006_Urban_no_change', '2000-2006_all_to_urban',
            '2006-2012_TWS no change', '2006-2012_TWS to Forest',
            '2006-2012_Urban_no_change', '2006-2012_all_to_urban',
            '2012-2018_TWS no change', '2012-2018_TWS to Forest',
            '2012-2018_Urban_no_change', '2012-2018_all_to_urban')) |>
  mutate(forest_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest no change`,
                                      lc_time_period == "2006-2012" ~ `2006-2012_Forest no change`,
                                      lc_time_period == "2012-2018" ~ `2012-2018_Forest no change`,
                                      TRUE ~ NA_real_),
    forest_to_tws = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest to TWS`,
                              lc_time_period == "2006-2012" ~ `2006-2012_Forest to TWS`,
                              lc_time_period == "2012-2018" ~ `2012-2018_Forest to TWS`,
                              TRUE ~ NA_real_)) |>
  select(-`2000-2006_Forest no change`, -`2006-2012_Forest no change`,
         -`2012-2018_Forest no change`, -`2000-2006_Forest to TWS`,
         -`2006-2012_Forest to TWS`, -`2012-2018_Forest to TWS`) |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort))

# Fit GLS with logged recorder effort
all_FTWS_beta_jne_model1 <- gls(beta_jne ~ forest_to_tws + forest_no_change + delta_recorder_effort + 
                                 log_recorder_effort + lc_time_period + temp_change + precip_change,
                               correlation = corExp(form = ~ x + y | time_numeric),
                               data = turnover_forest_tws_15km_all,
                               method = "REML")

# Save model
save(all_FTWS_beta_jne_model1,
     file = here("data", "models", "final", "all_FTWS_beta_jne_model1.RData"))

# Fit GLS with interaction term
all_FTWS_beta_jne_model2_interaction <-  gls(beta_jne ~ forest_to_tws * temp_change +
                                               forest_to_tws * precip_change + 
                                               forest_no_change +
                                               delta_recorder_effort + 
                                               log_recorder_effort + 
                                               lc_time_period,
                                             correlation = corExp(form = ~ x + y | time_numeric),
                                             data = turnover_forest_tws_15km_coords_time,
                                             method = "REML")


# Save model
save(all_FTWS_beta_jne_model2_interaction,
     file = here("data", "models", "exploratory", "all_FTWS_beta_jne_model2_interaction.RData"))

# Compare AICs between models
AICtab(all_FTWS_beta_jne_model1, all_FTWS_beta_jne_model2_interaction, base = TRUE)

# Get summary of final model
model_summary <- summary(all_FTWS_beta_jne_model1)

# Create dataframe of coefficients
coef_df_all_FTWS <- data.frame(term = names(model_summary$tTable[, "Value"]),
                      estimate = model_summary$tTable[, "Value"],
                      std.error = model_summary$tTable[, "Std.Error"],
                      statistic = model_summary$tTable[, "t-value"],
                      p.value = model_summary$tTable[, "p-value"])

# Add significance stars
coef_df_all_FTWS <- coef_df_all_FTWS |>
  mutate(significance = case_when( p.value < 0.001 ~ "***",
                                   p.value < 0.01 ~ "**",
                                   p.value < 0.05 ~ "*",
                                   p.value < 0.1 ~ ".",
                                   TRUE ~ ""))

# Set into nice table format
coef_df_all_FTWS_table <- coef_df_all_FTWS |>
  select(term, estimate, std.error, statistic, p.value, significance) |>
  flextable() |>
  set_header_labels(term = "Predictor", estimate = "Estimate", std.error = "SE",
                    statistic = "t-value", p.value = "p-value", significance = "") |>
  colformat_double(j = c("estimate", "std.error", "statistic"), digits = 3) |>
  colformat_double(j = "p.value", digits = 4) |>
  autofit()

# Display table
coef_df_all_FTWS_table

## 2.2. Vascular Plants --------------------------------------------------------

# Select only Forest -> TWS columns
turnover_forest_tws_15km_all <- vascular_plants_turnover_with_climate |>
  select(-c('2000-2006_TWS no change', '2000-2006_TWS to Forest',
            '2000-2006_Urban_no_change', '2000-2006_all_to_urban',
            '2006-2012_TWS no change', '2006-2012_TWS to Forest',
            '2006-2012_Urban_no_change', '2006-2012_all_to_urban',
            '2012-2018_TWS no change', '2012-2018_TWS to Forest',
            '2012-2018_Urban_no_change', '2012-2018_all_to_urban')) |>
  mutate(forest_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest no change`,
                                      lc_time_period == "2006-2012" ~ `2006-2012_Forest no change`,
                                      lc_time_period == "2012-2018" ~ `2012-2018_Forest no change`,
                                      TRUE ~ NA_real_),
         forest_to_tws = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest to TWS`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_Forest to TWS`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_Forest to TWS`,
                                   TRUE ~ NA_real_)) |>
  select(-`2000-2006_Forest no change`, -`2006-2012_Forest no change`,
         -`2012-2018_Forest no change`, -`2000-2006_Forest to TWS`,
         -`2006-2012_Forest to TWS`, -`2012-2018_Forest to TWS`) |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort))

# Fit GLS without interaction
plants_FTWS_beta_jne_model1 <- gls(beta_jne ~ forest_to_tws + forest_no_change +
                                    delta_recorder_effort + log_recorder_effort 
                                  + lc_time_period + temp_change + precip_change,
                                  correlation = corExp(form = ~ x + y | time_numeric),
                                  data = turnover_forest_tws_15km_plants,
                                  method = "REML")

# Save model
save(plants_FTWS_beta_jne_model1,
     file = here("data", "models", "final", "plants_FTWS_beta_jne_model1.RData"))

# Fit GLS with interaction term
plants_FTWS_beta_jne_model2_interaction <-  gls(beta_jne ~ forest_to_tws * temp_change +
                                               forest_to_tws * precip_change + 
                                               forest_no_change +
                                               delta_recorder_effort + 
                                               log_recorder_effort + 
                                               lc_time_period,
                                             correlation = corExp(form = ~ x + y | time_numeric),
                                             data = turnover_forest_tws_15km_plants,
                                             method = "REML")


# Save model
save(plants_FTWS_beta_jne_model2_interaction,
     file = here("data", "models", "exploratory", "plants_FTWS_beta_jne_model2_interaction.RData"))

# Compare models based on AIC
AICtab(plants_FTWS_beta_jne_model1, plants_FTWS_beta_jne_model2_interaction, base = TRUE)

# Get summary of final model
model_summary_plants <- summary(plants_FTWS_beta_jne_model1)

# Create dataframe of coefficients
coef_df_plants_FTWS <- data.frame(term = names(model_summary_plants$tTable[, "Value"]),
                               estimate = model_summary_plants$tTable[, "Value"],
                               std.error = model_summary_plants$tTable[, "Std.Error"],
                               statistic = model_summary_plants$tTable[, "t-value"],
                               p.value = model_summary_plants$tTable[, "p-value"])

# Add significance stars
coef_df_plants_FTWS <- coef_df_plants_FTWS |>
  mutate(significance = case_when( p.value < 0.001 ~ "***",
                                   p.value < 0.01 ~ "**",
                                   p.value < 0.05 ~ "*",
                                   p.value < 0.1 ~ ".",
                                   TRUE ~ ""))

# Set into nice table format
coef_df_plants_FTWS_table <- coef_df_plants_FTWS |>
  select(term, estimate, std.error, statistic, p.value, significance) |>
  flextable() |>
  set_header_labels(term = "Predictor", estimate = "Estimate", std.error = "SE",
                    statistic = "t-value", p.value = "p-value", significance = "") |>
  colformat_double(j = c("estimate", "std.error", "statistic"), digits = 3) |>
  colformat_double(j = "p.value", digits = 4) |>
  autofit()

# Display table
coef_df_plants_FTWS_table








