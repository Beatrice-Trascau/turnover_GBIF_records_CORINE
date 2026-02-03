##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 6.2_total_dissimilarity_model_setup_all_transitions
# This script contains code which sets up GLS models exploring the impact of 
# land cover transitions on beta_jac for:
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

# Total number of 100m x 100m pixels in a 15km x 15km cell
# 15,000m / 100m = 150 pixels per side
# 150 x 150 = 22,500 total pixels per cell
total_pixels_per_cell <- 22500

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
         log_recorder_effort = log(recorder_effort),
         forest_no_change_prop = forest_no_change / total_pixels_per_cell,
         forest_to_tws_prop = forest_to_tws / total_pixels_per_cell)

# Fit GLS with logged recorder effort
all_FTWS_beta_JAC_model1 <- gls(beta_jac ~ forest_to_tws_prop + forest_no_change_prop + delta_recorder_effort + 
                                  log_recorder_effort + lc_time_period + temp_change + precip_change,
                                correlation = corExp(form = ~ x + y | time_numeric),
                                data = turnover_forest_tws_15km_all,
                                method = "REML")

# Save model
save(all_FTWS_beta_JAC_model1,
     file = here("data", "models", "final", "all_FTWS_beta_jac_model1.RData"))

# Fit GLS with interaction term
all_FTWS_beta_JAC_model2_interaction <-  gls(beta_jac ~ forest_to_tws_prop * temp_change +
                                               forest_to_tws_prop * precip_change + 
                                               forest_no_change_prop +
                                               delta_recorder_effort + 
                                               log_recorder_effort + 
                                               lc_time_period,
                                             correlation = corExp(form = ~ x + y | time_numeric),
                                             data = turnover_forest_tws_15km_all,
                                             method = "REML")


# Save model
save(all_FTWS_beta_JAC_model2_interaction,
     file = here("data", "models", "exploratory", "all_FTWS_beta_JAC_model2_interaction.RData"))

# Compare AICs between models
AICtab(all_FTWS_beta_JAC_model1, all_FTWS_beta_JAC_model2_interaction, base = TRUE)

## 2.2. Vascular Plants --------------------------------------------------------

# Select only Forest -> TWS columns
turnover_forest_tws_15km_plants <- vascular_plants_turnover_with_climate |>
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
         log_recorder_effort = log(recorder_effort),
         forest_no_change_prop = forest_no_change / total_pixels_per_cell,
         forest_to_tws_prop = forest_to_tws / total_pixels_per_cell)

# Fit GLS without interaction
plants_FTWS_beta_JAC_model1 <- gls(beta_jac ~ forest_to_tws_prop + forest_no_change_prop +
                                     delta_recorder_effort + log_recorder_effort 
                                   + lc_time_period + temp_change + precip_change,
                                   correlation = corExp(form = ~ x + y | time_numeric),
                                   data = turnover_forest_tws_15km_plants,
                                   method = "REML")

# Save model
save(plants_FTWS_beta_JAC_model1,
     file = here("data", "models", "final", "plants_FTWS_beta_JAC_model1.RData"))

# Fit GLS with interaction term
plants_FTWS_beta_JAC_model2_interaction <-  gls(beta_jac ~ forest_to_tws_prop * temp_change +
                                                  forest_to_tws_prop * precip_change + 
                                                  forest_no_change_prop +
                                                  delta_recorder_effort + 
                                                  log_recorder_effort + 
                                                  lc_time_period,
                                                correlation = corExp(form = ~ x + y | time_numeric),
                                                data = turnover_forest_tws_15km_plants,
                                                method = "REML")

# Save model
save(plants_FTWS_beta_JAC_model2_interaction,
     file = here("data", "models", "exploratory", "plants_FTWS_beta_JAC_model2_interaction.RData"))

# Compare models based on AIC
AICtab(plants_FTWS_beta_JAC_model1, plants_FTWS_beta_JAC_model2_interaction, base = TRUE)

## 2.3. Birds ------------------------------------------------------------------

# Select only Forest -> TWS columns
turnover_forest_tws_15km_birds <- birds_turnover_with_climate |>
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
         log_recorder_effort = log(recorder_effort),
         forest_no_change_prop = forest_no_change / total_pixels_per_cell,
         forest_to_tws_prop = forest_to_tws / total_pixels_per_cell)

# Fit GLS without interaction
birds_FTWS_beta_JAC_model1 <- gls(beta_jac ~ forest_to_tws_prop + forest_no_change_prop + 
                                    delta_recorder_effort + log_recorder_effort +
                                    lc_time_period + temp_change + precip_change,
                                  correlation = corExp(form = ~ x + y | time_numeric),
                                  data = turnover_forest_tws_15km_birds,
                                  method = "REML")

# Save model
save(birds_FTWS_beta_JAC_model1,
     file = here("data", "models", "final", "birds_FTWS_beta_JAC_model1.RData"))


# Fit GLS with interaction
birds_FTWS_beta_JAC_model2_interaction <-  gls(beta_jac ~ forest_to_tws_prop * temp_change +
                                                 forest_to_tws_prop * precip_change + 
                                                 forest_no_change_prop +
                                                 delta_recorder_effort + 
                                                 log_recorder_effort + 
                                                 lc_time_period,
                                               correlation = corExp(form = ~ x + y | time_numeric),
                                               data = turnover_forest_tws_15km_birds,
                                               method = "REML")


# Save model
save(birds_FTWS_beta_JAC_model2_interaction,
     file = here("data", "models", "exploratory", "birds_FTWS_beta_JAC_model2_interaction.RData"))

# Compare models based on AIC
AICtab(birds_FTWS_beta_JAC_model1, birds_FTWS_beta_JAC_model2_interaction, base = TRUE)

# 3. TWS -> FOREST NESTEDNESS (BETA_JNE) ---------------------------------------

## 3.1. All occurrences --------------------------------------------------------

# Select only TWS -> Forest columns
turnover_tws_forest_15km_all <- all_periods_turnover_with_climate |>
  select(-c('2000-2006_Forest no change', '2000-2006_Forest to TWS',
            '2000-2006_Urban_no_change', '2000-2006_all_to_urban',
            '2006-2012_Forest no change', '2006-2012_Forest to TWS',
            '2006-2012_Urban_no_change', '2006-2012_all_to_urban',
            '2012-2018_Forest no change', '2012-2018_Forest to TWS',
            '2012-2018_Urban_no_change', '2012-2018_all_to_urban')) |>
  mutate(tws_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS no change`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_TWS no change`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_TWS no change`,
                                   TRUE ~ NA_real_),
         tws_to_forest = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS to Forest`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_TWS to Forest`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_TWS to Forest`,
                                   TRUE ~ NA_real_)) |>
  select(-`2000-2006_TWS no change`, -`2006-2012_TWS no change`,
         -`2012-2018_TWS no change`, -`2000-2006_TWS to Forest`,
         -`2006-2012_TWS to Forest`, -`2012-2018_TWS to Forest`) |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort),
         tws_no_change_prop = tws_no_change / total_pixels_per_cell,
         tws_to_forest_prop = tws_to_forest / total_pixels_per_cell)

# Fit GLS without interaction 
all_TWSF_beta_JAC_model1 <- gls(beta_jac ~ tws_to_forest_prop + tws_no_change_prop + 
                                  delta_recorder_effort + log_recorder_effort + 
                                  lc_time_period + temp_change + precip_change,
                                correlation = corExp(form = ~ x + y | time_numeric),
                                data = turnover_tws_forest_15km_all,
                                method = "REML")

# Save model
save(all_TWSF_beta_JAC_model1,
     file = here("data", "models", "final", "all_TWSF_beta_JAC_model1.RData"))

# Fit GLS with interaction term
all_TWSF_beta_JAC_model2_interaction <-  gls(beta_jac ~ tws_to_forest_prop * temp_change +
                                               tws_to_forest_prop * precip_change + 
                                               tws_no_change_prop +
                                               delta_recorder_effort + 
                                               log_recorder_effort + 
                                               lc_time_period,
                                             correlation = corExp(form = ~ x + y | time_numeric),
                                             data = turnover_tws_forest_15km_all,
                                             method = "REML")


# Save model
save(all_TWSF_beta_JAC_model2_interaction,
     file = here("data", "models", "exploratory", "all_TWSF_beta_JAC_model2_interaction.RData"))

# Compare AICs between models
AICtab(all_TWSF_beta_JAC_model1, all_TWSF_beta_JAC_model2_interaction, base = TRUE)

## 3.2. Vascular Plants --------------------------------------------------------

# Select only TWS -> Forest columns
turnover_tws_forest_15km_plants <- vascular_plants_turnover_with_climate |>
  select(-c('2000-2006_Forest no change', '2000-2006_Forest to TWS',
            '2000-2006_Urban_no_change', '2000-2006_all_to_urban',
            '2006-2012_Forest no change', '2006-2012_Forest to TWS',
            '2006-2012_Urban_no_change', '2006-2012_all_to_urban',
            '2012-2018_Forest no change', '2012-2018_Forest to TWS',
            '2012-2018_Urban_no_change', '2012-2018_all_to_urban')) |>
  mutate(tws_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS no change`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_TWS no change`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_TWS no change`,
                                   TRUE ~ NA_real_),
         tws_to_forest = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS to Forest`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_TWS to Forest`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_TWS to Forest`,
                                   TRUE ~ NA_real_)) |>
  select(-`2000-2006_TWS no change`, -`2006-2012_TWS no change`,
         -`2012-2018_TWS no change`, -`2000-2006_TWS to Forest`,
         -`2006-2012_TWS to Forest`, -`2012-2018_TWS to Forest`) |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort),
         tws_no_change_prop = tws_no_change / total_pixels_per_cell,
         tws_to_forest_prop = tws_to_forest / total_pixels_per_cell)


# Fit GLS without interaction
plants_TWSF_beta_JAC_model1 <- gls(beta_jac ~ tws_to_forest_prop + tws_no_change_prop + 
                                     delta_recorder_effort + log_recorder_effort +
                                     lc_time_period + temp_change + precip_change,
                                   correlation = corExp(form = ~ x + y | time_numeric),
                                   data = turnover_tws_forest_15km_plants,
                                   method = "REML")

# Save model output
save(plants_TWSF_beta_JAC_model1,
     file = here("data", "models", "final", "plants_TWSF_beta_JAC_model1.RData"))

# Fit GLS with interaction term
plants_TWSF_beta_JAC_model2_interaction <-  gls(beta_jac ~ tws_to_forest_prop * temp_change +
                                                  tws_to_forest_prop * precip_change + 
                                                  tws_no_change_prop +
                                                  delta_recorder_effort + 
                                                  log_recorder_effort + 
                                                  lc_time_period,
                                                correlation = corExp(form = ~ x + y | time_numeric),
                                                data = turnover_tws_forest_15km_plants,
                                                method = "REML")

# Save model
save(plants_TWSF_beta_JAC_model2_interaction,
     file = here("data", "models", "exploratory", "plants_TWSF_beta_JAC_model2_interaction.RData"))

# Compare models based on AIC
AICtab(plants_TWSF_beta_JAC_model1, plants_TWSF_beta_JAC_model2_interaction, base = TRUE)

## 3.3. Birds ------------------------------------------------------------------

# Select only TWS -> Forest columns
turnover_tws_forest_15km_birds <- birds_turnover_with_climate |>
  select(-c('2000-2006_Forest no change', '2000-2006_Forest to TWS',
            '2000-2006_Urban_no_change', '2000-2006_all_to_urban',
            '2006-2012_Forest no change', '2006-2012_Forest to TWS',
            '2006-2012_Urban_no_change', '2006-2012_all_to_urban',
            '2012-2018_Forest no change', '2012-2018_Forest to TWS',
            '2012-2018_Urban_no_change', '2012-2018_all_to_urban')) |>
  mutate(tws_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS no change`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_TWS no change`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_TWS no change`,
                                   TRUE ~ NA_real_),
         tws_to_forest = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS to Forest`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_TWS to Forest`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_TWS to Forest`,
                                   TRUE ~ NA_real_)) |>
  select(-`2000-2006_TWS no change`, -`2006-2012_TWS no change`,
         -`2012-2018_TWS no change`, -`2000-2006_TWS to Forest`,
         -`2006-2012_TWS to Forest`, -`2012-2018_TWS to Forest`) |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort),
         tws_no_change_prop = tws_no_change / total_pixels_per_cell,
         tws_to_forest_prop = tws_to_forest / total_pixels_per_cell)


# Fit GLS without interaction
birds_TWSF_beta_JAC_model1 <- gls(beta_jac ~ tws_to_forest_prop + tws_no_change_prop + 
                                    delta_recorder_effort + log_recorder_effort +
                                    lc_time_period + temp_change + precip_change,
                                  correlation = corExp(form = ~ x + y | time_numeric),
                                  data = turnover_tws_forest_15km_birds,
                                  method = "REML")

# Save model output
save(birds_TWSF_beta_JAC_model1,
     file = here("data", "models", "final", "birds_TWSF_beta_JAC_model1.RData"))

# Fit GLS with interaction
birds_TWSF_beta_JAC_model2_interaction <-  gls(beta_jac ~ tws_to_forest_prop * temp_change +
                                                 tws_to_forest_prop * precip_change + 
                                                 tws_no_change_prop +
                                                 delta_recorder_effort + 
                                                 log_recorder_effort + 
                                                 lc_time_period,
                                               correlation = corExp(form = ~ x + y | time_numeric),
                                               data = turnover_tws_forest_15km_birds,
                                               method = "REML")


# Save model
save(birds_TWSF_beta_JAC_model2_interaction,
     file = here("data", "models", "exploratory", "birds_TWSF_beta_JAC_model2_interaction.RData"))

# Compare models based on AIC
AICtab(birds_TWSF_beta_JAC_model1, birds_TWSF_beta_JAC_model2_interaction, base = TRUE)

# 4. ALL -> URBAN NESTEDNESS (BETA_JNE) ----------------------------------------

## 4.1. All occurrences --------------------------------------------------------

# Select only All -> Urban columns
turnover_all_urban_15km_all <- all_periods_turnover_with_climate |>
  select(-c('2000-2006_Forest no change', '2000-2006_Forest to TWS',
            '2000-2006_TWS no change', '2000-2006_TWS to Forest',
            '2006-2012_Forest no change', '2006-2012_Forest to TWS',
            '2006-2012_TWS no change', '2006-2012_TWS to Forest',
            '2012-2018_Forest no change', '2012-2018_Forest to TWS',
            '2012-2018_TWS no change', '2012-2018_TWS to Forest')) |>
  mutate(urban_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Urban_no_change`,
                                     lc_time_period == "2006-2012" ~ `2006-2012_Urban_no_change`,
                                     lc_time_period == "2012-2018" ~ `2012-2018_Urban_no_change`,
                                     TRUE ~ NA_real_),
         all_to_urban = case_when(lc_time_period == "2000-2006" ~ `2000-2006_all_to_urban`,
                                  lc_time_period == "2006-2012" ~ `2006-2012_all_to_urban`,
                                  lc_time_period == "2012-2018" ~ `2012-2018_all_to_urban`,
                                  TRUE ~ NA_real_)) |>
  select(-`2000-2006_Urban_no_change`, -`2000-2006_all_to_urban`,
         -`2006-2012_Urban_no_change`, -`2006-2012_all_to_urban`,
         -`2012-2018_Urban_no_change`, -`2012-2018_all_to_urban`) |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort),
         urban_no_change_prop = urban_no_change / total_pixels_per_cell,
         all_to_urban_prop = all_to_urban / total_pixels_per_cell)

# Fit GLS without interaction
all_urban_beta_JAC_model1 <- gls(beta_jac ~ all_to_urban_prop + urban_no_change_prop +
                                   delta_recorder_effort + log_recorder_effort + 
                                   lc_time_period + temp_change + precip_change,
                                 correlation = corExp(form = ~ x + y | time_numeric),
                                 data = turnover_all_urban_15km_all,
                                 method = "REML")

# Save model
save(all_urban_beta_JAC_model1,
     file = here("data", "models", "final", "all_urban_beta_JAC_model1.RData"))

# Fit GLS with interaction 
all_urban_beta_JAC_model2_interaction <-  gls(beta_jac ~ all_to_urban_prop * temp_change +
                                                all_to_urban_prop * precip_change + 
                                                urban_no_change_prop +
                                                delta_recorder_effort + 
                                                log_recorder_effort + 
                                                lc_time_period,
                                              correlation = corExp(form = ~ x + y | time_numeric),
                                              data = turnover_all_urban_15km_all,
                                              method = "REML")


# Save model
save(all_urban_beta_JAC_model2_interaction,
     file = here("data", "models", "exploratory", "all_urban_beta_JAC_model2_interaction.RData"))

# Compare AICs between models
AICtab(all_urban_beta_JAC_model1, all_urban_beta_JAC_model2_interaction, base = TRUE)

## 4.2. Vascular Plants --------------------------------------------------------

# Select only All -> Urban columns
turnover_all_urban_15km_plants <- vascular_plants_turnover_with_climate |>
  select(-c('2000-2006_Forest no change', '2000-2006_Forest to TWS',
            '2000-2006_TWS no change', '2000-2006_TWS to Forest',
            '2006-2012_Forest no change', '2006-2012_Forest to TWS',
            '2006-2012_TWS no change', '2006-2012_TWS to Forest',
            '2012-2018_Forest no change', '2012-2018_Forest to TWS',
            '2012-2018_TWS no change', '2012-2018_TWS to Forest')) |>
  mutate(urban_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Urban_no_change`,
                                     lc_time_period == "2006-2012" ~ `2006-2012_Urban_no_change`,
                                     lc_time_period == "2012-2018" ~ `2012-2018_Urban_no_change`,
                                     TRUE ~ NA_real_),
         all_to_urban = case_when(lc_time_period == "2000-2006" ~ `2000-2006_all_to_urban`,
                                  lc_time_period == "2006-2012" ~ `2006-2012_all_to_urban`,
                                  lc_time_period == "2012-2018" ~ `2012-2018_all_to_urban`,
                                  TRUE ~ NA_real_)) |>
  select(-`2000-2006_Urban_no_change`, -`2000-2006_all_to_urban`,
         -`2006-2012_Urban_no_change`, -`2006-2012_all_to_urban`,
         -`2012-2018_Urban_no_change`, -`2012-2018_all_to_urban`) |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort),
         urban_no_change_prop = urban_no_change / total_pixels_per_cell,
         all_to_urban_prop = all_to_urban / total_pixels_per_cell)

# Fit GLS without interaction
plants_urban_beta_JAC_model1 <- gls(beta_jac ~ all_to_urban_prop + urban_no_change_prop +
                                      delta_recorder_effort + log_recorder_effort + 
                                      lc_time_period + temp_change + precip_change,
                                    correlation = corExp(form = ~ x + y | time_numeric),
                                    data = turnover_all_urban_15km_plants,
                                    method = "REML")

# Save model 
save(plants_urban_beta_JAC_model1,
     file = here("data", "models", "final", "plants_urban_beta_JAC_model1.RData"))

# Fit GLS with interaction term
plants_urban_beta_JAC_model2_interaction <-  gls(beta_jac ~ all_to_urban_prop * temp_change +
                                                   all_to_urban_prop * precip_change + 
                                                   urban_no_change_prop +
                                                   delta_recorder_effort + 
                                                   log_recorder_effort + 
                                                   lc_time_period,
                                                 correlation = corExp(form = ~ x + y | time_numeric),
                                                 data = turnover_all_urban_15km_plants,
                                                 method = "REML")

# Save model
save(plants_urban_beta_JAC_model2_interaction,
     file = here("data", "models", "exploratory", "plants_urban_beta_JAC_model2_interaction.RData"))

# Compare models based on AIC
AICtab(plants_urban_beta_JAC_model1, plants_urban_beta_JAC_model2_interaction, base = TRUE)

## 4.3. Birds ------------------------------------------------------------------

# Select only All -> Urban columns
turnover_all_urban_15km_birds <- birds_turnover_with_climate |>
  select(-c('2000-2006_Forest no change', '2000-2006_Forest to TWS',
            '2000-2006_TWS no change', '2000-2006_TWS to Forest',
            '2006-2012_Forest no change', '2006-2012_Forest to TWS',
            '2006-2012_TWS no change', '2006-2012_TWS to Forest',
            '2012-2018_Forest no change', '2012-2018_Forest to TWS',
            '2012-2018_TWS no change', '2012-2018_TWS to Forest')) |>
  mutate(urban_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Urban_no_change`,
                                     lc_time_period == "2006-2012" ~ `2006-2012_Urban_no_change`,
                                     lc_time_period == "2012-2018" ~ `2012-2018_Urban_no_change`,
                                     TRUE ~ NA_real_),
         all_to_urban = case_when(lc_time_period == "2000-2006" ~ `2000-2006_all_to_urban`,
                                  lc_time_period == "2006-2012" ~ `2006-2012_all_to_urban`,
                                  lc_time_period == "2012-2018" ~ `2012-2018_all_to_urban`,
                                  TRUE ~ NA_real_)) |>
  select(-`2000-2006_Urban_no_change`, -`2000-2006_all_to_urban`,
         -`2006-2012_Urban_no_change`, -`2006-2012_all_to_urban`,
         -`2012-2018_Urban_no_change`, -`2012-2018_all_to_urban`) |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort),
         urban_no_change_prop = urban_no_change / total_pixels_per_cell,
         all_to_urban_prop = all_to_urban / total_pixels_per_cell)

# Fit GLS model without interaction
birds_urban_beta_JAC_model1 <- gls(beta_jac ~ all_to_urban_prop + urban_no_change_prop + 
                                     delta_recorder_effort + log_recorder_effort + 
                                     lc_time_period + temp_change + precip_change,
                                   correlation = corExp(form = ~ x + y | time_numeric),
                                   data = turnover_all_urban_15km_birds,
                                   method = "REML")

# Save model
save(birds_urban_beta_JAC_model1,
     file = here("data", "models", "final", "birds_urban_beta_JAC_model1.RData"))

# Fit GLS with interaction
birds_urban_beta_JAC_model2_interaction <-  gls(beta_jac ~ all_to_urban_prop * temp_change +
                                                  all_to_urban_prop * precip_change + 
                                                  urban_no_change_prop +
                                                  delta_recorder_effort + 
                                                  log_recorder_effort + 
                                                  lc_time_period,
                                                correlation = corExp(form = ~ x + y | time_numeric),
                                                data = turnover_all_urban_15km_birds,
                                                method = "REML")


# Save model
save(birds_urban_beta_JAC_model2_interaction,
     file = here("data", "models", "exploratory", "birds_urban_beta_JAC_model2_interaction.RData"))

# Compare models based on AIC
AICtab(birds_urban_beta_JAC_model1, birds_urban_beta_JAC_model2_interaction, base = TRUE)

# 5. CREATE COMBINED COEFFICIENT PLOT ------------------------------------------

# Function to create coefficient plot for a single model
create_coef_plot_beta_jac <- function(model, land_cover_type, show_y_axis = TRUE) {
  
  # Get model summary
  model_summary <- summary(model)
  
  # Create a dataframe of the coefficients
  coef_df <- data.frame(
    term = names(model_summary$tTable[, "Value"]),
    estimate = model_summary$tTable[, "Value"],
    std.error = model_summary$tTable[, "Std.Error"],
    statistic = model_summary$tTable[, "t-value"],
    p.value = model_summary$tTable[, "p-value"]
  ) |>
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      significance = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        p.value < 0.1 ~ ".",
        TRUE ~ ""
      ),
      effect_type = case_when(
        term == "(Intercept)" ~ "Intercept",
        p.value < 0.05 & estimate < 0 ~ "Negative (sig.)",
        p.value < 0.05 & estimate > 0 ~ "Positive (sig.)",
        TRUE ~ "Non-significant"
      )
    )
  
  # Clean term names based on land cover type
  if (land_cover_type == "ftws") {
    coef_df <- coef_df |>
      mutate(term_clean = case_when(
        term == "(Intercept)" ~ "Intercept",
        term == "forest_to_tws_prop" ~ "Forest → TWS",
        term == "forest_no_change_prop" ~ "Forest (no change)",
        term == "delta_recorder_effort" ~ "ΔRecorder effort",
        term == "log_recorder_effort" ~ "log(Recorder effort)",
        term == "lc_time_period2006-2012" ~ "Period: 2006-2012",
        term == "lc_time_period2012-2018" ~ "Period: 2012-2018",
        term == "temp_change" ~ "Temperature change",
        term == "precip_change" ~ "Precipitation change",
        TRUE ~ term
      ))
    level_order <- rev(c("Forest → TWS", "Forest (no change)", "Temperature change",
                         "Precipitation change", "ΔRecorder effort", "log(Recorder effort)",
                         "Period: 2006-2012", "Period: 2012-2018"))
  } else if (land_cover_type == "twsf") {
    coef_df <- coef_df |>
      mutate(term_clean = case_when(
        term == "(Intercept)" ~ "Intercept",
        term == "tws_to_forest_prop" ~ "TWS → Forest",
        term == "tws_no_change_prop" ~ "TWS (no change)",
        term == "delta_recorder_effort" ~ "ΔRecorder effort",
        term == "log_recorder_effort" ~ "log(Recorder effort)",
        term == "lc_time_period2006-2012" ~ "Period: 2006-2012",
        term == "lc_time_period2012-2018" ~ "Period: 2012-2018",
        term == "temp_change" ~ "Temperature change",
        term == "precip_change" ~ "Precipitation change",
        TRUE ~ term
      ))
    level_order <- rev(c("TWS → Forest", "TWS (no change)", "Temperature change",
                         "Precipitation change", "ΔRecorder effort", "log(Recorder effort)",
                         "Period: 2006-2012", "Period: 2012-2018"))
  } else {  # urban
    coef_df <- coef_df |>
      mutate(term_clean = case_when(
        term == "(Intercept)" ~ "Intercept",
        term == "all_to_urban_prop" ~ "All → Urban",
        term == "urban_no_change_prop" ~ "Urban (no change)",
        term == "delta_recorder_effort" ~ "ΔRecorder effort",
        term == "log_recorder_effort" ~ "log(Recorder effort)",
        term == "lc_time_period2006-2012" ~ "Period: 2006-2012",
        term == "lc_time_period2012-2018" ~ "Period: 2012-2018",
        term == "temp_change" ~ "Temperature change",
        term == "precip_change" ~ "Precipitation change",
        TRUE ~ term
      ))
    level_order <- rev(c("All → Urban", "Urban (no change)", "Temperature change",
                         "Precipitation change", "ΔRecorder effort", "log(Recorder effort)",
                         "Period: 2006-2012", "Period: 2012-2018"))
  }
  
  # Split into intercept and effects
  intercept_df <- coef_df |> filter(term == "(Intercept)")
  effects_df <- coef_df |>
    filter(term != "(Intercept)") |>
    mutate(term_clean = factor(term_clean, levels = level_order))
  
  # Define colors - ensure proper factor levels
  effect_colors <- c(
    "Positive (sig.)" = "#FF9800",
    "Negative (sig.)" = "#9C27B0",
    "Non-significant" = "grey60",
    "Intercept" = "grey30"
  )
  
  # Create intercept plot (top)
  plot_intercept <- ggplot(intercept_df, aes(x = estimate, y = term_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.2, linewidth = 0.7, color = "grey30") +
    geom_point(aes(color = effect_type), size = 3, shape = 15) +
    geom_text(aes(x = conf.high, label = significance),
              hjust = -0.3, vjust = 0.5, size = 3.5, fontface = "bold") +
    scale_color_manual(values = effect_colors, guide = "none") +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 10, color = "black", face = "bold"),
      axis.text.x = element_text(size = 9, color = "black"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      plot.margin = margin(t = 5, r = 10, b = 2, l = 5)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)))
  
  # Create effects plot (bottom)
  plot_effects <- ggplot(effects_df, aes(x = estimate, y = term_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.2, linewidth = 0.7, color = "grey30") +
    geom_point(aes(color = effect_type), size = 3) +
    geom_text(aes(x = conf.high, label = significance),
              hjust = -0.3, vjust = 0.5, size = 3.5, fontface = "bold") +
    scale_color_manual(
      values = effect_colors,
      name = "Effect type",
      breaks = c("Positive (sig.)", "Negative (sig.)", "Non-significant"),
      labels = c("Positive (p < 0.05)", "Negative (p < 0.05)", "Non-significant")
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 9, color = "black"),
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      plot.margin = margin(t = 2, r = 10, b = 5, l = 5)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)))
  
  # Show or hide y-axis labels
  if (!show_y_axis) {
    plot_intercept <- plot_intercept +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    plot_effects <- plot_effects +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  # Combine intercept and effects for this model
  combined <- plot_intercept / plot_effects +
    plot_layout(heights = c(1, 4))
  
  return(combined)
}

# Create individual plots for each model
# Row 1: Forest → TWS
ftws_all_plot <- create_coef_plot_beta_jac(all_FTWS_beta_JAC_model1, "ftws", show_y_axis = TRUE)
ftws_plants_plot <- create_coef_plot_beta_jac(plants_FTWS_beta_JAC_model1, "ftws", show_y_axis = FALSE)
ftws_birds_plot <- create_coef_plot_beta_jac(birds_FTWS_beta_JAC_model1, "ftws", show_y_axis = FALSE)

# Row 2: TWS → Forest
twsf_all_plot <- create_coef_plot_beta_jac(all_TWSF_beta_JAC_model1, "twsf", show_y_axis = TRUE)
twsf_plants_plot <- create_coef_plot_beta_jac(plants_TWSF_beta_JAC_model1, "twsf", show_y_axis = FALSE)
twsf_birds_plot <- create_coef_plot_beta_jac(birds_TWSF_beta_JAC_model1, "twsf", show_y_axis = FALSE)

# Row 3: All → Urban
urban_all_plot <- create_coef_plot_beta_jac(all_urban_beta_JAC_model1, "urban", show_y_axis = TRUE)
urban_plants_plot <- create_coef_plot_beta_jac(plants_urban_beta_JAC_model1, "urban", show_y_axis = FALSE)
urban_birds_plot <- create_coef_plot_beta_jac(birds_urban_beta_JAC_model1, "urban", show_y_axis = FALSE)

# Extract legend from the effects plot (which has the legend configured)
legend <- get_legend(
  ggplot(data.frame(
    x = c(1, 2, 3),
    y = c(1, 2, 3),
    effect_type = factor(c("Positive (sig.)", "Negative (sig.)", "Non-significant"),
                         levels = c("Positive (sig.)", "Negative (sig.)", "Non-significant"))
  ), aes(x = x, y = y, color = effect_type)) +
    geom_point(size = 3) +
    scale_color_manual(
      values = c(
        "Positive (sig.)" = "#FF9800",
        "Negative (sig.)" = "#9C27B0",
        "Non-significant" = "grey60"
      ),
      name = "Effect type",
      breaks = c("Positive (sig.)", "Negative (sig.)", "Non-significant"),
      labels = c("Positive (p < 0.05)", "Negative (p < 0.05)", "Non-significant")
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    )
)

# Remove legends from individual plots
ftws_all_plot <- ftws_all_plot & theme(legend.position = "none")
ftws_plants_plot <- ftws_plants_plot & theme(legend.position = "none")
ftws_birds_plot <- ftws_birds_plot & theme(legend.position = "none")
twsf_all_plot <- twsf_all_plot & theme(legend.position = "none")
twsf_plants_plot <- twsf_plants_plot & theme(legend.position = "none")
twsf_birds_plot <- twsf_birds_plot & theme(legend.position = "none")
urban_all_plot <- urban_all_plot & theme(legend.position = "none")
urban_plants_plot <- urban_plants_plot & theme(legend.position = "none")
urban_birds_plot <- urban_birds_plot & theme(legend.position = "none")

# Create column headers using wrap_elements with grid
col1_header <- wrap_elements(grid::textGrob("All Occurrences", 
                                            gp = grid::gpar(fontsize = 12, fontface = "bold")))
col2_header <- wrap_elements(grid::textGrob("Vascular Plants", 
                                            gp = grid::gpar(fontsize = 12, fontface = "bold")))
col3_header <- wrap_elements(grid::textGrob("Birds", 
                                            gp = grid::gpar(fontsize = 12, fontface = "bold")))

# Combine all plots into a 3x3 grid with headers
combined_plot <- (
  (col1_header | col2_header | col3_header) /
    (ftws_all_plot | ftws_plants_plot | ftws_birds_plot) /
    (twsf_all_plot | twsf_plants_plot | twsf_birds_plot) /
    (urban_all_plot | urban_plants_plot | urban_birds_plot) /
    legend
) +
  plot_layout(heights = c(0.05, 1, 1, 1, 0.08)) +
  plot_annotation(
    title = "Beta_jac (Jaccard Dissimilarity) - Model Coefficients Across Land Cover Changes",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

# Save the plot
ggsave(
  filename = here("figures", "beta_jac_combined_coefficients.png"),
  plot = combined_plot,
  width = 18,
  height = 16,
  dpi = 600,
  bg = "white"
)

# END OF SCRIPT ----------------------------------------------------------------