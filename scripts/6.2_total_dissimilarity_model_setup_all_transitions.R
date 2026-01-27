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
all_FTWS_beta_JAC_model1 <- gls(beta_jac ~ forest_to_tws + forest_no_change + delta_recorder_effort + 
                                  log_recorder_effort + lc_time_period + temp_change + precip_change,
                                correlation = corExp(form = ~ x + y | time_numeric),
                                data = turnover_forest_tws_15km_all,
                                method = "REML")

# Save model
save(all_FTWS_beta_JAC_model1,
     file = here("data", "models", "final", "all_FTWS_beta_jac_model1.RData"))

# Fit GLS with interaction term
all_FTWS_beta_JAC_model2_interaction <-  gls(beta_jac ~ forest_to_tws * temp_change +
                                               forest_to_tws * precip_change + 
                                               forest_no_change +
                                               delta_recorder_effort + 
                                               log_recorder_effort + 
                                               lc_time_period,
                                             correlation = corExp(form = ~ x + y | time_numeric),
                                             data = turnover_forest_tws_15km_coords_time,
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
         log_recorder_effort = log(recorder_effort))

# Fit GLS without interaction
plants_FTWS_beta_JAC_model1 <- gls(beta_jac ~ forest_to_tws + forest_no_change +
                                     delta_recorder_effort + log_recorder_effort 
                                   + lc_time_period + temp_change + precip_change,
                                   correlation = corExp(form = ~ x + y | time_numeric),
                                   data = turnover_forest_tws_15km_plants,
                                   method = "REML")

# Save model
save(plants_FTWS_beta_JAC_model1,
     file = here("data", "models", "final", "plants_FTWS_beta_JAC_model1.RData"))

# Fit GLS with interaction term
plants_FTWS_beta_JAC_model2_interaction <-  gls(beta_jac ~ forest_to_tws * temp_change +
                                                  forest_to_tws * precip_change + 
                                                  forest_no_change +
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
         log_recorder_effort = log(recorder_effort))

# Fit GLS without interaction
birds_FTWS_beta_JAC_model1 <- gls(beta_jac ~ forest_to_tws + forest_no_change + 
                                    delta_recorder_effort + log_recorder_effort +
                                    lc_time_period + temp_change + precip_change,
                                  correlation = corExp(form = ~ x + y | time_numeric),
                                  data = turnover_forest_tws_15km_birds,
                                  method = "REML")

# Save model
save(birds_FTWS_beta_JAC_model1,
     file = here("data", "models", "final", "birds_FTWS_beta_JAC_model1.RData"))


# Fit GLS with interaction
birds_FTWS_beta_JAC_model2_interaction <-  gls(beta_jac ~ forest_to_tws * temp_change +
                                                 forest_to_tws * precip_change + 
                                                 forest_no_change +
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
         log_recorder_effort = log(recorder_effort))

# Fit GLS without interaction 
all_TWSF_beta_JAC_model1 <- gls(beta_jac ~ tws_to_forest + tws_no_change + 
                                  delta_recorder_effort + log_recorder_effort + 
                                  lc_time_period + temp_change + precip_change,
                                correlation = corExp(form = ~ x + y | time_numeric),
                                data = turnover_tws_forest_15km_all,
                                method = "REML")

# Save model
save(all_TWSF_beta_JAC_model1,
     file = here("data", "models", "final", "all_TWSF_beta_JAC_model1.RData"))

# Fit GLS with interaction term
all_TWSF_beta_JAC_model2_interaction <-  gls(beta_jac ~ tws_to_forest * temp_change +
                                               tws_to_forest * precip_change + 
                                               tws_no_change +
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
         log_recorder_effort = log(recorder_effort))


# Fit GLS without interaction
plants_TWSF_beta_JAC_model1 <- gls(beta_jac ~ tws_to_forest + tws_no_change + 
                                     delta_recorder_effort + log_recorder_effort +
                                     lc_time_period + temp_change + precip_change,
                                   correlation = corExp(form = ~ x + y | time_numeric),
                                   data = turnover_tws_forest_15km_plants,
                                   method = "REML")

# Save model output
save(plants_TWSF_beta_JAC_model1,
     file = here("data", "models", "final", "plants_TWSF_beta_JAC_model1.RData"))

# Fit GLS with interaction term
plants_TWSF_beta_JAC_model2_interaction <-  gls(beta_jac ~ tws_to_forest * temp_change +
                                                  tws_to_forest * precip_change + 
                                                  tws_no_change +
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
         log_recorder_effort = log(recorder_effort))


# Fit GLS without interaction
birds_TWSF_beta_JAC_model1 <- gls(beta_jac ~ tws_to_forest + tws_no_change + 
                                    delta_recorder_effort + log_recorder_effort +
                                    lc_time_period + temp_change + precip_change,
                                  correlation = corExp(form = ~ x + y | time_numeric),
                                  data = turnover_tws_forest_15km_birds,
                                  method = "REML")

# Save model output
save(birds_TWSF_beta_JAC_model1,
     file = here("data", "models", "final", "birds_TWSF_beta_JAC_model1.RData"))

# Fit GLS with interaction
birds_TWSF_beta_JAC_model2_interaction <-  gls(beta_jac ~ tws_to_forest * temp_change +
                                                 tws_to_forest * precip_change + 
                                                 tws_no_change +
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
         log_recorder_effort = log(recorder_effort))

# Fit GLS without interaction
all_urban_beta_JAC_model1 <- gls(beta_jac ~ all_to_urban + urban_no_change +
                                   delta_recorder_effort + log_recorder_effort + 
                                   lc_time_period + temp_change + precip_change,
                                 correlation = corExp(form = ~ x + y | time_numeric),
                                 data = turnover_all_urban_15km_all,
                                 method = "REML")

# Save model
save(all_urban_beta_JAC_model1,
     file = here("data", "models", "final", "all_urban_beta_JAC_model1.RData"))

# Fit GLS with interaction 
all_urban_beta_JAC_model2_interaction <-  gls(beta_jac ~ all_to_urban * temp_change +
                                                all_to_urban * precip_change + 
                                                urban_no_change +
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
         log_recorder_effort = log(recorder_effort))

# Fit GLS without interaction
plants_urban_beta_JAC_model1 <- gls(beta_jac ~ all_to_urban + urban_no_change +
                                      delta_recorder_effort + log_recorder_effort + 
                                      lc_time_period + temp_change + precip_change,
                                    correlation = corExp(form = ~ x + y | time_numeric),
                                    data = turnover_all_urban_15km_plants,
                                    method = "REML")

# Save model 
save(plants_urban_beta_JAC_model1,
     file = here("data", "models", "final", "plants_urban_beta_JAC_model1.RData"))

# Fit GLS with interaction term
plants_urban_beta_JAC_model2_interaction <-  gls(beta_jac ~ all_to_urban * temp_change +
                                                   all_to_urban * precip_change + 
                                                   urban_no_change +
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
         log_recorder_effort = log(recorder_effort))

# Fit GLS model without interaction
birds_urban_beta_JAC_model1 <- gls(beta_jac ~ all_to_urban + urban_no_change + 
                                     delta_recorder_effort + log_recorder_effort + 
                                     lc_time_period + temp_change + precip_change,
                                   correlation = corExp(form = ~ x + y | time_numeric),
                                   data = turnover_all_urban_15km_birds,
                                   method = "REML")

# Save model
save(birds_urban_beta_JAC_model1,
     file = here("data", "models", "final", "birds_urban_beta_JAC_model1.RData"))

# Fit GLS with interaction
birds_urban_beta_JAC_model2_interaction <-  gls(beta_jac ~ all_to_urban * temp_change +
                                                  all_to_urban * precip_change + 
                                                  urban_no_change +
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

# END OF SCRIPT ----------------------------------------------------------------