##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 5.1_model_setup_turnover_forest_tws
# This script contains code which sets up the models exploring the impact of 
# Forest -> TWS land cover transition on temporal turnover (beta_jtu) for all 
# occurrences, plant-only occurrences and bird-only occurrences
##----------------------------------------------------------------------------##

library(here)
source(here("scripts", "0_setup.R"))

# 1. LOAD DATA -----------------------------------------------------------------

# Turnover data for plants
load(here("data", "derived_data", 
          "vascular_plants_turnover_all_land_cover_climate_15km.rda"))

# Turnover data for birds
load(here("data", "derived_data", 
          "bird_turnover_all_land_cover_climate_15km.rda"))

# 2. FOREST -> TWS TURNOVER (BETA_JTU) -----------------------------------------

# Check column names
colnames(vascular_plants_turnover_with_climate)
colnames(birds_turnover_with_climate)

# Total number of 100m x 100m pixels in a 15km x 15km cell
# 15,000m / 100m = 150 pixels per side
# 150 x 150 = 22,500 total pixels per cell
total_pixels_per_cell <- 22500

## 2.1. Vascular Plants --------------------------------------------------------

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
  # prepare data for GLS (i.e. remove rows with missing x or y and categorise time periods)
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort),
         forest_no_change_prop = forest_no_change / total_pixels_per_cell,
         forest_to_tws_prop = forest_to_tws / total_pixels_per_cell)

# Check transformation went ok
head(turnover_forest_tws_15km_plants)

# Check to make sure no 0 got logged
any(!is.finite(turnover_forest_tws_15km_plants$log_recorder_effort)) # FALSE!

# Check if there are any cells with recorder effort = 0
a <- turnover_forest_tws_15km_plants |>
  filter(recorder_effort == 0)
a #0 - Good!

# Fit GLS with logged recorder effort - no interaction
plants_FTWS_turnover_model1 <- gls(beta_jtu ~ forest_to_tws_prop + forest_no_change_prop +
                                     delta_recorder_effort + log_recorder_effort +
                                     lc_time_period + temp_change + precip_change,
                                   correlation = corExp(form = ~ x + y | time_numeric),
                                   data = turnover_forest_tws_15km_plants,
                                   method = "REML")

# Fit GLS with logged recorder effort and interaction
plants_FTWS_turnover_model2_interaction <- gls(beta_jtu ~ forest_to_tws_prop * temp_change +
                                                 forest_to_tws_prop * precip_change +
                                                 forest_no_change_prop +
                                                 delta_recorder_effort +
                                                 log_recorder_effort +
                                                 lc_time_period,
                                               correlation = corExp(form = ~ x + y | time_numeric),
                                               data = turnover_forest_tws_15km_plants,
                                               method = "REML")

# Compare models based on AIC
AICtab(plants_FTWS_turnover_model1, plants_FTWS_turnover_model2_interaction, base = TRUE)
# AIC     dAIC    df
# plants_FTWS_turnover_model1             -1939.7     0.0 11
# plants_FTWS_turnover_model2_interaction -1928.6    11.1 13

# Save models (higher AIC will go into exploratory, lower AIC will go into final)
save(plants_FTWS_turnover_model1,
     file = here("data", "models", "final", "plants_FTWS_turnover_model1.RData"))
save(plants_FTWS_turnover_model2_interaction,
     file = here("data", "models", "exploratory", "plants_FTWS_turnover_model2_interaction.RData"))

# Extract residuals
plant_residuals_gls <- residuals(plants_FTWS_turnover_model1, type = "normalized")

# Basic residual plots
par(mfrow = c(2, 2))

# Residuals vs fitted
plot(fitted(plants_FTWS_turnover_model1), plant_residuals_gls,
     xlab = "Fitted Values", ylab = "Normalized Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# QQ plot
qqnorm(plant_residuals_gls, main = "Normal Q-Q Plot")
qqline(plant_residuals_gls, col = "red")

# Residuals vs predictors
plot(turnover_forest_tws_15km_plants$forest_to_tws_prop, plant_residuals_gls,
     xlab = "Forest to TWS", ylab = "Normalized Residuals",
     main = "Residuals vs Forest to TWS")
abline(h = 0, col = "red", lty = 2)

plot(turnover_forest_tws_15km_plants$log_recorder_effort, plant_residuals_gls,
     xlab = "Recorder Effort", ylab = "Normalized Residuals",
     main = "Residuals vs Recorder Effort")
abline(h = 0, col = "red", lty = 2)

## 2.2. Birds ------------------------------------------------------------------

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

# Check to make sure no 0 got logged
any(!is.finite(turnover_forest_tws_15km_birds$log_recorder_effort)) # FALSE!

# Check if there are any cells with recorder effort = 0
b <- turnover_forest_tws_15km_birds |>
  filter(recorder_effort == 0)
b #0 - Good!

# Fit GLS with logged recorder effort - no interaction
birds_FTWS_turnover_model1 <- gls(beta_jtu ~ forest_to_tws_prop + forest_no_change_prop +
                                    delta_recorder_effort + log_recorder_effort +
                                    lc_time_period + temp_change + precip_change,
                                  correlation = corExp(form = ~ x + y | time_numeric),
                                  data = turnover_forest_tws_15km_birds,
                                  method = "REML")

# Fit GLS with interaction term
birds_FTWS_turnover_model2_interaction <- gls(beta_jtu ~ forest_to_tws_prop * temp_change +
                                                forest_to_tws_prop * precip_change +
                                                forest_no_change_prop +
                                                delta_recorder_effort +
                                                log_recorder_effort +
                                                lc_time_period,
                                              correlation = corExp(form = ~ x + y | time_numeric),
                                              data = turnover_forest_tws_15km_birds,
                                              method = "REML")

# Compare models based on AIC
AICtab(birds_FTWS_turnover_model1, birds_FTWS_turnover_model2_interaction, base = TRUE)

# Save models
save(birds_FTWS_turnover_model1,
     file = here("data", "models", "final", "birds_FTWS_turnover_model1.RData"))
save(birds_FTWS_turnover_model2_interaction,
     file = here("data", "models", "exploratory", "birds_FTWS_turnover_model2_interaction.RData"))

# Extract residuals
bird_residuals_gls <- residuals(birds_FTWS_turnover_model1, type = "normalized")

par(mfrow = c(2, 2))

# Residuals vs fitted
plot(fitted(birds_FTWS_turnover_model1), bird_residuals_gls,
     xlab = "Fitted Values", ylab = "Normalized Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# QQ plot
qqnorm(bird_residuals_gls, main = "Normal Q-Q Plot")
qqline(bird_residuals_gls, col = "red")

# Residuals vs land-cover change
plot(turnover_forest_tws_15km_birds$forest_to_tws_prop, bird_residuals_gls,
     xlab = "Forest to TWS", ylab = "Normalized Residuals",
     main = "Residuals vs Forest to TWS")
abline(h = 0, col = "red", lty = 2)

# Residuals vs recorder effort
plot(turnover_forest_tws_15km_birds$log_recorder_effort, bird_residuals_gls,
     xlab = "Recorder Effort", ylab = "Normalized Residuals",
     main = "Residuals vs Recorder Effort")
abline(h = 0, col = "red", lty = 2)

# 3. TWS -> FOREST -------------------------------------------------------------

## 3.1. Vascular Plants --------------------------------------------------------

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

# Check to make sure no 0 got logged
any(!is.finite(turnover_tws_forest_15km_plants$log_recorder_effort)) # FALSE!

# Check if there are any cells with recorder effort = 0
c <- turnover_tws_forest_15km_plants |>
  filter(recorder_effort == 0)
c #0 - Good!

# Fit GLS with logged recorder effort and no interaction
plants_TWSF_turnover_model1 <- gls(beta_jtu ~ tws_to_forest_prop + tws_no_change_prop +
                                     delta_recorder_effort + log_recorder_effort +
                                     lc_time_period + temp_change + precip_change,
                                   correlation = corExp(form = ~ x + y | time_numeric),
                                   data = turnover_tws_forest_15km_plants,
                                   method = "REML")

# Fit GLS with logged recorder effort and interaction
plants_TWSF_turnover_model2_interaction <- gls(beta_jtu ~ tws_to_forest_prop * temp_change +
                                                 tws_to_forest_prop * precip_change +
                                                 tws_no_change_prop +
                                                 delta_recorder_effort +
                                                 log_recorder_effort +
                                                 lc_time_period,
                                               correlation = corExp(form = ~ x + y | time_numeric),
                                               data = turnover_tws_forest_15km_plants,
                                               method = "REML")

# Compare models
AICtab(plants_TWSF_turnover_model1, plants_TWSF_turnover_model2_interaction, base = TRUE)

# Save models
save(plants_TWSF_turnover_model1,
     file = here("data", "models", "final", "plants_TWSF_turnover_model1.RData"))
save(plants_TWSF_turnover_model2_interaction,
     file = here("data", "models", "exploratory", "plants_TWSF_turnover_model2_interaction.RData"))

# Extract model residuals
plant_residuals_gls <- residuals(plants_TWSF_turnover_model1, type = "normalized")

# Residuals vs fitted
plot(fitted(plants_TWSF_turnover_model1), plant_residuals_gls,
     xlab = "Fitted Values", ylab = "Normalized Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# QQ plot
qqnorm(plant_residuals_gls, main = "Normal Q-Q Plot")
qqline(plant_residuals_gls, col = "red")

# Residuals vs land-cover change
plot(turnover_tws_forest_15km_plants$tws_to_forest_prop, plant_residuals_gls,
     xlab = "TWS to Forest", ylab = "Normalized Residuals",
     main = "Residuals vs TWS to Forest")
abline(h = 0, col = "red", lty = 2)

# Residuals vs recorder effort
plot(turnover_tws_forest_15km_plants$log_recorder_effort, plant_residuals_gls,
     xlab = "Recorder Effort", ylab = "Normalized Residuals",
     main = "Residuals vs Recorder Effort")
abline(h = 0, col = "red", lty = 2)


# 3. PLANT OCCURRENCES ONLY ----------------------------------------------------


## 3.1. Prepare plant data for analysis ----------------------------------------

# Select only Forest -> TWS columns
# check column names
colnames(vascular_plants_turnover_with_climate)

# Also rename the columns for easier manipulation of df
plants_turnover_forest_tws_15km <- vascular_plants_turnover_with_climate |>
  select(-c('2000-2006_TWS no change', '2000-2006_TWS to Forest',
            '2000-2006_Urban_no_change', '2000-2006_all_to_urban',
            '2006-2012_TWS no change', '2006-2012_TWS to Forest',
            '2006-2012_Urban_no_change', '2006-2012_all_to_urban',
            '2012-2018_TWS no change', '2012-2018_TWS to Forest',
            '2012-2018_Urban_no_change', '2012-2018_all_to_urban')) |>
  # determine which rows belong to which time period
  mutate(forest_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest no change`,
                                      lc_time_period == "2006-2012" ~ `2006-2012_Forest no change`,
                                      lc_time_period == "2012-2018" ~ `2012-2018_Forest no change`,
                                      TRUE ~ NA_real_),
         forest_to_tws = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest to TWS`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_Forest to TWS`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_Forest to TWS`,
                                   TRUE ~ NA_real_)) |>
  # convert land-cover changes to proportions
  mutate(forest_no_change_prop = forest_no_change / total_pixels_per_cell,
         forest_to_tws_prop = forest_to_tws / total_pixels_per_cell) |>
  # remove columns no longer required
  select(-`2000-2006_Forest no change`, -`2006-2012_Forest no change`,
         -`2012-2018_Forest no change`,-`2000-2006_Forest to TWS`,
         -`2006-2012_Forest to TWS`, -`2012-2018_Forest to TWS`)

# Check transformation went ok
head(plants_turnover_forest_tws_15km)

# Remove rows that might have NA for x or y and categorise time periods for GLS
plants_turnover_forest_tws_15km_coords_time <- plants_turnover_forest_tws_15km |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3))

## 3.2. Run GLS on the plant data ----------------------------------------------

# Check if there are any cells with recorder effort = 0
b <- plants_turnover_forest_tws_15km_coords_time |>
  filter(recorder_effort == 0)
length(a) #0 - Good!

# Log transform recorder effort values
plants_turnover_forest_tws_15km_coords_time <- plants_turnover_forest_tws_15km_coords_time |>
  mutate(log_recorder_effort = log(recorder_effort))

# Check log transformed values
summary(plants_turnover_forest_tws_15km_coords_time$recorder_effort)
any(!is.finite(plants_turnover_forest_tws_15km_coords_time$recorder_effort)) # FALSE = no infinite values - Good!

# Define GLS
plants_FTWS_model1_gls <- gls(beta_jtu ~ forest_to_tws_prop + forest_no_change_prop +
                    delta_recorder_effort + log_recorder_effort + lc_time_period + temp_change + precip_change,
                  correlation = corExp(form = ~ x + y | time_numeric),
                  data = plants_turnover_forest_tws_15km_coords_time,
                  method = "REML")

# Model validation revealed that this model is ok to use
  # despite some deviances, we will use this model to keep things consistent

# Save model output to file
save(plants_FTWS_model1_gls,
     file = here("data", "models", "final",
                 "plants_FTWS_model1_gls.RData"))

## 3.3. Plant GLS with interaction ---------------------------------------------

# Define GLS
plants_FTWS_model2_gls_interaction <- gls(beta_jtu ~ forest_to_tws_prop * temp_change +
                                            forest_to_tws_prop * precip_change + 
                                            forest_no_change_prop +
                                                  delta_recorder_effort + 
                                                  log_recorder_effort + 
                                                  lc_time_period,
                                                correlation = corExp(form = ~ x + y | time_numeric),
                                                data = plants_turnover_forest_tws_15km_coords_time,
                                                method = "REML")
# Model validaiton looked ok enough

# Compare models based on AIC
AICtab(plants_FTWS_model1_gls, plants_FTWS_model2_gls_interaction, base = TRUE)
# model without interaction preferred => save interaction model in exploratory folder

# Save model
save(plants_FTWS_model2_gls_interaction,
     file = here("data", "models", "exploratory",
                 "plants_FTWS_model2_gls_interactionn.RData"))

## 3.4. Get model summary ------------------------------------------------------

# Get summary of final model
model_summary <- summary(plants_FTWS_model1_gls)

# Create dataframe of coefficients
coef_df <- data.frame(term = names(model_summary$tTable[, "Value"]),
                      estimate = model_summary$tTable[, "Value"],
                      std.error = model_summary$tTable[, "Std.Error"],
                      statistic = model_summary$tTable[, "t-value"],
                      p.value = model_summary$tTable[, "p-value"])

# Add significance stars
coef_df <- coef_df |>
  mutate(significance = case_when( p.value < 0.001 ~ "***",
                                   p.value < 0.01 ~ "**",
                                   p.value < 0.05 ~ "*",
                                   p.value < 0.1 ~ ".",
                                   TRUE ~ ""))

coef_table <- coef_df |>
  select(term, estimate, std.error, statistic, p.value, significance) |>
  flextable() |>
  set_header_labels(term = "Predictor", estimate = "Estimate", std.error = "SE",
                    statistic = "t-value", p.value = "p-value", significance = "") |>
  colformat_double(j = c("estimate", "std.error", "statistic"), digits = 3) |>
  colformat_double(j = "p.value", digits = 4) |>
  autofit()

# Display table
coef_table

# 4. BIRD OCCURRENCES ONLY -----------------------------------------------------

## 4.1. Prepare bird data for analysis -----------------------------------------

# Check column names
colnames(birds_turnover_with_climate)

# Also rename the columns for easier manipulation of df
birds_turnover_forest_tws_15km <- birds_turnover_with_climate |>
  select(-c('2000-2006_TWS no change', '2000-2006_TWS to Forest',
            '2000-2006_Urban_no_change', '2000-2006_all_to_urban',
            '2006-2012_TWS no change', '2006-2012_TWS to Forest',
            '2006-2012_Urban_no_change', '2006-2012_all_to_urban',
            '2012-2018_TWS no change', '2012-2018_TWS to Forest',
            '2012-2018_Urban_no_change', '2012-2018_all_to_urban')) |>
  # determine which rows belong to which time period
  mutate(forest_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest no change`,
                                      lc_time_period == "2006-2012" ~ `2006-2012_Forest no change`,
                                      lc_time_period == "2012-2018" ~ `2012-2018_Forest no change`,
                                      TRUE ~ NA_real_),
         forest_to_tws = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest to TWS`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_Forest to TWS`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_Forest to TWS`,
                                   TRUE ~ NA_real_)) |>
  # convert land-cover changes to proportions
  mutate(forest_no_change_prop = forest_no_change / total_pixels_per_cell,
         forest_to_tws_prop = forest_to_tws / total_pixels_per_cell) |>
  # remove columns no longer required
  select(-`2000-2006_Forest no change`, -`2006-2012_Forest no change`, 
         -`2012-2018_Forest no change`,-`2000-2006_Forest to TWS`, 
         -`2006-2012_Forest to TWS`, -`2012-2018_Forest to TWS`)

# Check transformation went ok
head(birds_turnover_forest_tws_15km)

# Remove rows that might have NA for x or y & categorise time periods
birds_turnover_forest_tws_15km_coords_time <- birds_turnover_forest_tws_15km |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3))

## 4.2. Run GLS on bird data ---------------------------------------------------

# Check if there are any cells with recorder effort = 0
c <- birds_turnover_forest_tws_15km_coords_time |>
  filter(recorder_effort == 0)
length(c) #0 - Good!

# Log transform recorder effort values
birds_turnover_forest_tws_15km_coords_time <- birds_turnover_forest_tws_15km_coords_time |>
  mutate(log_recorder_effort = log(recorder_effort))

# Check log transformed values
summary(birds_turnover_forest_tws_15km_coords_time$recorder_effort)
any(!is.finite(birds_turnover_forest_tws_15km_coords_time$recorder_effort)) # FALSE = no infinite values - Good!

# Define GLS
birds_FTWS_model1_gls <- gls(beta_jtu ~ forest_to_tws_prop + forest_no_change_prop +
                           delta_recorder_effort + log_recorder_effort + lc_time_period + temp_change + precip_change,
                         correlation = corExp(form = ~ x + y | time_numeric),
                         data = birds_turnover_forest_tws_15km_coords_time,
                         method = "REML")

# Save model output to file
save(birds_FTWS_model1_gls,
     file = here("data", "models", "final",
                 "birds_FTWS_model1_gls.RData"))

## 4.3. Bird GLS with interaction ----------------------------------------------

# Define GLS
birds_FTWS_model2_gls_interaction <- gls(beta_jtu ~ forest_to_tws_prop * temp_change +
                                           forest_to_tws_prop * precip_change + 
                                           forest_no_change_prop +
                                            delta_recorder_effort + 
                                            log_recorder_effort + 
                                            lc_time_period,
                                          correlation = corExp(form = ~ x + y | time_numeric),
                                          data = birds_turnover_forest_tws_15km_coords_time,
                                          method = "REML")
# Model validation looked ok enough

# Compare models based on AIC
AICtab(birds_FTWS_model1_gls, birds_FTWS_model2_gls_interaction, base = TRUE)
# model without interaction preferred => save interaction model in exploratory folder

# Save model
save(birds_FTWS_model2_gls_interaction,
     file = here("data", "models", "exploratory",
                 "birds_FTWS_model2_gls_interaction.RData"))

## 4.4. Get model summary ------------------------------------------------------

# Get summary of final model
model_summary <- summary(birds_FTWS_model1_gls)

# Create dataframe of coefficients
coef_df <- data.frame(term = names(model_summary$tTable[, "Value"]),
                      estimate = model_summary$tTable[, "Value"],
                      std.error = model_summary$tTable[, "Std.Error"],
                      statistic = model_summary$tTable[, "t-value"],
                      p.value = model_summary$tTable[, "p-value"])

# Add significance stars
coef_df <- coef_df |>
  mutate(significance = case_when( p.value < 0.001 ~ "***",
                                   p.value < 0.01 ~ "**",
                                   p.value < 0.05 ~ "*",
                                   p.value < 0.1 ~ ".",
                                   TRUE ~ ""))

coef_table <- coef_df |>
  select(term, estimate, std.error, statistic, p.value, significance) |>
  flextable() |>
  set_header_labels(term = "Predictor", estimate = "Estimate", std.error = "SE",
                    statistic = "t-value", p.value = "p-value", significance = "") |>
  colformat_double(j = c("estimate", "std.error", "statistic"), digits = 3) |>
  colformat_double(j = "p.value", digits = 4) |>
  autofit()

# Display table
coef_table

# 5. PLOT MODEL OUTPUTS --------------------------------------------------------

# Use function defined in section 14.1 of the 0_setup.R script

# Create plots for each model
plot_all <- create_coef_plot_forest(FTWS_turnover_model5_gls_log, 
                             "a) All occurrences", 
                             show_y_axis = TRUE)

plot_plants <- create_coef_plot_forest(plants_FTWS_model1_gls, 
                                "b) Vascular plants", 
                                show_y_axis = FALSE)

plot_birds <- create_coef_plot_forest(birds_FTWS_model1_gls, 
                               "c) Birds", 
                               show_y_axis = FALSE)

# Combine all three panels horizontally
combined_plot <- plot_all | plot_plants | plot_birds

# Create a dummy plot just to extract the legend
dummy_effects <- data.frame(term_clean = "dummy",estimate = 0,
                            conf.low = -0.1, conf.high = 0.1,
                            effect_type = c("Positive (sig.)", "Negative (sig.)", "Non-significant"))

effect_colors <- c("Negative (sig.)" = "#9C27B0",
                   "Positive (sig.)" = "#FF9800",
                   "Non-significant" = "grey60")

legend_plot <- ggplot(dummy_effects, aes(x = estimate, y = term_clean, color = effect_type)) +
  geom_point() +
  scale_color_manual(
    values = effect_colors,
    name = NULL,
    breaks = c("Positive (sig.)", "Negative (sig.)", "Non-significant"),
    labels = c("Positive effect (p < 0.05)", "Negative effect (p < 0.05)", "Non-significant")
  ) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 9))

# Extract legend
legend <- get_legend(legend_plot)

# Final combined plot with legend
final_plot <- combined_plot / legend +
  plot_layout(heights = c(20, 1)) +
  plot_annotation(theme = theme(plot.title = element_text(size = 13, 
                                                          face = "bold", 
                                                          hjust = 0.5, 
                                                          margin = margin(b = 10))))

# Display
print(final_plot)

# Save
ggsave(filename = here("figures", "Figure6abc_betaJTU_FTWS_models.png"),
       plot = final_plot, width = 14, height = 7, dpi = 300, bg = "white")

ggsave(filename = here("figures", "Figure6abc_betaJTU_FTWS_models.pdf"),
       plot = final_plot, width = 14, height = 7, bg = "white")

# 6. PREDICTION PLOTS ----------------------------------------------------------

## 6.1. All occurrences --------------------------------------------------------

# Get the data and model coefficients
data_clean <- turnover_forest_tws_15km_coords_time
coefs <- coef(FTWS_turnover_model5_gls_log)

# Define colors
line_color <- "#2196F3"

# Plot 1: forest_to_tws 

# Get range of values
pred_range_1 <- seq(min(data_clean$forest_to_tws, na.rm = TRUE),
                    max(data_clean$forest_to_tws, na.rm = TRUE),
                    length.out = 100)

# Manual prediction calculation
pred_1 <- coefs["(Intercept)"] + 
  coefs["forest_to_tws"] * pred_range_1 +
  coefs["forest_no_change"] * mean(data_clean$forest_no_change, na.rm = TRUE) +
  coefs["delta_recorder_effort"] * mean(data_clean$delta_recorder_effort, na.rm = TRUE) +
  coefs["log_recorder_effort"] * mean(data_clean$log_recorder_effort, na.rm = TRUE) +
  coefs["temp_change"] * mean(data_clean$temp_change, na.rm = TRUE) +
  coefs["precip_change"] * mean(data_clean$precip_change, na.rm = TRUE) +
  coefs["lc_time_period2006-2012"] * 1  # Reference level = 2006-2012

# Plot prediction for forest to tws
plot1 <- ggplot(data.frame(x = pred_range_1, y = pred_1), aes(x = x, y = y)) +
  geom_line(color = line_color, linewidth = 1.2) +
  labs(x = "Forest → TWS (proportion)", y = "Predicted beta_jtu", title = "a)") +
  theme_classic() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3))

# Plot 2: forest_no_change 

# Get range of values for forest_no_change 
pred_range_2 <- seq(min(data_clean$forest_no_change, na.rm = TRUE),
                    max(data_clean$forest_no_change, na.rm = TRUE),
                    length.out = 100)

# Manual calculation of predictors
pred_2 <- coefs["(Intercept)"] + 
  coefs["forest_to_tws"] * mean(data_clean$forest_to_tws, na.rm = TRUE) +
  coefs["forest_no_change"] * pred_range_2 +
  coefs["delta_recorder_effort"] * mean(data_clean$delta_recorder_effort, na.rm = TRUE) +
  coefs["log_recorder_effort"] * mean(data_clean$log_recorder_effort, na.rm = TRUE) +
  coefs["temp_change"] * mean(data_clean$temp_change, na.rm = TRUE) +
  coefs["precip_change"] * mean(data_clean$precip_change, na.rm = TRUE) +
  coefs["lc_time_period2006-2012"] * 1

# Plot prediction for forest no change
plot2 <- ggplot(data.frame(x = pred_range_2, y = pred_2), aes(x = x, y = y)) +
  geom_line(color = line_color, linewidth = 1.2) +
  labs(x = "Forest no change (proportion)", y = NULL, title = "b)") +
  theme_classic() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3))

# Plot 3: log_recorder_effort 

# Get range of values
pred_range_3 <- seq(min(data_clean$log_recorder_effort, na.rm = TRUE),
                    max(data_clean$log_recorder_effort, na.rm = TRUE),
                    length.out = 100)

# Manual calculation of predictor
pred_3 <- coefs["(Intercept)"] + 
  coefs["forest_to_tws"] * mean(data_clean$forest_to_tws, na.rm = TRUE) +
  coefs["forest_no_change"] * mean(data_clean$forest_no_change, na.rm = TRUE) +
  coefs["delta_recorder_effort"] * mean(data_clean$delta_recorder_effort, na.rm = TRUE) +
  coefs["log_recorder_effort"] * pred_range_3 +
  coefs["temp_change"] * mean(data_clean$temp_change, na.rm = TRUE) +
  coefs["precip_change"] * mean(data_clean$precip_change, na.rm = TRUE) +
  coefs["lc_time_period2006-2012"] * 1

# Plot predictor for log_recorder_effort
plot3 <- ggplot(data.frame(x = pred_range_3, y = pred_3), aes(x = x, y = y)) +
  geom_line(color = line_color, linewidth = 1.2) +
  labs(x = "log(Recorder effort)", y = NULL, title = "c)") +
  theme_classic() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3))

# Plot 4: delta_recorder_effort 

# Get range of values for delta recorder effort
pred_range_4 <- seq(min(data_clean$delta_recorder_effort, na.rm = TRUE),
                    max(data_clean$delta_recorder_effort, na.rm = TRUE),
                    length.out = 100)

# Manual calculation of predictor
pred_4 <- coefs["(Intercept)"] + 
  coefs["forest_to_tws"] * mean(data_clean$forest_to_tws, na.rm = TRUE) +
  coefs["forest_no_change"] * mean(data_clean$forest_no_change, na.rm = TRUE) +
  coefs["delta_recorder_effort"] * pred_range_4 +
  coefs["log_recorder_effort"] * mean(data_clean$log_recorder_effort, na.rm = TRUE) +
  coefs["temp_change"] * mean(data_clean$temp_change, na.rm = TRUE) +
  coefs["precip_change"] * mean(data_clean$precip_change, na.rm = TRUE) +
  coefs["lc_time_period2006-2012"] * 1

# Plot predicted values for delta recorder effort
plot4 <- ggplot(data.frame(x = pred_range_4, y = pred_4), aes(x = x, y = y)) +
  geom_line(color = line_color, linewidth = 1.2) +
  labs(x = "ΔRecorder effort", y = "Predicted beta_jtu", title = "d)") +
  theme_classic() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3))

# Plot 5: temp_change 

# Get range of temperature change values
pred_range_5 <- seq(min(data_clean$temp_change, na.rm = TRUE),
                    max(data_clean$temp_change, na.rm = TRUE),
                    length.out = 100)

# Manual calculation of predictors
pred_5 <- coefs["(Intercept)"] + 
  coefs["forest_to_tws"] * mean(data_clean$forest_to_tws, na.rm = TRUE) +
  coefs["forest_no_change"] * mean(data_clean$forest_no_change, na.rm = TRUE) +
  coefs["delta_recorder_effort"] * mean(data_clean$delta_recorder_effort, na.rm = TRUE) +
  coefs["log_recorder_effort"] * mean(data_clean$log_recorder_effort, na.rm = TRUE) +
  coefs["temp_change"] * pred_range_5 +
  coefs["precip_change"] * mean(data_clean$precip_change, na.rm = TRUE) +
  coefs["lc_time_period2006-2012"] * 1

# Plot predicted values for temperature change
plot5 <- ggplot(data.frame(x = pred_range_5, y = pred_5), aes(x = x, y = y)) +
  geom_line(color = line_color, linewidth = 1.2) +
  labs(x = "Temperature change (°C)", y = NULL, title = "e)") +
  theme_classic() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3))

# Plot 6: precip_change

# Get range of values for precipitation changes
pred_range_6 <- seq(min(data_clean$precip_change, na.rm = TRUE),
                    max(data_clean$precip_change, na.rm = TRUE),
                    length.out = 100)

# Manual calculation of predictions
pred_6 <- coefs["(Intercept)"] + 
  coefs["forest_to_tws"] * mean(data_clean$forest_to_tws, na.rm = TRUE) +
  coefs["forest_no_change"] * mean(data_clean$forest_no_change, na.rm = TRUE) +
  coefs["delta_recorder_effort"] * mean(data_clean$delta_recorder_effort, na.rm = TRUE) +
  coefs["log_recorder_effort"] * mean(data_clean$log_recorder_effort, na.rm = TRUE) +
  coefs["temp_change"] * mean(data_clean$temp_change, na.rm = TRUE) +
  coefs["precip_change"] * pred_range_6 +
  coefs["lc_time_period2006-2012"] * 1

# Plot predicted values for precipitation change
plot6 <- ggplot(data.frame(x = pred_range_6, y = pred_6), aes(x = x, y = y)) +
  geom_line(color = line_color, linewidth = 1.2) +
  labs(x = "Precipitation change (mm/yr)", y = NULL, title = "f)") +
  theme_classic() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3))

# Combine all plots 

combined_plot <- (plot1 | plot2 | plot3) / (plot4 | plot5 | plot6)

# Display
print(combined_plot)

# 7. GLS WITH DIFFERENT RECORDING EFFORT MEASUREMENT ---------------------------

# Since delta recording effort seems to have a +ve effect on the turnover while
# recording effort (total_occ_before + total_occ_after) has a -ve effect on turnover
# which seems like a counter intuitive thing, I decided to run an additional model
# where recorder effort is calculated as mean(total_occ_before, total_occ_after)
# and check if the results are the same

# Calculate mean recorder effort
turnover_forest_tws_15km_coords_time <- turnover_forest_tws_15km_coords_time |>
  mutate(recorder_effort_mean = (total_occ_before + total_occ_after) / 2,
         log_recorder_effort_mean = log(recorder_effort_mean))

# Check for any issues (zeros, NAs)
summary(turnover_forest_tws_15km_coords_time$recorder_effort_mean)
any(!is.finite(turnover_forest_tws_15km_coords_time$log_recorder_effort_mean))

# Fit the model with mean recorder effort
FTWS_turnover_model_mean_effort <- gls(beta_jtu ~ forest_to_tws + forest_no_change +
                                         delta_recorder_effort + log_recorder_effort_mean + 
                                         lc_time_period + temp_change + precip_change,
                                       correlation = corExp(form = ~ x + y | time_numeric),
                                       data = turnover_forest_tws_15km_coords_time,
                                       method = "REML")

# Compare with original model
AICtab(FTWS_turnover_model5_gls_log, FTWS_turnover_model_mean_effort, base = TRUE)

# Look at the summaries side by side
summary(FTWS_turnover_model5_gls_log)
summary(FTWS_turnover_model_mean_effort)

# Compare coefficients
coef_sum <- coef(FTWS_turnover_model5_gls_log)
coef_mean <- coef(FTWS_turnover_model_mean_effort)

comparison <- data.frame(predictor = names(coef_sum),
                         sum_model = coef_sum,
                         mean_model = coef_mean,
                         difference = coef_mean - coef_sum)
# END OF SCRIPT ----------------------------------------------------------------