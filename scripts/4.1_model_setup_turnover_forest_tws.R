##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 4.1_model_setup_turnover_forest_tws
# This script contains code which sets up the models exploring the impact of 
# Forest -> TWS land cover transition on temporal turnover for all occurrences,
# plant-only occurrences and bird-only occurrences
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# Turnover data for all occurrences
load(here("data", "derived_data", 
          "all_periods_turnover_all_land_cover_chanegs_15km.rda"))

# Turnover data for plants
load(here("data", "derived_data", 
          "vascular_plants_all_periods_turnover_all_land_cover_chanegs_15km.rda"))

# Turnover data for birds
load(here("data", "derived_data", 
          "birds_all_periods_turnover_all_land_cover_chanegs_15km.rda"))

# 2. ALL OCCURRENCES -----------------------------------------------------------

## 2.1. Prepare data for analysis ----------------------------------------------

# Select only Forest -> TWS columns 
# check column names
colnames(all_periods_turnover_all_land_cover_chanegs_15km)

# Also rename the columns for easier manipulation of df
turnover_forest_tws_15km <- all_periods_turnover_all_land_cover_chanegs_15km |>
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
  # remove columns no longer required
  select(-`2000-2006_Forest no change`, -`2006-2012_Forest no change`, 
         -`2012-2018_Forest no change`,-`2000-2006_Forest to TWS`, 
         -`2006-2012_Forest to TWS`, -`2012-2018_Forest to TWS`)

# Check transformation went ok
head(turnover_forest_tws_15km)

# Prepare data for GLS: remove rows with missing x or y and categorise time periods
turnover_forest_tws_15km_coords_time <- turnover_forest_tws_15km |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3))

## 2.2. GLS with raw data ---------------------------------------------------------

# Fit gls
model1_gls <- gls(JDI ~ forest_to_tws + forest_no_change + 
                    delta_recorder_effort + recorder_effort + lc_time_period,
                  correlation = corExp(form = ~ x + y | time_numeric),
                  data = turnover_forest_tws_15km_coords_time,
                  method = "REML")

# Save model output to file
save(model1_gls, 
     file = here("data", "models", 
                 "gls_model1_all_occurrences_turnover_results.RData"))

# Get summary
summary(model1_gls)

# Extract correlation structure parameters
print(model1_gls$modelStruct$corStruct)

# Get range parameter
range_param <- coef(model1_gls$modelStruct$corStruct, unconstrained = FALSE)

# Extract residuals
residuals_gls <- residuals(model1_gls, type = "normalized")

# Basic residual plots
par(mfrow = c(2, 2))

# Residuals vs fitted
plot(fitted(model1_gls), residuals_gls,
     xlab = "Fitted Values", ylab = "Normalized Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# QQ plot
qqnorm(residuals_gls, main = "Normal Q-Q Plot")
qqline(residuals_gls, col = "red")

# Residuals vs Forest to TWS 
plot(turnover_forest_tws_15km_coords_time$forest_to_tws, residuals_gls,
     xlab = "Forest to TWS", ylab = "Normalized Residuals",
     main = "Residuals vs Forest to TWS")
abline(h = 0, col = "red", lty = 2)

# Residuals vs Recorder Effort
plot(turnover_forest_tws_15km_coords_time$recorder_effort, residuals_gls,
     xlab = "Recorder Effort", ylab = "Normalized Residuals",
     main = "Residuals vs Recorder Effort")
abline(h = 0, col = "red", lty = 2)

## 2.3. GLS with logged recorder effort ----------------------------------------

# Check if there are any cells with recorder effort = 0
a <- turnover_forest_tws_15km_coords_time |>
  filter(recorder_effort == 0)
length(a) #0 - Good!

# Log transform recorder effort values
turnover_forest_tws_15km_coords_time <- turnover_forest_tws_15km_coords_time |>
  mutate(log_recorder_effort = log(recorder_effort))

# Check log transformed values
summary(turnover_forest_tws_15km_coords_time$recorder_effort)
any(!is.finite(turnover_forest_tws_15km_coords_time$recorder_effort)) # FALSE = no infinite values - Good!


# Define GLS
model2_gls <- gls(JDI ~ forest_to_tws + forest_no_change + 
                    delta_recorder_effort + log_recorder_effort + lc_time_period,
                  correlation = corExp(form = ~ x + y | time_numeric),  
                  data = turnover_forest_tws_15km_coords_time,
                  method = "REML")

# Save model output to file
save(model2_gls, 
     file = here("data", "models", 
                 "gls_model2_all_occurrences_turnover_results.RData"))

# Get summary
model2_gls_summary <- summary(model2_gls)

# Extract correlation structure parameters
print(model2_gls$modelStruct$corStruct)

# Get range parameter
range_param <- coef(model2_gls$modelStruct$corStruct, unconstrained = FALSE)

# Extract residuals
residuals_gls <- residuals(model2_gls, type = "normalized")

# Basic residual plots
par(mfrow = c(2, 2))

# Residuals vs fitted
plot(fitted(model2_gls), residuals_gls,
     xlab = "Fitted Values", ylab = "Normalized Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# QQ plot
qqnorm(residuals_gls, main = "Normal Q-Q Plot")
qqline(residuals_gls, col = "red")

# Residuals vs predictors
plot(turnover_forest_tws_15km_coords_time$forest_to_tws, residuals_gls,
     xlab = "Forest to TWS", ylab = "Normalized Residuals",
     main = "Residuals vs Forest to TWS")
abline(h = 0, col = "red", lty = 2)

plot(turnover_forest_tws_15km_coords_time$log_recorder_effort, residuals_gls,
     xlab = "Recorder Effort", ylab = "Normalized Residuals",
     main = "Residuals vs Recorder Effort")
abline(h = 0, col = "red", lty = 2)

## 2.4. Plot model output ------------------------------------------------------

# Create dataframe of coeficients
model2_coef_df <- data.frame(term = names(model_summary$tTable[, "Value"]),
                             estimate = model_summary$tTable[, "Value"],
                             std.error = model_summary$tTable[, "Std.Error"],
                             statistic = model_summary$tTable[, "t-value"],
                             p.value = model_summary$tTable[, "p-value"])

# Remove the intercept
model2_coef_df_no_intercept <- model2_coef_df[model2_coef_df$term != "(Intercept)", ]

# Create coefficient plot
figure6_a <- ggplot(model2_coef_df_no_intercept, aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_y_discrete(labels = c("forest_to_tws" = "Forest to TWS",
                              "forest_no_change" = "Forest No Change", 
                              "delta_recorder_effort" = "ΔRecorder Effort",
                              "log_recorder_effort" = "log Recorder Effort",
                              "lc_time_period2006-2012" = "2006-2012 Time Period",
                              "lc_time_period2012-2018" = "2012-2018 Time Period"),
                   limits = c("lc_time_period2012-2018",
                              "lc_time_period2006-2012", 
                              "log_recorder_effort",
                              "delta_recorder_effort",
                              "forest_no_change",
                              "forest_to_tws")) +
  labs(x = "Estimate ± 95% CI", y = NULL) +
  theme_classic()

# Save figure as .png
ggsave(filename = here("figures", "Figure6a_gls_model_output_all_occurrences_turnover.png"),
       width = 12, height = 8, dpi = 300)

# Save figure as .svg
ggsave(filename = here("figures", "Figure6a_gls_model_output_all_occurrences_turnover.svg"),
       width = 12, height = 8, dpi = 300)

## 2.5. GAM with delta recorder effort and spatial smooth ----------------------

# Convert lc_time_period to factor
turnover_forest_tws_15km_coords_time <- turnover_forest_tws_15km_coords_time |>
  mutate(lc_time_period = as.factor(lc_time_period))

# Modify JDI values so that they do not touch [0,1] - to fit a beta GAM
# get N
N <- nrow(turnover_forest_tws_15km_coords_time)
# calculate new JDI values
turnover_forest_tws_15km_coords_time <- turnover_forest_tws_15km_coords_time |>
  mutate(JDI_beta = (JDI * (N - 1) + 0.5) / N)

# Define model
model5_gam <- gam(JDI_beta ~ forest_to_tws + forest_no_change + s(delta_recorder_effort) +
                           log_recorder_effort + lc_time_period +
                           s(x, y, by = lc_time_period, k = 120),
                         data = turnover_forest_tws_15km_coords_time,
                         family = betar(link = "logit"),
                         method = "REML")

# View model summary
summary(model5_gam)

# Save model output
# save(plants_model5_gam, 
#      file = here("data", "models", 
#                  "gam_model5_plant_occurrences_turnover_forest_tws_results.RData"))

# Run diagnostics
par(mfrow = c(2, 2))
gam.check(model5_gam)

# 3. PLANT OCCURRENCES ONLY ----------------------------------------------------

## 3.1. Prepare plant data for analysis ----------------------------------------

# Select only Forest -> TWS columns
# check column names
colnames(vascular_plants_all_periods_turnover_all_land_cover_chanegs_15km)

# Also rename the columns for easier manipulation of df
plants_turnover_forest_tws_15km <- vascular_plants_all_periods_turnover_all_land_cover_chanegs_15km |>
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
plants_model3_gls <- gls(JDI ~ forest_to_tws + forest_no_change + 
                    delta_recorder_effort + log_recorder_effort + lc_time_period,
                  correlation = corExp(form = ~ x + y | time_numeric),  
                  data = plants_turnover_forest_tws_15km_coords_time,
                  method = "REML")

# Save model output to file
save(plants_model3_gls, 
     file = here("data", "models", 
                 "gls_model3_plant_occurrences_turnover_forest_tws_results.RData"))

# Get summary
plants_model3_gls_summary <- summary(plants_model3_gls)

# Extract correlation structure parameters
print(plants_model3_gls$modelStruct$corStruct)

# Get range parameter
range_param <- coef(plants_model3_gls$modelStruct$corStruct, unconstrained = FALSE)

# Extract residuals
plant_residuals_gls <- residuals(plants_model3_gls, type = "normalized")

# Basic residual plots
par(mfrow = c(2, 2))

# Residuals vs fitted
plot(fitted(plants_model3_gls), plant_residuals_gls,
     xlab = "Fitted Values", ylab = "Normalized Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# QQ plot
qqnorm(plant_residuals_gls, main = "Normal Q-Q Plot")
qqline(plant_residuals_gls, col = "red")

# Residuals vs predictors
plot(plants_turnover_forest_tws_15km_coords_time$forest_to_tws, plant_residuals_gls,
     xlab = "Forest to TWS", ylab = "Normalized Residuals",
     main = "Residuals vs Forest to TWS")
abline(h = 0, col = "red", lty = 2)

plot(plants_turnover_forest_tws_15km_coords_time$log_recorder_effort, plant_residuals_gls,
     xlab = "Recorder Effort", ylab = "Normalized Residuals",
     main = "Residuals vs Recorder Effort")
abline(h = 0, col = "red", lty = 2)

## 3.3. GAM with spatial smooth ------------------------------------------------

# Convert lc_time_period to factor
plants_turnover_forest_tws_15km_coords_time <- plants_turnover_forest_tws_15km_coords_time |>
  mutate(lc_time_period = as.factor(lc_time_period))

# Modify JDI values so that they do not touch [0,1] - to fit a beta GAM
# get N
N <- nrow(plants_turnover_forest_tws_15km_coords_time)
# calculate new JDI values
plants_turnover_forest_tws_15km_coords_time <- plants_turnover_forest_tws_15km_coords_time |>
  mutate(JDI_beta = (JDI * (N - 1) + 0.5) / N)

# Fit GAM with spatial smooth that is separate per time period
plants_model4_gam <- gam(JDI_beta ~ forest_to_tws + forest_no_change + delta_recorder_effort +
                           log_recorder_effort + lc_time_period +
                           s(x, y, by = lc_time_period, k = 120),
                         data = plants_turnover_forest_tws_15km_coords_time,
                         family = betar(link = "logit"),
                         method = "REML")

# View model summary
summary(plants_model4_gam)

# Save model output
save(plants_model4_gam, 
     file = here("data", "models", 
                 "gam_model4_plant_occurrences_turnover_forest_tws_results.RData"))

# Run diagnostics
par(mfrow = c(2, 2))
gam.check(plants_model4_gam)

## 3.4. GAM with forest -> tws and spatial smooth ------------------------------

# Define model
plants_model5_gam <- gam(JDI_beta ~ s(forest_to_tws) + forest_no_change + delta_recorder_effort +
                           log_recorder_effort + lc_time_period +
                           s(x, y, by = lc_time_period, k = 120),
                         data = plants_turnover_forest_tws_15km_coords_time,
                         family = betar(link = "logit"),
                         method = "REML")

# View model summary
summary(plants_model5_gam)

# Save model output
save(plants_model5_gam, 
     file = here("data", "models", 
                 "gam_model5_plant_occurrences_turnover_forest_tws_results.RData"))

# Run diagnostics
par(mfrow = c(2, 2))
gam.check(plants_model5_gam)

# 4. BIRD OCCURRENCES ONLY -----------------------------------------------------

## 4.1. Prepare bird data for analysis -----------------------------------------

# Check column names
colnames(birds_all_periods_turnover_all_land_cover_chanegs_15km)

# Also rename the columns for easier manipulation of df
birds_turnover_forest_tws_15km <- birds_all_periods_turnover_all_land_cover_chanegs_15km |>
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
length(a) #0 - Good!

# Log transform recorder effort values
birds_turnover_forest_tws_15km_coords_time <- birds_turnover_forest_tws_15km_coords_time |>
  mutate(log_recorder_effort = log(recorder_effort))

# Check log transformed values
summary(birds_turnover_forest_tws_15km_coords_time$recorder_effort)
any(!is.finite(birds_turnover_forest_tws_15km_coords_time$recorder_effort)) # FALSE = no infinite values - Good!

# Define GLS
birds_model6_gls <- gls(JDI ~ forest_to_tws + forest_no_change + 
                           delta_recorder_effort + log_recorder_effort + lc_time_period,
                         correlation = corExp(form = ~ x + y | time_numeric),  
                         data = birds_turnover_forest_tws_15km_coords_time,
                         method = "REML")

# Save model output to file
save(birds_model6_gls, 
     file = here("data", "models", 
                 "gls_model6_birds_occurrences_turnover_forest_tws_results.RData"))

# Get summary
birds_model6_gls_summary <- summary(birds_model6_gls)

# Extract correlation structure parameters
print(birds_model6_gls$modelStruct$corStruct)

# Get range parameter
birds_range_param <- coef(birds_model6_gls$modelStruct$corStruct, unconstrained = FALSE)

# Extract residuals
birds_residuals_gls <- residuals(birds_model6_gls, type = "normalized")

# Basic residual plots
par(mfrow = c(2, 2))

# Residuals vs fitted
plot(fitted(birds_model6_gls), birds_residuals_gls,
     xlab = "Fitted Values", ylab = "Normalized Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# QQ plot
qqnorm(birds_residuals_gls, main = "Normal Q-Q Plot")
qqline(birds_residuals_gls, col = "red")

# Residuals vs predictors
plot(birds_turnover_forest_tws_15km_coords_time$forest_to_tws, birds_residuals_gls,
     xlab = "Forest to TWS", ylab = "Normalized Residuals",
     main = "Residuals vs Forest to TWS")
abline(h = 0, col = "red", lty = 2)

plot(birds_turnover_forest_tws_15km_coords_time$log_recorder_effort, birds_residuals_gls,
     xlab = "Recorder Effort", ylab = "Normalized Residuals",
     main = "Residuals vs Recorder Effort")
abline(h = 0, col = "red", lty = 2)

## 4.3. GAM with spatial smooth ------------------------------------------------

# Convert lc_time_period to factor
birds_turnover_forest_tws_15km_coords_time <- birds_turnover_forest_tws_15km_coords_time |>
  mutate(lc_time_period = as.factor(lc_time_period))

# Modify JDI values so that they do not touch [0,1] - to fit a beta GAM
# get N
N <- nrow(birds_turnover_forest_tws_15km_coords_time)
# calculate new JDI values
birds_turnover_forest_tws_15km_coords_time <- birds_turnover_forest_tws_15km_coords_time |>
  mutate(JDI_beta = (JDI * (N - 1) + 0.5) / N)

# Fit GAM with spatial smooth that is separate per time period
birds_model7_gam <- gam(JDI_beta ~ forest_to_tws + forest_no_change + delta_recorder_effort +
                           log_recorder_effort + lc_time_period +
                           s(x, y, by = lc_time_period, k = 120),
                         data = birds_turnover_forest_tws_15km_coords_time,
                         family = betar(link = "logit"),
                         method = "REML")

# View model summary
summary(birds_model7_gam)

# Save model output
save(birds_model7_gam, 
     file = here("data", "models", 
                 "gam_model7_birds_occurrences_turnover_forest_tws_results.RData"))

# Run diagnostics
par(mfrow = c(2, 2))
gam.check(birds_model7_gam)

# END OF SCRIPT ----------------------------------------------------------------