##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 4.1_model_setup_turnover_forest_tws
# This script contains code which sets up the models exploring the imapct of 
# Forest -> TWS land cover transition on temporal turnover
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# Turnover data
load(here("data", "derived_data", 
          "all_periods_turnover_all_land_cover_chanegs_15km.rda"))

# 2. PREPARE DATA FOR ANALYSIS -------------------------------------------------

## 2.1. Select only Forest -> TWS columns --------------------------------------

# Check column names
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

## 2.2. Prepare data for GLS ---------------------------------------------------

# Remove rows that might have NA for x or y
turnover_forest_tws_15km_coords <- turnover_forest_tws_15km |>
  filter(!is.na(x) & !is.na(y))

# Categorise time periods
turnover_forest_tws_15km_coords_time <- turnover_forest_tws_15km_coords |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3))

# 3. GLS with raw data ---------------------------------------------------------

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

# 4. GLS WITH LOG RECORDER EFFORT ----------------------------------------------

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
summary(model2_gls)

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

# END OF SCRIPT ----------------------------------------------------------------