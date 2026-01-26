##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 5.1_model_setup_turnover_forest_tws
# This script contains code which sets up the models exploring the impact of 
# Forest -> TWS land cover transition on temporal turnover for all occurrences,
# plant-only occurrences and bird-only occurrences
##----------------------------------------------------------------------------##

library(here)
source(here("scripts", "0_setup.R"))

# 1. LOAD DATA -----------------------------------------------------------------

# Turnover data for all occurrences
load(here("data", "derived_data", 
          "all_periods_turnover_all_land_cover_climate_15km.rda"))

# Turnover data for plants
load(here("data", "derived_data", 
          "vascular_plants_turnover_all_land_cover_climate_15km.rda"))

# Turnover data for birds
load(here("data", "derived_data", 
          "bird_turnover_all_land_cover_climate_15km.rda"))

# 2. ALL OCCURRENCES -----------------------------------------------------------

## 2.1. Prepare data for analysis ----------------------------------------------

# Select only Forest -> TWS columns
# check column names
colnames(all_periods_turnover_with_climate)

# Also rename the columns for easier manipulation of df
turnover_forest_tws_15km <- all_periods_turnover_with_climate |>
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

## 2.2. Beta regression --------------------------------------------------------

# Get N
N <- nrow(turnover_forest_tws_15km)

# Transform JDI values so they do not touch [0,1] anymore
turnover_forest_tws_15km <- turnover_forest_tws_15km |>
  mutate(JDI_beta = (JDI * (N - 1) + 0.5) / N)

# Run model
FTWS_turnover_model1_beta_regression <- betareg::betareg(JDI_beta ~ forest_to_tws + forest_no_change +
                               delta_recorder_effort + recorder_effort + lc_time_period,
                             data = turnover_forest_tws_15km)

# Save model
save(FTWS_turnover_model1_beta_regression,
     file = here("data", "models", "exploratory",
                 "FTWS_turnover_model1_beta_regression.RData"))


## 2.3. GLMM with SSB ID as random effect --------------------------------------

# Run model
FTWS_turnover_model2_GLMM <- glmmTMB(JDI_beta ~ forest_to_tws + forest_no_change +
                            delta_recorder_effort + recorder_effort + lc_time_period +
                            (1|ssb_id),
                          family = beta_family(),
                          data = turnover_forest_tws_15km)

# Save model
save(FTWS_turnover_model2_GLMM,
     file = here("data", "models", "exploratory",
                 "FTWS_turnover_model2_GLMM.RData"))

## 2.4. Ordered beta regression ------------------------------------------------

# Run model
# FTWS_turnover_model3_ordered_beta <- ordbetareg(JDI ~ forest_to_tws + forest_no_change +
#                                                   delta_recorder_effort + recorder_effort + lc_time_period +
#                                                   (1|ssb_id),
#                                                 data = turnover_forest_tws_15km,
#                                                 cores = 4,
#                                                 iter = 4000, # double the iterations
#                                                 chains = 4,
#                                                 control = list(adapt_delta = 0.95,  # more conservative sampling
#                                                                max_treedepth = 12))

# Save model
# save(FTWS_turnover_model3_ordered_beta,
#      file = here("data", "models", "exploratory",
#                  "FTWS_turnover_model3_ordered_beta.RData"))

## 2.5. GLS with raw data and exponential structure ----------------------------

# Fit gls
FTWS_turnover_model4_gls_raw <- gls(JDI ~ forest_to_tws + forest_no_change +
                    delta_recorder_effort + recorder_effort + lc_time_period,
                  correlation = corExp(form = ~ x + y | time_numeric),
                  data = turnover_forest_tws_15km_coords_time,
                  method = "REML")

# Save model output to file
save(FTWS_turnover_model4_gls_raw,
     file = here("data", "models", "exploratory",
                 "FTWS_turnover_model4_gls_raw.RData"))


## 2.6. GLS with logged recorder effort ----------------------------------------

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
FTWS_turnover_model5_gls_log <- gls(JDI ~ forest_to_tws + forest_no_change +
                    delta_recorder_effort + log_recorder_effort + lc_time_period + temp_change + precip_change,
                  correlation = corExp(form = ~ x + y | time_numeric),
                  data = turnover_forest_tws_15km_coords_time,
                  method = "REML")

# Model validation looks ok => THIS IS THE FINAL MODEL WE WILL USE

# Save model output to file
save(FTWS_turnover_model5_gls_log,
     file = here("data", "models", "final",
                 "FTWS_turnover_model5_gls_log.RData"))

## 2.4. Plot model output ------------------------------------------------------

# Get summary of the model
model_summary <- summary(FTWS_turnover_model5_gls_log)

# Create dataframe of coeficients
FTWS_turnover_model5_gls_log_coef_df <- data.frame(term = names(model_summary$tTable[, "Value"]),
                                                   estimate = model_summary$tTable[, "Value"],
                                                   std.error = model_summary$tTable[, "Std.Error"],
                                                   statistic = model_summary$tTable[, "t-value"],
                                                   p.value = model_summary$tTable[, "p-value"])

# Remove the intercept
FTWS_turnover_model5_gls_log_coef_df_no_intercept <- FTWS_turnover_model5_gls_log_coef_df[FTWS_turnover_model5_gls_log_coef_df$term != "(Intercept)", ]

# Create coefficient plot
figure6_a <- ggplot(FTWS_turnover_model5_gls_log_coef_df_no_intercept, aes(x = estimate, y = term)) +
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
plants_FTWS_model1_gls <- gls(JDI ~ forest_to_tws + forest_no_change +
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

## 3.3. Plot model output ------------------------------------------------------

# Get model summary
plants_FTWS_model1_gls_summary <- summary(plants_FTWS_model1_gls)

# Create dataframe of coeficients
plants_FTWS_model1_gls_coef_df <- data.frame(term = names(plants_FTWS_model1_gls_summary$tTable[, "Value"]),
                                             estimate = plants_FTWS_model1_gls_summary$tTable[, "Value"],
                                             std.error = plants_FTWS_model1_gls_summary$tTable[, "Std.Error"],
                                             statistic = plants_FTWS_model1_gls_summary$tTable[, "t-value"],
                                             p.value = plants_FTWS_model1_gls_summary$tTable[, "p-value"])

# Remove the intercept
plants_FTWS_model1_gls_coef_df_no_intercept <- plants_FTWS_model1_gls_coef_df[plants_FTWS_model1_gls_coef_df$term != "(Intercept)", ]

# Create coefficient plot
figure6_b <- ggplot(plants_FTWS_model1_gls_coef_df_no_intercept,
                    aes(x = estimate, y = term)) +
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
ggsave(filename = here("figures", "Figure6b_gls_model_output_plants_turnover.png"),
       width = 12, height = 8, dpi = 300)

# Save figure as .svg
ggsave(filename = here("figures", "Figure6b_gls_model_output_plants_turnover.svg"),
       width = 12, height = 8, dpi = 300)

## 3.4. GAM with spatial smooth ------------------------------------------------

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
plants_FTWS_model2_GAM <- gam(JDI_beta ~ forest_to_tws + forest_no_change + delta_recorder_effort +
                           log_recorder_effort + lc_time_period +
                           s(x, y, by = lc_time_period, k = 120),
                         data = plants_turnover_forest_tws_15km_coords_time,
                         family = betar(link = "logit"),
                         method = "REML")

# Save model output
save(plants_FTWS_model2_GAM,
     file = here("data", "models", "exploratory",
                 "plants_FTWS_model2_GAM.RData"))

## 3.5. GAM with forest -> tws and spatial smooth ------------------------------

# Define model
plants_FTWS_model3_GAM_extra_smooth <- gam(JDI_beta ~ s(forest_to_tws) + forest_no_change + delta_recorder_effort +
                                             log_recorder_effort + lc_time_period +
                                             s(x, y, by = lc_time_period, k = 120),
                                           data = plants_turnover_forest_tws_15km_coords_time,
                                           family = betar(link = "logit"),
                                           method = "REML")


# Save model output
save(plants_FTWS_model3_GAM_extra_smooth,
     file = here("data", "models", "exploratory",
                 "plants_FTWS_model3_GAM_extra_smooth.RData"))

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
birds_FTWS_model1_gls <- gls(JDI ~ forest_to_tws + forest_no_change +
                           delta_recorder_effort + log_recorder_effort + lc_time_period + temp_change + precip_change,
                         correlation = corExp(form = ~ x + y | time_numeric),
                         data = birds_turnover_forest_tws_15km_coords_time,
                         method = "REML")

# Save model output to file
# save(birds_FTWS_model1_gls, 
#      file = here("data", "models", "final",
#                  "birds_FTWS_model1_gls.RData"))

## 4.3. Plot model output ------------------------------------------------------

# Get summary
# model_summary_birds <- summary(birds_FTWS_model1_gls)
# 
# # Create dataframe of coeficients
# birds_FTWS_model1_gls_coef_df <- data.frame(term = names(model_summary_birds$tTable[, "Value"]),
#                                             estimate = model_summary_birds$tTable[, "Value"],
#                                             std.error = model_summary_birds$tTable[, "Std.Error"],
#                                             statistic = model_summary_birds$tTable[, "t-value"],
#                                             p.value = model_summary_birds$tTable[, "p-value"])
# 
# # Remove the intercept
# birds_FTWS_model1_gls_coef_df_no_intercept <- birds_FTWS_model1_gls_coef_df[birds_FTWS_model1_gls_coef_df$term != "(Intercept)", ]
# 
# # Create coefficient plot
# figure6_c <- ggplot(birds_FTWS_model1_gls_coef_df_no_intercept, aes(x = estimate, y = term)) +
#   geom_point() +
#   geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
#                      xmax = estimate + 1.96 * std.error)) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   scale_y_discrete(labels = c("forest_to_tws" = "Forest to TWS",
#                               "forest_no_change" = "Forest No Change", 
#                               "delta_recorder_effort" = "ΔRecorder Effort",
#                               "log_recorder_effort" = "log Recorder Effort",
#                               "lc_time_period2006-2012" = "2006-2012 Time Period",
#                               "lc_time_period2012-2018" = "2012-2018 Time Period"),
#                    limits = c("lc_time_period2012-2018",
#                               "lc_time_period2006-2012", 
#                               "log_recorder_effort",
#                               "delta_recorder_effort",
#                               "forest_no_change",
#                               "forest_to_tws")) +
#   labs(x = "Estimate ± 95% CI", y = NULL) +
#   theme_classic()
# 
# # Display plot
# figure6_c
# 
# # Save figure as .png
# ggsave(filename = here("figures", "Figure6c_gls_model_output_birds_turnover.png"),
#        width = 12, height = 8, dpi = 300)
# 
# # Save figure as .svg
# ggsave(filename = here("figures", "Figure6c_gls_model_output_birds_turnover.svg"),
#        width = 12, height = 8, dpi = 300)

## 4.4. Weighted GLS -----------------------------------------------------------

# Define model
birds_FTWS_model2_weighted_gls <- gls(JDI ~ forest_to_tws + forest_no_change +
                                        delta_recorder_effort + log_recorder_effort + lc_time_period,
                                      weights = varPower(form = ~ fitted(.)),
                                      correlation = corExp(form = ~ x + y | time_numeric),
                                      data = birds_turnover_forest_tws_15km_coords_time,
                                      method = "REML")

# Save model
save(birds_FTWS_model2_weighted_gls, 
     file = here("data", "models", "exploratory",
                 "birds_FTWS_model2_weighted_gls.RData"))

## 4.5. GAM with spatial smooth ------------------------------------------------

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
birds_FTWS_model3_GAM <- gam(JDI_beta ~ forest_to_tws + forest_no_change + delta_recorder_effort +
                           log_recorder_effort + lc_time_period +
                           s(x, y, by = lc_time_period, k = 120),
                         data = birds_turnover_forest_tws_15km_coords_time,
                         family = betar(link = "logit"),
                         method = "REML")

# View model summary
summary(birds_FTWS_model3_GAM)

# Save model output
save(birds_FTWS_model3_GAM, 
     file = here("data", "models", "exploratory",
                 "birds_FTWS_model3_GAM.RData"))

# END OF SCRIPT ----------------------------------------------------------------