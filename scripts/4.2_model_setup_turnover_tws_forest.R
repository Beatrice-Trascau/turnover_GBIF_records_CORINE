##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 4.2_model_setup_turnover_tws_forest
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

## 2.1. Select only Forest -> TWS columns --------------------------------------

# Check column names
colnames(all_periods_turnover_all_land_cover_chanegs_15km)

# Also rename the columns for easier manipulation of df
turnover_tws_forest_15km <- all_periods_turnover_all_land_cover_chanegs_15km |>
  select(-c('2000-2006_Forest no change', '2000-2006_Forest to TWS',
            '2000-2006_Urban_no_change', '2000-2006_all_to_urban',
            '2006-2012_Forest no change', '2006-2012_Forest to TWS',
            '2006-2012_Urban_no_change', '2006-2012_all_to_urban',
            '2012-2018_Forest no change', '2012-2018_Forest to TWS',
            '2012-2018_Urban_no_change', '2012-2018_all_to_urban')) |>
  # determine which rows belong to which time period
  mutate(tws_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS no change`,
                                      lc_time_period == "2006-2012" ~ `2006-2012_TWS no change`,
                                      lc_time_period == "2012-2018" ~ `2012-2018_TWS no change`,
                                      TRUE ~ NA_real_),
         tws_to_forest = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS to Forest`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_TWS to Forest`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_TWS to Forest`,
                                   TRUE ~ NA_real_)) |>
  # remove columns no longer required
  select(-`2000-2006_TWS no change`, -`2006-2012_TWS no change`, 
         -`2012-2018_TWS no change`,-`2000-2006_TWS to Forest`, 
         -`2006-2012_TWS to Forest`, -`2012-2018_TWS to Forest`)

# Check transformation went ok
head(turnover_tws_forest_15km)

## 2.2. Beta regression --------------------------------------------------------

# JDI is in [0,1] but beta regression requires that values do not touch 0 and 1
# so we will transform the JDI values according to this formula: 
# ( Y * (N - 1) + 0.5 ) / N, Y = response variable, N = sample size

# Get N
N <- nrow(turnover_tws_forest_15km)

# Calculate new JDI values
turnover_tws_forest_15km <- turnover_tws_forest_15km |>
  mutate(JDI_beta = (JDI * (N - 1) + 0.5) / N)

# Check the rest of the values are what you expect them to be
glimpse(turnover_tws_forest_15km) # everything looks ok

# Run model
TWSF_turnover_model1_beta_regression <- betareg::betareg(JDI_beta ~ tws_to_forest + tws_no_change + 
                               delta_recorder_effort + recorder_effort + lc_time_period,
                             data = turnover_tws_forest_15km)

# Save model output
save(TWSF_turnover_model1_beta_regression,
     file = here("data", "models", "exploratory", 
                 "TWSF_turnover_model1_beta_regression.RData"))

## 2.3. GLS with logged recorder effort ----------------------------------------

# Prepare data for GLS: remove rows with missing x or y and categorise time periods
turnover_tws_forest_15km_coords_time <- turnover_tws_forest_15km |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3))

# Check if there are any cells with recorder effort = 0
b <- turnover_tws_forest_15km_coords_time |>
  filter(recorder_effort == 0)
length(b) #0 - Good!

# Log transform recorder effort values
turnover_tws_forest_15km_coords_time <- turnover_tws_forest_15km_coords_time |>
  mutate(log_recorder_effort = log(recorder_effort))

# Check log transformed values
summary(turnover_tws_forest_15km_coords_time$log_recorder_effort)
any(!is.finite(turnover_tws_forest_15km_coords_time$log_recorder_effort)) # FALSE = no infinite values - Good!

# Define GLS
TWSF_turnover_model2_gls_log <- gls(JDI ~ tws_to_forest + tws_no_change + 
                                      delta_recorder_effort + log_recorder_effort + lc_time_period,
                                    correlation = corExp(form = ~ x + y | time_numeric),  
                                    data = turnover_tws_forest_15km_coords_time,
                                    method = "REML")

# Save model output to file
save(TWSF_turnover_model2_gls_log, 
     file = here("data", "models", "exploratory",
                 "TWSF_turnover_model2_gls_log.RData"))

# Model passed the validation - diagnostic plots looked good!
  # save model in the final folder as well
save(TWSF_turnover_model2_gls_log, 
     file = here("data", "models", "final",
                 "TWSF_turnover_model2_gls_log.RData"))

## 2.4. Plot model output ------------------------------------------------------

# Get summary of the model
model_summary <- summary(TWSF_turnover_model2_gls_log)

# Create dataframe of coeficients
TWSF_turnover_model2_gls_log_coef_df <- data.frame(term = names(model_summary$tTable[, "Value"]),
                                                   estimate = model_summary$tTable[, "Value"],
                                                   std.error = model_summary$tTable[, "Std.Error"],
                                                   statistic = model_summary$tTable[, "t-value"],
                                                   p.value = model_summary$tTable[, "p-value"])

# Remove the intercept
TWSF_turnover_model2_gls_log_coef_df_no_intercept <- TWSF_turnover_model2_gls_log_coef_df[TWSF_turnover_model2_gls_log_coef_df$term != "(Intercept)", ]

# Create coefficient plot
figure7_a <- ggplot(TWSF_turnover_model2_gls_log_coef_df_no_intercept, aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_y_discrete(labels = c("tws_to_forest" = "TWS to Forest",
                              "tws_no_change" = "TWS No Change", 
                              "delta_recorder_effort" = "ΔRecorder Effort",
                              "log_recorder_effort" = "log Recorder Effort",
                              "lc_time_period2006-2012" = "2006-2012 Time Period",
                              "lc_time_period2012-2018" = "2012-2018 Time Period"),
                   limits = c("lc_time_period2012-2018",
                              "lc_time_period2006-2012", 
                              "log_recorder_effort",
                              "delta_recorder_effort",
                              "tws_no_change",
                              "tws_to_forest")) +
  labs(x = "Estimate ± 95% CI", y = NULL) +
  theme_classic()

# Save figure as .png
ggsave(filename = here("figures", "Figure7a_TWSF_gls_model_output_all_occurrences_turnover.png"),
       width = 12, height = 8, dpi = 300)

# Save figure as .svg
ggsave(filename = here("figures", "Figure7a_TWSF_gls_model_output_all_occurrences_turnover.svg"),
       width = 12, height = 8, dpi = 300)

# 3. PLANT OCCURRENCES ONLY ----------------------------------------------------

## 3.1. Prepare plant data for analysis ----------------------------------------

# Select only TWS -> Forest columns
# check column names
colnames(vascular_plants_all_periods_turnover_all_land_cover_chanegs_15km)

# Also rename the columns for easier manipulation of df
plants_turnover_tws_forest_15km <- vascular_plants_all_periods_turnover_all_land_cover_chanegs_15km |>
  select(-c('2000-2006_Forest no change', '2000-2006_Forest to TWS',
            '2000-2006_Urban_no_change', '2000-2006_all_to_urban',
            '2006-2012_Forest no change', '2006-2012_Forest to TWS',
            '2006-2012_Urban_no_change', '2006-2012_all_to_urban',
            '2012-2018_Forest no change', '2012-2018_Forest to TWS',
            '2012-2018_Urban_no_change', '2012-2018_all_to_urban')) |>
  # determine which rows belong to which time period
  mutate(tws_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS no change`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_TWS no change`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_TWS no change`,
                                   TRUE ~ NA_real_),
         tws_to_forest = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS to Forest`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_TWS to Forest`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_TWS to Forest`,
                                   TRUE ~ NA_real_)) |>
  # remove columns no longer required
  select(-`2000-2006_TWS no change`, -`2006-2012_TWS no change`, 
         -`2012-2018_TWS no change`,-`2000-2006_TWS to Forest`, 
         -`2006-2012_TWS to Forest`, -`2012-2018_TWS to Forest`)

# Check transformation went ok
head(plants_turnover_tws_forest_15km)

# Remove rows that might have NA for x or y and categorise time periods for GLS
plants_turnover_tws_forest_15km_coords_time <- plants_turnover_tws_forest_15km |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort))

## 3.2. GLS on Plant data ------------------------------------------------------

# Check if the recorder effort values were logged correctly
any(!is.finite(plants_turnover_tws_forest_15km_coords_time$recorder_effort)) # FALSE = no infinite values - Good!

# Define model
plants_TWSF_model1_gls <- gls(JDI ~ tws_to_forest + tws_no_change + 
                                delta_recorder_effort + log_recorder_effort + lc_time_period,
                              correlation = corExp(form = ~ x + y | time_numeric),  
                              data = plants_turnover_tws_forest_15km_coords_time,
                              method = "REML")

# Save model output
save(plants_TWSF_model1_gls, 
     file = here("data", "models", "exploratory", "plants_TWSF_model1_gls.RData"))

# Model diagnostics (see 4.2.qmd) revealed that model is acceptable to use 
  # saving it in the final folder as well
save(plants_TWSF_model1_gls, 
     file = here("data", "models", "final", "plants_TWSF_model1_gls.RData"))

## 3.3. Plot model output ------------------------------------------------------

# Get model summary
plants_TWSF_model1_gls_summary <- summary(plants_TWSF_model1_gls)

# Create dataframe of coeficients
plants_TWSF_model1_gls_summary_coef_df <- data.frame(term = names(plants_TWSF_model1_gls_summary$tTable[, "Value"]),
                                             estimate = plants_TWSF_model1_gls_summary$tTable[, "Value"],
                                             std.error = plants_TWSF_model1_gls_summary$tTable[, "Std.Error"],
                                             statistic = plants_TWSF_model1_gls_summary$tTable[, "t-value"],
                                             p.value = plants_TWSF_model1_gls_summary$tTable[, "p-value"])

# Remove the intercept
plants_TWSF_model1_gls_summary_coef_df_no_intercept <- plants_TWSF_model1_gls_summary_coef_df[plants_TWSF_model1_gls_summary_coef_df$term != "(Intercept)", ]

# Create coefficient plot
figure7_b <- ggplot(plants_TWSF_model1_gls_coef_df_no_intercept, aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_y_discrete(labels = c("tws_to_forest" = "TWS to Forest",
                              "tws_no_change" = "TWS No Change", 
                              "delta_recorder_effort" = "ΔRecorder Effort",
                              "log_recorder_effort" = "log Recorder Effort",
                              "lc_time_period2006-2012" = "2006-2012 Time Period",
                              "lc_time_period2012-2018" = "2012-2018 Time Period"),
                   limits = c("lc_time_period2012-2018",
                              "lc_time_period2006-2012", 
                              "log_recorder_effort",
                              "delta_recorder_effort",
                              "tws_no_change",
                              "tws_to_forest")) +
  labs(x = "Estimate ± 95% CI", y = NULL) +
  theme_classic()

# Save figure as .png
ggsave(filename = here("figures", "Figure7b_TWSF_gls_model_output_plants_turnover.png"),
       width = 12, height = 8, dpi = 300)

# Save figure as .svg
ggsave(filename = here("figures", "Figure7b_TWS_gls_model_output_plants_turnover.svg"),
       width = 12, height = 8, dpi = 300)

# 4. BIRD OCCURRENCES ONLY -----------------------------------------------------

## 4.1. Prepare bird data for analysis -----------------------------------------

# Check column names
colnames(birds_all_periods_turnover_all_land_cover_chanegs_15km)

# Also rename the columns for easier manipulation of df
birds_turnover_tws_forest_15km <- birds_all_periods_turnover_all_land_cover_chanegs_15km |>
  select(-c('2000-2006_Forest no change', '2000-2006_Forest to TWS',
            '2000-2006_Urban_no_change', '2000-2006_all_to_urban',
            '2006-2012_Forest no change', '2006-2012_Forest to TWS',
            '2006-2012_Urban_no_change', '2006-2012_all_to_urban',
            '2012-2018_Forest no change', '2012-2018_Forest to TWS',
            '2012-2018_Urban_no_change', '2012-2018_all_to_urban')) |>
  # determine which rows belong to which time period
  mutate(tws_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS no change`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_TWS no change`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_TWS no change`,
                                   TRUE ~ NA_real_),
         tws_to_forest = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS to Forest`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_TWS to Forest`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_TWS to Forest`,
                                   TRUE ~ NA_real_)) |>
  # remove columns no longer required
  select(-`2000-2006_TWS no change`, -`2006-2012_TWS no change`, 
         -`2012-2018_TWS no change`,-`2000-2006_TWS to Forest`, 
         -`2006-2012_TWS to Forest`, -`2012-2018_TWS to Forest`)

# Check transformation went ok
head(birds_turnover_tws_forest_15km)

# Remove rows that might have NA for x or y & categorise time periods
birds_turnover_tws_forest_15km_coords_time <- birds_turnover_tws_forest_15km |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort))

## 4.2. GLS on bird data -------------------------------------------------------

# Check if the recorder effort values were logged correctly
any(!is.finite(birds_turnover_tws_forest_15km_coords_time$recorder_effort)) # FALSE = no infinite values - Good!

# Define model
birds_TWSF_model1_gls <- gls(JDI ~ tws_to_forest + tws_no_change + 
                                delta_recorder_effort + log_recorder_effort + lc_time_period,
                              correlation = corExp(form = ~ x + y | time_numeric),  
                              data = birds_turnover_tws_forest_15km_coords_time,
                              method = "REML")

# Save model output
save(birds_TWSF_model1_gls, 
     file = here("data", "models", "exploratory", "birds_TWSF_model1_gls.RData"))

# Model diagnostics (see 4.2.qmd) revealed that model is acceptable to use 
  # saving it in the final folder as well
save(birds_TWSF_model1_gls, 
     file = here("data", "models", "final", "birds_TWSF_model1_gls.RData"))

## 4.3. Plot model output ------------------------------------------------------

# Get summary of the model
birds_summary <- summary(birds_TWSF_model1_gls)

# Create dataframe of coeficients
birds_TWSF_model1_gls_coef_df <- data.frame(term = names(birds_summary$tTable[, "Value"]),
                                            estimate = birds_summary$tTable[, "Value"],
                                            std.error = birds_summary$tTable[, "Std.Error"],
                                            statistic = birds_summary$tTable[, "t-value"],
                                            p.value = birds_summary$tTable[, "p-value"])

# Remove the intercept
birds_TWSF_model1_gls_coef_df_no_intercept <- birds_TWSF_model1_gls_coef_df[birds_TWSF_model1_gls_coef_df$term != "(Intercept)", ]

# Create coefficient plot
figure7_c <- ggplot(birds_TWSF_model1_gls_coef_df_no_intercept, aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_y_discrete(labels = c("tws_to_forest" = "TWS to Forest",
                              "tws_no_change" = "TWS No Change", 
                              "delta_recorder_effort" = "ΔRecorder Effort",
                              "log_recorder_effort" = "log Recorder Effort",
                              "lc_time_period2006-2012" = "2006-2012 Time Period",
                              "lc_time_period2012-2018" = "2012-2018 Time Period"),
                   limits = c("lc_time_period2012-2018",
                              "lc_time_period2006-2012", 
                              "log_recorder_effort",
                              "delta_recorder_effort",
                              "tws_no_change",
                              "tws_to_forest")) +
  labs(x = "Estimate ± 95% CI", y = NULL) +
  theme_classic()

# Display plot
figure7_c

# Save figure as .png
ggsave(filename = here("figures", "Figure7c_TWSF_gls_model_output_birds_turnover.png"),
       width = 12, height = 8, dpi = 300)

# Save figure as .svg
ggsave(filename = here("figures", "Figure7c_TWSF_gls_model_output_birds_turnover.svg"),
       width = 12, height = 8, dpi = 300)

# END OF SCRIPT ----------------------------------------------------------------