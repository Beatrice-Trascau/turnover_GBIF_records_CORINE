##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 4.3_model_setup_turnover_all_urban
# This script contains code which sets up the models exploring the imapct of 
# Forest -> TWS land cover transition on temporal turnover
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

# 2. PREPARE DATA FOR ANALYSIS -------------------------------------------------

## 2.1. Select only Forest -> TWS columns --------------------------------------

# Check column names
colnames(all_periods_turnover_all_land_cover_chanegs_15km)

# Also rename the columns for easier manipulation of df
turnover_all_urban_15km <- all_periods_turnover_all_land_cover_chanegs_15km |>
  select(-c('2000-2006_Forest no change', '2000-2006_Forest to TWS',
            '2000-2006_TWS no change', '2000-2006_TWS to Forest',
            '2006-2012_Forest no change', '2006-2012_Forest to TWS',
            '2006-2012_TWS no change', '2006-2012_TWS to Forest',
            '2012-2018_Forest no change', '2012-2018_Forest to TWS',
            '2012-2018_TWS no change', '2012-2018_TWS to Forest')) |>
  # determine which rows belong to which time period
  mutate(urban_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Urban_no_change`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_Urban_no_change`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_Urban_no_change`,
                                   TRUE ~ NA_real_),
         all_to_urban = case_when(lc_time_period == "2000-2006" ~ `2000-2006_all_to_urban`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_all_to_urban`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_all_to_urban`,
                                   TRUE ~ NA_real_)) |>
  # remove columns no longer required
  select(-`2000-2006_Urban_no_change`, -`2006-2012_Urban_no_change`, 
         -`2012-2018_Urban_no_change`,-`2000-2006_all_to_urban`, 
         -`2006-2012_all_to_urban`, -`2012-2018_all_to_urban`)

# Check transformation went ok
head(turnover_all_urban_15km)

## 2.2. Beta regression --------------------------------------------------------

# JDI is in [0,1] but beta regression requires that values do not touch 0 and 1
# so we will transform the JDI values according to this formula: 
# ( Y * (N - 1) + 0.5 ) / N, Y = response variable, N = sample size

# Get N
N <- nrow(turnover_all_urban_15km)

# Calculate new JDI values
turnover_all_urban_15km <- turnover_all_urban_15km |>
  mutate(JDI_beta = (JDI * (N - 1) + 0.5) / N)

# Check the rest of the values are what you expect them to be
glimpse(turnover_all_urban_15km) # everything looks ok

# Run model
AllUrban_turnover_model1_beta_regression <- betareg::betareg(JDI_beta ~ all_to_urban + urban_no_change + 
                               delta_recorder_effort + recorder_effort + lc_time_period,
                             data = turnover_all_urban_15km)
# Save model
save(AllUrban_turnover_model1_beta_regression,
     file = here("data", "models", "exploratory", 
                 "AllUrban_turnover_model1_beta_regression.RData"))

## 2.3. GLS with logged recorder effort ----------------------------------------

# Prepare data for GLS: remove rows with missing x or y and categorise time periods
turnover_all_urban_15km_coords_time <- turnover_all_urban_15km |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3))

# Check if there are any cells with recorder effort = 0
c <- turnover_all_urban_15km_coords_time |>
  filter(recorder_effort == 0)
length(c) #0 - Good!

# Log transform recorder effort values
turnover_all_urban_15km_coords_time <- turnover_all_urban_15km_coords_time |>
  mutate(log_recorder_effort = log(recorder_effort))

# Check log transformed values
summary(turnover_all_urban_15km_coords_time$log_recorder_effort)
any(!is.finite(turnover_all_urban_15km_coords_time$log_recorder_effort)) # FALSE = no infinite values - Good!

# Define GLS
AllUrban_turnover_model2_gls_log <- gls(JDI ~ all_to_urban + urban_no_change + 
                                      delta_recorder_effort + log_recorder_effort + lc_time_period,
                                    correlation = corExp(form = ~ x + y | time_numeric),  
                                    data = turnover_all_urban_15km_coords_time,
                                    method = "REML")

# Save model output to file
save(AllUrban_turnover_model2_gls_log, 
     file = here("data", "models", "exploratory",
                 "AllUrban_turnover_model2_gls_log.RData"))

# Model passed the validation - diagnostic plots looked good!
  # save model in the final folder as well
save(AllUrban_turnover_model2_gls_log, 
     file = here("data", "models", "final",
                 "AllUrban_turnover_model2_gls_log.RData"))

## 2.4. Plot model output ------------------------------------------------------

# Get summary of the model
model_summary <- summary(AllUrban_turnover_model2_gls_log)

# Create dataframe of coeficients
AllUrban_turnover_model2_gls_log_coef_df <- data.frame(term = names(model_summary$tTable[, "Value"]),
                                                   estimate = model_summary$tTable[, "Value"],
                                                   std.error = model_summary$tTable[, "Std.Error"],
                                                   statistic = model_summary$tTable[, "t-value"],
                                                   p.value = model_summary$tTable[, "p-value"])

# Remove the intercept
AllUrban_turnover_model2_gls_log_coef_df_no_intercept <- AllUrban_turnover_model2_gls_log_coef_df[AllUrban_turnover_model2_gls_log_coef_df$term != "(Intercept)", ]

# Create coefficient plot
figure8_a <- ggplot(AllUrban_turnover_model2_gls_log_coef_df_no_intercept, aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_y_discrete(labels = c("all_to_urban" = "TWS to Forest",
                              "urban_no_change" = "TWS No Change", 
                              "delta_recorder_effort" = "ΔRecorder Effort",
                              "log_recorder_effort" = "log Recorder Effort",
                              "lc_time_period2006-2012" = "2006-2012 Time Period",
                              "lc_time_period2012-2018" = "2012-2018 Time Period"),
                   limits = c("lc_time_period2012-2018",
                              "lc_time_period2006-2012", 
                              "log_recorder_effort",
                              "delta_recorder_effort",
                              "urban_no_change",
                              "all_to_urban")) +
  labs(x = "Estimate ± 95% CI", y = NULL) +
  theme_classic()

# Save figure as .png
ggsave(filename = here("figures", "Figure8a_AllUrban_gls_model_output_all_occurrences_turnover.png"),
       width = 12, height = 8, dpi = 300)

# Save figure as .svg
ggsave(filename = here("figures", "Figure8a_AllUrban_gls_model_output_all_occurrences_turnover.svg"),
       width = 12, height = 8, dpi = 300)

# 3. PLANT OCCURRENCES ONLY ----------------------------------------------------

## 3.1. Prepare plant data for analysis ----------------------------------------

# Select only All -> Urban columns
  # check column names
colnames(vascular_plants_all_periods_turnover_all_land_cover_chanegs_15km)

# Also rename the columns for easier manipulation of df
plants_turnover_all_urban_15km <- vascular_plants_all_periods_turnover_all_land_cover_chanegs_15km |>
  select(-c('2000-2006_Forest no change', '2000-2006_Forest to TWS',
            '2000-2006_TWS no change', '2000-2006_TWS to Forest',
            '2006-2012_Forest no change', '2006-2012_Forest to TWS',
            '2006-2012_TWS no change', '2006-2012_TWS to Forest',
            '2012-2018_Forest no change', '2012-2018_Forest to TWS',
            '2012-2018_TWS no change', '2012-2018_TWS to Forest')) |>
  # determine which rows belong to which time period
  mutate(urban_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Urban_no_change`,
                                     lc_time_period == "2006-2012" ~ `2006-2012_Urban_no_change`,
                                     lc_time_period == "2012-2018" ~ `2012-2018_Urban_no_change`,
                                     TRUE ~ NA_real_),
         all_to_urban = case_when(lc_time_period == "2000-2006" ~ `2000-2006_all_to_urban`,
                                  lc_time_period == "2006-2012" ~ `2006-2012_all_to_urban`,
                                  lc_time_period == "2012-2018" ~ `2012-2018_all_to_urban`,
                                  TRUE ~ NA_real_)) |>
  # remove columns no longer required
  select(-`2000-2006_Urban_no_change`, -`2006-2012_Urban_no_change`, 
         -`2012-2018_Urban_no_change`,-`2000-2006_all_to_urban`, 
         -`2006-2012_all_to_urban`, -`2012-2018_all_to_urban`)

# Check transformation went ok
head(plants_turnover_all_urban_15km)

# Remove rows that might have NA for x or y and categorise time periods for GLS
plants_turnover_all_urban_15km_coords_time <- plants_turnover_all_urban_15km |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort))

## 3.2. GLS with logged recorder effort ----------------------------------------

#Check if the recorder effort values were logged correctly
any(!is.finite(plants_turnover_all_urban_15km_coords_time$recorder_effort)) # FALSE = no infinite values - Good!

# Define model
plants_AllUrban_model1_gls <- gls(JDI ~ all_to_urban + urban_no_change + 
                                delta_recorder_effort + log_recorder_effort + lc_time_period,
                              correlation = corExp(form = ~ x + y | time_numeric),  
                              data = plants_turnover_all_urban_15km_coords_time,
                              method = "REML")

# Save model output
save(plants_AllUrban_model1_gls, 
     file = here("data", "models", "exploratory", "plants_AllUrban_model1_gls.RData"))

# Model diagnostics (see 4.2.qmd) revealed that model is acceptable to use 
  # saving it in the final folder as well
save(plants_AllUrban_model1_gls, 
     file = here("data", "models", "final", "plants_AllUrban_model1_gls.RData"))

## 3.4. Plot model output ------------------------------------------------------

# Get summary of the model
plants_model_summary <- summary(plants_AllUrban_model1_gls)

# Create dataframe of coeficients
plants_AllUrban_model1_glslog_coef_df <- data.frame(term = names(plants_model_summary$tTable[, "Value"]),
                                                       estimate = plants_model_summary$tTable[, "Value"],
                                                       std.error = plants_model_summary$tTable[, "Std.Error"],
                                                       statistic = plants_model_summary$tTable[, "t-value"],
                                                       p.value = plants_model_summary$tTable[, "p-value"])

# Remove the intercept
plants_AllUrban_model1_glslog_coef_df_no_intercept <- plants_AllUrban_model1_glslog_coef_df[plants_AllUrban_model1_glslog_coef_df$term != "(Intercept)", ]

# Create coefficient plot
figure8_b <- ggplot(plants_AllUrban_model1_glslog_coef_df_no_intercept, aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_y_discrete(labels = c("all_to_urban" = "TWS to Forest",
                              "urban_no_change" = "TWS No Change", 
                              "delta_recorder_effort" = "ΔRecorder Effort",
                              "log_recorder_effort" = "log Recorder Effort",
                              "lc_time_period2006-2012" = "2006-2012 Time Period",
                              "lc_time_period2012-2018" = "2012-2018 Time Period"),
                   limits = c("lc_time_period2012-2018",
                              "lc_time_period2006-2012", 
                              "log_recorder_effort",
                              "delta_recorder_effort",
                              "urban_no_change",
                              "all_to_urban")) +
  labs(x = "Estimate ± 95% CI", y = NULL) +
  theme_classic()

# Save figure as .png
ggsave(filename = here("figures", "Figure8b_AllUrban_gls_model_output_plants_turnover.png"),
       width = 12, height = 8, dpi = 300)

# Save figure as .svg
ggsave(filename = here("figures", "Figure8a_AllUrban_gls_model_output_plants_turnover.svg"),
       width = 12, height = 8, dpi = 300)

# 4. BIRD OCCURRENCES ONLY -----------------------------------------------------

## 4.1. Prepare bird data for analysis -----------------------------------------

# Check column names
colnames(birds_all_periods_turnover_all_land_cover_chanegs_15km)

# Also rename the columns for easier manipulation of df
birds_turnover_all_urban_15km <- birds_all_periods_turnover_all_land_cover_chanegs_15km |>
  select(-c('2000-2006_Forest no change', '2000-2006_Forest to TWS',
            '2000-2006_TWS no change', '2000-2006_TWS to Forest',
            '2006-2012_Forest no change', '2006-2012_Forest to TWS',
            '2006-2012_TWS no change', '2006-2012_TWS to Forest',
            '2012-2018_Forest no change', '2012-2018_Forest to TWS',
            '2012-2018_TWS no change', '2012-2018_TWS to Forest')) |>
  # determine which rows belong to which time period
  mutate(urban_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Urban_no_change`,
                                     lc_time_period == "2006-2012" ~ `2006-2012_Urban_no_change`,
                                     lc_time_period == "2012-2018" ~ `2012-2018_Urban_no_change`,
                                     TRUE ~ NA_real_),
         all_to_urban = case_when(lc_time_period == "2000-2006" ~ `2000-2006_all_to_urban`,
                                  lc_time_period == "2006-2012" ~ `2006-2012_all_to_urban`,
                                  lc_time_period == "2012-2018" ~ `2012-2018_all_to_urban`,
                                  TRUE ~ NA_real_)) |>
  # remove columns no longer required
  select(-`2000-2006_Urban_no_change`, -`2006-2012_Urban_no_change`, 
         -`2012-2018_Urban_no_change`,-`2000-2006_all_to_urban`, 
         -`2006-2012_all_to_urban`, -`2012-2018_all_to_urban`)

# Check transformation went ok
head(birds_turnover_all_urban_15km)

# Remove rows that might have NA for x or y & categorise time periods
birds_turnover_all_urban_15km_coords_time <- birds_turnover_all_urban_15km |>
  filter(!is.na(x) & !is.na(y)) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort))

## 4.2. GLS on bird data -------------------------------------------------------

# Check if the recorder effort values were logged correctly
any(!is.finite(birds_turnover_all_urban_15km_coords_time$recorder_effort)) # FALSE = no infinite values - Good!

# Define model
birds_AllUrban_model1_gls <- gls(JDI ~ all_to_urban + urban_no_change + 
                               delta_recorder_effort + log_recorder_effort + lc_time_period,
                             correlation = corExp(form = ~ x + y | time_numeric),  
                             data = birds_turnover_all_urban_15km_coords_time,
                             method = "REML")

# Save model output
save(birds_AllUrban_model1_gls, 
     file = here("data", "models", "exploratory", "birds_AllUrban_model1_gls.RData"))

# Model diagnostics (see 4.2.qmd) revealed that model is acceptable to use 
  # saving it in the final folder as well
save(birds_AllUrban_model1_gls, 
     file = here("data", "models", "final", "birds_AllUrban_model1_gls.RData"))

## 4.3. Plot model output ------------------------------------------------------

# Get summary of the model
model_summary <- summary(birds_AllUrban_model1_gls)

# Create dataframe of coeficients
birds_AllUrban_model1_gls_coef_df <- data.frame(term = names(model_summary$tTable[, "Value"]),
                                                estimate = model_summary$tTable[, "Value"],
                                                std.error = model_summary$tTable[, "Std.Error"],
                                                statistic = model_summary$tTable[, "t-value"],
                                                p.value = model_summary$tTable[, "p-value"])

# Remove the intercept
birds_AllUrban_model1_gls_coef_df_no_intercept <- birds_AllUrban_model1_gls_coef_df[birds_AllUrban_model1_gls_coef_df$term != "(Intercept)", ]

# Create coefficient plot
figure8_c <- ggplot(birds_AllUrban_model1_gls_coef_df_no_intercept, aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_y_discrete(labels = c("all_to_urban" = "All to Forest",
                              "urban_no_change" = "Forest No Change", 
                              "delta_recorder_effort" = "ΔRecorder Effort",
                              "log_recorder_effort" = "log Recorder Effort",
                              "lc_time_period2006-2012" = "2006-2012 Time Period",
                              "lc_time_period2012-2018" = "2012-2018 Time Period"),
                   limits = c("lc_time_period2012-2018",
                              "lc_time_period2006-2012", 
                              "log_recorder_effort",
                              "delta_recorder_effort",
                              "urban_no_change",
                              "all_to_urban")) +
  labs(x = "Estimate ± 95% CI", y = NULL) +
  theme_classic()

# Display plot
figure8_c

# Save figure as .png
ggsave(filename = here("figures", "Figure8c_AllUrban_gls_model_output_birds_turnover.png"),
       width = 12, height = 8, dpi = 300)

# Save figure as .svg
ggsave(filename = here("figures", "Figure8c_AllUrban_gls_model_output_birds_turnover.svg"),
       width = 12, height = 8, dpi = 300)

# END OF SCRIPT ----------------------------------------------------------------