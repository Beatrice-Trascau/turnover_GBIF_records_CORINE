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
FTWS_turnover_model4_gls_raw <- gls(beta_jtu ~ forest_to_tws + forest_no_change +
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
FTWS_turnover_model5_gls_log <- gls(beta_jtu ~ forest_to_tws + forest_no_change +
                    delta_recorder_effort + log_recorder_effort + lc_time_period + temp_change + precip_change,
                  correlation = corExp(form = ~ x + y | time_numeric),
                  data = turnover_forest_tws_15km_coords_time,
                  method = "REML")

# Model validation looks ok => THIS IS THE FINAL MODEL WE WILL USE

# Save model output to file
save(FTWS_turnover_model5_gls_log,
     file = here("data", "models", "final",
                 "FTWS_turnover_model5_gls_log.RData"))

## 2.7. GLS with logged recorder effort and interaction ------------------------

# Define GLS
FTWS_turnover_model6_gls_log_interaction <- gls(beta_jtu ~ forest_to_tws * temp_change +
                                                  forest_to_tws * precip_change + 
                                                  forest_no_change +
                                                  delta_recorder_effort + 
                                                  log_recorder_effort + 
                                                  lc_time_period,
                                                correlation = corExp(form = ~ x + y | time_numeric),
                                                data = turnover_forest_tws_15km_coords_time,
                                                method = "REML")
# Model validaiton looked ok enough

# Compare models based on AIC
AICtab(FTWS_turnover_model5_gls_log, FTWS_turnover_model6_gls_log_interaction, base = TRUE)
# model without interaction preferred => save interaction model in exploratory folder

# Save model
save(FTWS_turnover_model6_gls_log_interaction,
     file = here("data", "models", "exploratory",
                 "FTWS_turnover_model6_gls_log_interaction.RData"))

## 2.8. Get model summary ------------------------------------------------------

# Get summary of final model
model_summary <- summary(FTWS_turnover_model5_gls_log)

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
plants_FTWS_model1_gls <- gls(beta_jtu ~ forest_to_tws + forest_no_change +
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
plants_FTWS_model2_gls_interaction <- gls(beta_jtu ~ forest_to_tws * temp_change +
                                                  forest_to_tws * precip_change + 
                                                  forest_no_change +
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
birds_FTWS_model1_gls <- gls(beta_jtu ~ forest_to_tws + forest_no_change +
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
birds_FTWS_model2_gls_interaction <- gls(beta_jtu ~ forest_to_tws * temp_change +
                                            forest_to_tws * precip_change + 
                                            forest_no_change +
                                            delta_recorder_effort + 
                                            log_recorder_effort + 
                                            lc_time_period,
                                          correlation = corExp(form = ~ x + y | time_numeric),
                                          data = plants_turnover_forest_tws_15km_coords_time,
                                          method = "REML")
# Model validaiton looked ok enough

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

# Function to create coefficient plot for a single model
create_coef_plot <- function(model, title, show_y_axis = TRUE) {
  
  # get model summary
  model_summary <- summary(model)
  
  # create a dataframe of the coefficients
  coef_df <- data.frame(term = names(model_summary$tTable[, "Value"]),
                        estimate = model_summary$tTable[, "Value"],
                        std.error = model_summary$tTable[, "Std.Error"],
                        statistic = model_summary$tTable[, "t-value"],
                        p.value = model_summary$tTable[, "p-value"]) |>
    mutate(conf.low = estimate - 1.96 * std.error,
           conf.high = estimate + 1.96 * std.error,
           significance = case_when(p.value < 0.001 ~ "***",
                                    p.value < 0.01 ~ "**", 
                                    p.value < 0.05 ~ "*",
                                    p.value < 0.1 ~ ".",
                                    TRUE ~ ""),
      effect_type = case_when(term == "(Intercept)" ~ "Intercept",
                              p.value < 0.05 & estimate < 0 ~ "Negative (sig.)",
                              p.value < 0.05 & estimate > 0 ~ "Positive (sig.)",
                              TRUE ~ "Non-significant"),
      term_clean = case_when(term == "(Intercept)" ~ "Intercept",
                             term == "forest_to_tws" ~ "Forest → TWS",
                             term == "forest_no_change" ~ "Forest (no change)",
                             term == "delta_recorder_effort" ~ "ΔRecorder effort",
                             term == "log_recorder_effort" ~ "log(Recorder effort)",
                             term == "lc_time_period2006-2012" ~ "Period: 2006-2012",
                             term == "lc_time_period2012-2018" ~ "Period: 2012-2018",
                             term == "temp_change" ~ "Temperature change",
                             term == "precip_change" ~ "Precipitation change",
                             TRUE ~ term))
  
  # split into intercept and effects
  intercept_df <- coef_df |> filter(term == "(Intercept)")
  effects_df <- coef_df |> 
    filter(term != "(Intercept)") |>
    mutate(term_clean = factor(term_clean, 
                               levels = rev(c("Forest → TWS",
                                              "Forest (no change)", 
                                              "Temperature change",
                                              "Precipitation change",
                                              "ΔRecorder effort",
                                              "log(Recorder effort)",
                                              "Period: 2006-2012",
                                              "Period: 2012-2018"))))
  
  # define colors
  effect_colors <- c( "Negative (sig.)" = "#9C27B0",
                      "Positive (sig.)" = "#FF9800",  
                      "Non-significant" = "grey60",
                      "Intercept" = "grey30")
  
  # create intercept plot (top)
  plot_intercept <- ggplot(intercept_df, aes(x = estimate, y = term_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                   height = 0.2, linewidth = 0.7, color = "grey30") +
    geom_point(aes(color = effect_type), size = 3, shape = 15) +
    geom_text(aes(x = conf.high, label = significance), 
              hjust = -0.3, vjust = 0.5, size = 3.5, fontface = "bold") +
    scale_color_manual(values = effect_colors, guide = "none") +
    labs(x = NULL, y = NULL, title = title) +
    theme_classic() +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 10, color = "black", face = "bold"),
          axis.text.x = element_text(size = 9, color = "black"),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
          plot.margin = margin(t = 5, r = 10, b = 2, l = 5)) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)))
  
  # create effects plot (bottom)
  plot_effects <- ggplot(effects_df, aes(x = estimate, y = term_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                   height = 0.2, linewidth = 0.7, color = "grey30") +
    geom_point(aes(color = effect_type), size = 3) +
    geom_text(aes(x = conf.high, label = significance), 
              hjust = -0.3, vjust = 0.5, size = 3.5, fontface = "bold") +
    scale_color_manual(values = effect_colors,
                       name = "Effect type",
                       breaks = c("Positive (sig.)", "Negative (sig.)", "Non-significant"),
                       labels = c("Positive (p < 0.05)", "Negative (p < 0.05)", "Non-significant")) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 9, color = "black"),
          legend.position = "none",
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
          plot.margin = margin(t = 2, r = 10, b = 5, l = 5)) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)))
  
  # show or hide y-axis labels
  if (!show_y_axis) {
    plot_intercept <- plot_intercept + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_effects <- plot_effects + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  # combine intercept and effects for this model
  combined <- plot_intercept / plot_effects + 
    plot_layout(heights = c(1, 4))
  
  return(combined)
}

# Create plots for each model
plot_all <- create_coef_plot(FTWS_turnover_model5_gls_log, 
                             "a) All occurrences", 
                             show_y_axis = TRUE)

plot_plants <- create_coef_plot(plants_FTWS_model1_gls, 
                                "b) Vascular plants", 
                                show_y_axis = FALSE)

plot_birds <- create_coef_plot(birds_FTWS_model1_gls, 
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

ggsave(filename = here("figures", "Figure_combined_models_3panel.pdf"),
       plot = final_plot, width = 14, height = 7, bg = "white")

# END OF SCRIPT ----------------------------------------------------------------