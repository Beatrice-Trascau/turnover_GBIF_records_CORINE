##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 5.2_model_setup_turnover_tws_forest
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

## 2.1. Select only Forest -> TWS columns --------------------------------------

# Check column names
colnames(all_periods_turnover_with_climate)

# Total number of 100m x 100m pixels in a 15km x 15km cell
# 15,000m / 100m = 150 pixels per side
# 150 x 150 = 22,500 total pixels per cell
total_pixels_per_cell <- 22500

# Also rename the columns for easier manipulation of df
turnover_tws_forest_15km <- all_periods_turnover_with_climate |>
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
  # convert land-cover changes to proportions
  mutate(tws_no_change_prop = tws_no_change / total_pixels_per_cell,
         tws_to_forest_prop = tws_to_forest / total_pixels_per_cell) |>
  # remove columns no longer required
  select(-`2000-2006_TWS no change`, -`2006-2012_TWS no change`, 
         -`2012-2018_TWS no change`,-`2000-2006_TWS to Forest`, 
         -`2006-2012_TWS to Forest`, -`2012-2018_TWS to Forest`)

# Check transformation went ok
head(turnover_tws_forest_15km)

## 2.2. GLS with logged recorder effort ----------------------------------------

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
TWSF_turnover_model1 <- gls(beta_jtu ~ tws_to_forest_prop + tws_no_change_prop + 
                                      delta_recorder_effort + log_recorder_effort + lc_time_period + temp_change + precip_change,
                                    correlation = corExp(form = ~ x + y | time_numeric),  
                                    data = turnover_tws_forest_15km_coords_time,
                                    method = "REML")

# Model passed the validation - diagnostic plots looked good!
  # save model in the final folder as well
save(TWSF_turnover_model1, 
     file = here("data", "models", "final",
                 "TWSF_turnover_model1.RData"))

## 2.3. GLS with interaction ---------------------------------------------------

# Define GLS
TWSF_turnover_model2_interaction <- gls(beta_jtu ~ tws_to_forest_prop * temp_change +
                                          tws_to_forest_prop * precip_change + 
                                          tws_no_change_prop +
                                                  delta_recorder_effort + 
                                                  log_recorder_effort + 
                                                  lc_time_period,
                                                correlation = corExp(form = ~ x + y | time_numeric),
                                                data = turnover_tws_forest_15km_coords_time,
                                                method = "REML")
# Model validaiton looked ok enough

# Compare models based on AIC
AICtab(TWSF_turnover_model1, TWSF_turnover_model2_interaction, base = TRUE)
# model without interaction preferred => save interaction model in exploratory folder

# Save model
save(TWSF_turnover_model2_interaction,
     file = here("data", "models", "exploratory",
                 "TWSF_turnover_model2_interaction.RData"))

## 2.3. Get model summary ------------------------------------------------------

# Get summary of final model
model_summary <- summary(TWSF_turnover_model1)

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

# Select only TWS -> Forest columns
# check column names
colnames(vascular_plants_turnover_with_climate)

# Also rename the columns for easier manipulation of df
plants_turnover_tws_forest_15km <- vascular_plants_turnover_with_climate |>
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
  # convert land-cover changes to proportions
  mutate(tws_no_change_prop = tws_no_change / total_pixels_per_cell,
         tws_to_forest_prop = tws_to_forest / total_pixels_per_cell) |>
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
plants_TWSF_model1_gls <- gls(beta_jtu ~ tws_to_forest_prop + tws_no_change_prop + 
                                delta_recorder_effort + log_recorder_effort + lc_time_period + temp_change + precip_change,
                              correlation = corExp(form = ~ x + y | time_numeric),  
                              data = plants_turnover_tws_forest_15km_coords_time,
                              method = "REML")


# Model diagnostics (see 4.2.qmd) revealed that model is acceptable to use 
  # saving it in the final folder as well
save(plants_TWSF_model1_gls, 
     file = here("data", "models", "final", "plants_TWSF_model1_gls.RData"))

## 3.3. Plant GLS with interaction ---------------------------------------------

# Define model
plants_TWSF_model2_gls_interaction <- gls(beta_jtu ~ tws_to_forest_prop * temp_change +
                                            tws_to_forest_prop * precip_change + 
                                            tws_no_change_prop +
                                            delta_recorder_effort + 
                                            log_recorder_effort + 
                                            lc_time_period,
                                          correlation = corExp(form = ~ x + y | time_numeric),
                                          data = plants_turnover_tws_forest_15km_coords_time,
                                          method = "REML")

# Model validaiton looked ok enough

# Compare models based on AIC
AICtab(plants_TWSF_model1_gls, plants_TWSF_model2_gls_interaction, base = TRUE)
# model without interaction preferred => save interaction model in exploratory folder

# Save model
save(plants_TWSF_model2_gls_interaction,
     file = here("data", "models", "exploratory",
                 "plants_TWSF_model2_gls_interaction.RData"))

## 3.4. Get plant model summary ------------------------------------------------

# Get summary of final model
model_summary <- summary(plants_TWSF_model1_gls)

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
birds_turnover_tws_forest_15km <- birds_turnover_with_climate |>
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
  # convert land-cover changes to proportions
  mutate(tws_no_change_prop = tws_no_change / total_pixels_per_cell,
         tws_to_forest_prop = tws_to_forest / total_pixels_per_cell) |>
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
birds_TWSF_model1_gls <- gls(beta_jtu ~ tws_to_forest_prop + tws_no_change_prop + 
                                delta_recorder_effort + log_recorder_effort + lc_time_period + temp_change + precip_change,
                              correlation = corExp(form = ~ x + y | time_numeric),  
                              data = birds_turnover_tws_forest_15km_coords_time,
                              method = "REML")

# Model diagnostics (see 4.2.qmd) revealed that model is acceptable to use 
  # saving it in the final folder as well
save(birds_TWSF_model1_gls, 
     file = here("data", "models", "final", "birds_TWSF_model1_gls.RData"))

## 4.3. Bird GLS with interaction ----------------------------------------------

# Define model
birds_TWSF_model2_gls_interaction <- gls(beta_jtu ~ tws_to_forest_prop * temp_change +
                                           tws_to_forest_prop * precip_change + 
                                           tws_no_change_prop +
                                           delta_recorder_effort + 
                                           log_recorder_effort + 
                                           lc_time_period,
                                         correlation = corExp(form = ~ x + y | time_numeric),
                                         data = birds_turnover_tws_forest_15km_coords_time,
                                         method = "REML")

# Model validation looked ok enough

# Compare models based on AIC
AICtab(birds_TWSF_model1_gls, birds_TWSF_model2_gls_interaction, base = TRUE)
# model without interaction preferred => save interaction model in exploratory folder

# Save model
save(birds_TWSF_model2_gls_interaction,
     file = here("data", "models", "exploratory",
                 "birds_TWSF_model2_gls_interaction.RData"))

## 4.4. Get model summary ------------------------------------------------------

# Get summary of final model
model_summary <- summary(birds_TWSF_model1_gls)

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

# Use function defined in section 14.2 of the 0_setup.R script

# Create plots for each model
plot_all <- create_coef_plot_tws(TWSF_turnover_model1, 
                             "a) All occurrences", 
                             show_y_axis = TRUE)

plot_plants <- create_coef_plot_tws(plants_TWSF_model1_gls, 
                                "b) Vascular plants", 
                                show_y_axis = FALSE)

plot_birds <- create_coef_plot_tws(birds_TWSF_model1_gls, 
                               "c) Birds", 
                               show_y_axis = FALSE)

# Combine all three panels horizontally
combined_plot <- plot_all | plot_plants | plot_birds

# Create a dummy plot just to extract the legend
dummy_effects <- data.frame(term_clean = "dummy",estimate = 0,
                            conf.low = -0.1, conf.high = 0.1,
                            effect_type = c("Positive (sig.)", "Negative (sig.)", "Non-significant"))

# Define colours to tell effects apart
effect_colors <- c("Negative (sig.)" = "#9C27B0",
                   "Positive (sig.)" = "#FF9800",
                   "Non-significant" = "grey60")

# Create legend
legend_plot <- ggplot(dummy_effects, aes(x = estimate, y = term_clean, color = effect_type)) +
  geom_point() +
  scale_color_manual(values = effect_colors,
                     name = NULL,
                     breaks = c("Positive (sig.)", "Negative (sig.)", "Non-significant"),
                     labels = c("Positive effect (p < 0.05)", "Negative effect (p < 0.05)", "Non-significant")) +
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
ggsave(filename = here("figures", "Figure6def_betaJTU_TWSF_models.png"),
       plot = final_plot, width = 14, height = 7, dpi = 300, bg = "white")

ggsave(filename = here("figures", "Figure6def_betaJTU_TWSF_models.pdf"),
       plot = final_plot, width = 14, height = 7, bg = "white")

# 6. PREDICTION PLOTS ----------------------------------------------------------

## 6.1. All occurrences --------------------------------------------------------

# Get the data and model coefficients
data_clean <- turnover_tws_forest_15km_coords_time
coefs <- coef(TWSF_turnover_model1)

# Define colors
line_color <- "#2196F3"

# Plot 1: tws_to_forest 

# Get range of values
pred_range_1 <- seq(min(data_clean$tws_to_forest, na.rm = TRUE),
                    max(data_clean$tws_to_forest, na.rm = TRUE),
                    length.out = 100)

# Manual prediction calculation
pred_1 <- coefs["(Intercept)"] + 
  coefs["tws_to_forest"] * pred_range_1 +
  coefs["tws_no_change"] * mean(data_clean$tws_no_change, na.rm = TRUE) +
  coefs["delta_recorder_effort"] * mean(data_clean$delta_recorder_effort, na.rm = TRUE) +
  coefs["log_recorder_effort"] * mean(data_clean$log_recorder_effort, na.rm = TRUE) +
  coefs["temp_change"] * mean(data_clean$temp_change, na.rm = TRUE) +
  coefs["precip_change"] * mean(data_clean$precip_change, na.rm = TRUE) +
  coefs["lc_time_period2006-2012"] * 1  # Reference level = 2006-2012

# Plot prediction for forest to tws
plot1 <- ggplot(data.frame(x = pred_range_1, y = pred_1), aes(x = x, y = y)) +
  geom_line(color = line_color, linewidth = 1.2) +
  labs(x = "TWS → Forest (proportion)", y = "Predicted beta_jtu", title = "a)") +
  theme_classic() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3))

# Plot 2: tws_no_change 

# Get range of values for tws_no_change 
pred_range_2 <- seq(min(data_clean$tws_no_change, na.rm = TRUE),
                    max(data_clean$tws_no_change, na.rm = TRUE),
                    length.out = 100)

# Manual calculation of predictors
pred_2 <- coefs["(Intercept)"] + 
  coefs["tws_to_forest"] * mean(data_clean$tws_to_forest, na.rm = TRUE) +
  coefs["tws_no_change"] * pred_range_2 +
  coefs["delta_recorder_effort"] * mean(data_clean$delta_recorder_effort, na.rm = TRUE) +
  coefs["log_recorder_effort"] * mean(data_clean$log_recorder_effort, na.rm = TRUE) +
  coefs["temp_change"] * mean(data_clean$temp_change, na.rm = TRUE) +
  coefs["precip_change"] * mean(data_clean$precip_change, na.rm = TRUE) +
  coefs["lc_time_period2006-2012"] * 1

# Plot prediction for forest no change
plot2 <- ggplot(data.frame(x = pred_range_2, y = pred_2), aes(x = x, y = y)) +
  geom_line(color = line_color, linewidth = 1.2) +
  labs(x = "TWS no change (proportion)", y = NULL, title = "b)") +
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
  coefs["tws_to_forest"] * mean(data_clean$tws_to_forest, na.rm = TRUE) +
  coefs["tws_no_change"] * mean(data_clean$tws_no_change, na.rm = TRUE) +
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
  coefs["tws_to_forest"] * mean(data_clean$tws_to_forest, na.rm = TRUE) +
  coefs["tws_no_change"] * mean(data_clean$tws_no_change, na.rm = TRUE) +
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
  coefs["tws_to_forest"] * mean(data_clean$tws_to_forest, na.rm = TRUE) +
  coefs["tws_no_change"] * mean(data_clean$tws_no_change, na.rm = TRUE) +
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
  coefs["tws_to_forest"] * mean(data_clean$tws_to_forest, na.rm = TRUE) +
  coefs["tws_no_change"] * mean(data_clean$tws_no_change, na.rm = TRUE) +
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

# END OF SCRIPT ----------------------------------------------------------------