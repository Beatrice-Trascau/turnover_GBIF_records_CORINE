##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 4.3_model_setup_turnover_all_urban
# This script contains code which sets up the models exploring the imapct of 
# Forest -> TWS land cover transition on temporal turnover
##----------------------------------------------------------------------------##

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


# 2. PREPARE DATA FOR ANALYSIS -------------------------------------------------

## 2.1. Select only Forest -> TWS columns --------------------------------------

# Check column names
colnames(all_periods_turnover_with_climate)

# Also rename the columns for easier manipulation of df
turnover_all_urban_15km <- all_periods_turnover_with_climate |>
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

## 2.2. GLS with logged recorder effort ----------------------------------------

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
AllUrban_turnover_model1_gls <- gls(beta_jtu ~ all_to_urban + urban_no_change + 
                                      delta_recorder_effort + log_recorder_effort + lc_time_period + temp_change + precip_change,
                                    correlation = corExp(form = ~ x + y | time_numeric),  
                                    data = turnover_all_urban_15km_coords_time,
                                    method = "REML")

# Model passed the validation - diagnostic plots looked good!
  # save model in the final folder as well
save(AllUrban_turnover_model1_gls, 
     file = here("data", "models", "final",
                 "AllUrban_turnover_model1_gls.RData"))

## 2.3. GLS with interaction term ----------------------------------------------

# Define GLS
AllUrban_turnover_model2_gls_interaction <- gls(beta_jtu ~ all_to_urban * temp_change +
                                                  all_to_urban * precip_change + 
                                                  urban_no_change +
                                                  delta_recorder_effort + 
                                                  log_recorder_effort + 
                                                  lc_time_period,
                                                correlation = corExp(form = ~ x + y | time_numeric),
                                                data = turnover_all_urban_15km_coords_time,
                                                method = "REML")
# Model validaiton looked ok enough

# Compare models based on AIC
AICtab(AllUrban_turnover_model1_gls, AllUrban_turnover_model2_gls_interaction, base = TRUE)
# model without interaction preferred => save interaction model in exploratory folder

# Save model
save(AllUrban_turnover_model2_gls_interaction,
     file = here("data", "models", "exploratory",
                 "AllUrban_turnover_model2_gls_interaction.RData"))

## 2.4. Get model summary ------------------------------------------------------

# Get summary of final model
model_summary <- summary(AllUrban_turnover_model1_gls)

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

# Select only All -> Urban columns
  # check column names
colnames(vascular_plants_turnover_with_climate)

# Also rename the columns for easier manipulation of df
plants_turnover_all_urban_15km <- vascular_plants_turnover_with_climate |>
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
plants_AllUrban_model1_gls <- gls(beta_jtu ~ all_to_urban + urban_no_change + 
                                delta_recorder_effort + log_recorder_effort + lc_time_period + temp_change + precip_change,
                              correlation = corExp(form = ~ x + y | time_numeric),  
                              data = plants_turnover_all_urban_15km_coords_time,
                              method = "REML")


# Model diagnostics (see 4.2.qmd) revealed that model is acceptable to use 
  # saving it in the final folder as well
save(plants_AllUrban_model1_gls, 
     file = here("data", "models", "final", "plants_AllUrban_model1_gls.RData"))

## 3.3. Plant GLS with interaction ---------------------------------------------

# Define GLS
plants_AllUrban_model2_gls_interaction <- gls(beta_jtu ~ all_to_urban * temp_change +
                                                all_to_urban * precip_change + 
                                                urban_no_change +
                                            delta_recorder_effort + 
                                            log_recorder_effort + 
                                            lc_time_period,
                                          correlation = corExp(form = ~ x + y | time_numeric),
                                          data = plants_turnover_all_urban_15km_coords_time,
                                          method = "REML")
# Model validaiton looked ok enough

# Compare models based on AIC
AICtab(plants_AllUrban_model1_gls, plants_AllUrban_model2_gls_interaction, base = TRUE)
# model without interaction preferred => save interaction model in exploratory folder

# Save model
save(plants_AllUrban_model2_gls_interaction,
     file = here("data", "models", "exploratory",
                 "plants_AllUrban_model2_gls_interaction.RData"))

## 3.4. Get model summary ------------------------------------------------------

# Get summary of final model
model_summary <- summary(plants_AllUrban_model1_gls)

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
birds_turnover_all_urban_15km <- birds_turnover_with_climate |>
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
birds_AllUrban_model1_gls <- gls(beta_jtu ~ all_to_urban + urban_no_change + 
                               delta_recorder_effort + log_recorder_effort + lc_time_period + temp_change + precip_change,
                             correlation = corExp(form = ~ x + y | time_numeric),  
                             data = birds_turnover_all_urban_15km_coords_time,
                             method = "REML")


# Model diagnostics (see 4.2.qmd) revealed that model is acceptable to use 
  # saving it in the final folder 
save(birds_AllUrban_model1_gls, 
     file = here("data", "models", "final", "birds_AllUrban_model1_gls.RData"))

## 4.3. Bird GLS with interaction ----------------------------------------------

# Define GLS
birds_AllUrban_model2_gls_interaction <- gls(beta_jtu ~ all_to_urban * temp_change +
                                               all_to_urban * precip_change + 
                                               urban_no_change +
                                           delta_recorder_effort + 
                                           log_recorder_effort + 
                                           lc_time_period,
                                         correlation = corExp(form = ~ x + y | time_numeric),
                                         data = birds_turnover_all_urban_15km_coords_time,
                                         method = "REML")
# Model validation looked ok enough

# Compare models based on AIC
AICtab(birds_AllUrban_model1_gls, birds_AllUrban_model2_gls_interaction, base = TRUE)
# model without interaction preferred => save interaction model in exploratory folder

# Save model
save(birds_AllUrban_model2_gls_interaction,
     file = here("data", "models", "exploratory",
                 "birds_AllUrban_model2_gls_interaction.RData"))

## 4.4. Get model summary ------------------------------------------------------ 

# Get summary of final model
model_summary <- summary(birds_AllUrban_model1_gls)

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

# Use function defined in section 14.3 of the 0_setup.R

# Create plots for each model
plot_all <- create_coef_plot_urban(AllUrban_turnover_model1_gls, 
                                   "a) All occurrences", 
                                   show_y_axis = TRUE)

plot_plants <- create_coef_plot_urban(plants_AllUrban_model1_gls, 
                                      "b) Vascular plants", 
                                      show_y_axis = FALSE)

plot_birds <- create_coef_plot_urban(birds_AllUrban_model1_gls, 
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
ggsave(filename = here("figures", "Figure6ghi_betaJTU_AllUrban_models.png"),
       plot = final_plot, width = 14, height = 7, dpi = 300, bg = "white")

ggsave(filename = here("figures", "Figure6ghi_betaJTU_AllUrban_models.pdf"),
       plot = final_plot, width = 14, height = 7, bg = "white")


# END OF SCRIPT ----------------------------------------------------------------