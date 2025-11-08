##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 5.2_species_richness_models
# This script models changes in species richness in relation to the three
# land cover changes
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# Turnover data
load(here("data", "derived_data", 
          "vascular_plants_all_periods_turnover_all_land_cover_chanegs_15km.rda"))
load(here("data", "derived_data", 
          "birds_all_periods_turnover_all_land_cover_chanegs_15km.rda"))

# Rename for easier handling
plants_turnover <- vascular_plants_all_periods_turnover_all_land_cover_chanegs_15km
birds_turnover <- birds_all_periods_turnover_all_land_cover_chanegs_15km

# Forest -> TWS raster
forest_tws_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_forest_tws_masked.tif"))

# TWS -> Forest raster
tws_forest_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_tws_forest_masked.tif"))

# All -> Urban raster
all_urban_15km <- rast(here("data", "derived_data", 
                            "clc_status_15km_all_urban_combined_masked.tif"))

# 2. PREPARE DATA FOR ANALYSIS -------------------------------------------------

## 2.1. Prepare vascular plants data -------------------------------------------

# Process all land cover change information
plants_all_changes <- plants_turnover |>
  # add species change metrics
  mutate(species_change = total_spp_after - total_spp_before,
         change_type = case_when(species_change > 0 ~ "Gain",
                                 species_change < 0 ~ "Loss",
                                 TRUE ~ "No Change"),
         abs_change = abs(species_change),
         taxonomic_group = "Vascular Plants") |>
  # extract land cover change information based on time period
  mutate(
    # forest to TWS transitions
    forest_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest no change`,
                                 lc_time_period == "2006-2012" ~ `2006-2012_Forest no change`,
                                 lc_time_period == "2012-2018" ~ `2012-2018_Forest no change`,
                                 TRUE ~ NA_real_),
    forest_to_tws = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest to TWS`,
                              lc_time_period == "2006-2012" ~ `2006-2012_Forest to TWS`,
                              lc_time_period == "2012-2018" ~ `2012-2018_Forest to TWS`,
                              TRUE ~ NA_real_),
    # TWS to Forest transitions
    tws_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS no change`,
                              lc_time_period == "2006-2012" ~ `2006-2012_TWS no change`,
                              lc_time_period == "2012-2018" ~ `2012-2018_TWS no change`,
                              TRUE ~ NA_real_),
    tws_to_forest = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS to Forest`,
                              lc_time_period == "2006-2012" ~ `2006-2012_TWS to Forest`,
                              lc_time_period == "2012-2018" ~ `2012-2018_TWS to Forest`,
                              TRUE ~ NA_real_),
    # all to Urban transitions
    urban_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Urban_no_change`,
                                lc_time_period == "2006-2012" ~ `2006-2012_Urban_no_change`,
                                lc_time_period == "2012-2018" ~ `2012-2018_Urban_no_change`,
                                TRUE ~ NA_real_),
    all_to_urban = case_when(lc_time_period == "2000-2006" ~ `2000-2006_all_to_urban`,
                             lc_time_period == "2006-2012" ~ `2006-2012_all_to_urban`,
                             lc_time_period == "2012-2018" ~ `2012-2018_all_to_urban`,
                             TRUE ~ NA_real_)) |>
  # create land cover change flags for filtering
  mutate(has_forest_to_tws = !is.na(forest_to_tws) & forest_to_tws > 0,
         has_tws_to_forest = !is.na(tws_to_forest) & tws_to_forest > 0,
         has_all_to_urban = !is.na(all_to_urban) & all_to_urban > 0) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort)) |>
  # Select relevant columns
  select(cell_ID, lc_time_period, x, y, taxonomic_group, species_change, 
         change_type, abs_change, total_spp_before, total_spp_after, recorder_effort,
         forest_no_change, forest_to_tws, has_forest_to_tws, delta_recorder_effort,
         tws_no_change, tws_to_forest, has_tws_to_forest, log_recorder_effort,
         urban_no_change, all_to_urban, has_all_to_urban, time_numeric)

## 2.2 Prepare birds data ------------------------------------------------------

birds_all_changes <- birds_turnover |>
  # Add species change metrics
  mutate(species_change = total_spp_after - total_spp_before,
         change_type = case_when(species_change > 0 ~ "Gain",
                                 species_change < 0 ~ "Loss",
                                 TRUE ~ "No Change"),
         abs_change = abs(species_change),
         taxonomic_group = "Birds") |>
  # extract land cover change information based on time period
  mutate(
    # forest to TWS transitions
    forest_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest no change`,
                                 lc_time_period == "2006-2012" ~ `2006-2012_Forest no change`,
                                 lc_time_period == "2012-2018" ~ `2012-2018_Forest no change`,
                                 TRUE ~ NA_real_),
    forest_to_tws = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest to TWS`,
                              lc_time_period == "2006-2012" ~ `2006-2012_Forest to TWS`,
                              lc_time_period == "2012-2018" ~ `2012-2018_Forest to TWS`,
                              TRUE ~ NA_real_),
    # TWS to Forest transitions
    tws_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS no change`,
                              lc_time_period == "2006-2012" ~ `2006-2012_TWS no change`,
                              lc_time_period == "2012-2018" ~ `2012-2018_TWS no change`,
                              TRUE ~ NA_real_),
    tws_to_forest = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS to Forest`,
                              lc_time_period == "2006-2012" ~ `2006-2012_TWS to Forest`,
                              lc_time_period == "2012-2018" ~ `2012-2018_TWS to Forest`,
                              TRUE ~ NA_real_),
    # all to Urban transitions
    urban_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Urban_no_change`,
                                lc_time_period == "2006-2012" ~ `2006-2012_Urban_no_change`,
                                lc_time_period == "2012-2018" ~ `2012-2018_Urban_no_change`,
                                TRUE ~ NA_real_),
    all_to_urban = case_when(lc_time_period == "2000-2006" ~ `2000-2006_all_to_urban`,
                             lc_time_period == "2006-2012" ~ `2006-2012_all_to_urban`,
                             lc_time_period == "2012-2018" ~ `2012-2018_all_to_urban`,
                             TRUE ~ NA_real_)) |>
  # create land cover change flags for filtering
  mutate(has_forest_to_tws = !is.na(forest_to_tws) & forest_to_tws > 0,
         has_tws_to_forest = !is.na(tws_to_forest) & tws_to_forest > 0,
         has_all_to_urban = !is.na(all_to_urban) & all_to_urban > 0) |>
  mutate(time_numeric = case_when(lc_time_period == "2000-2006" ~ 1,
                                  lc_time_period == "2006-2012" ~ 2,
                                  lc_time_period == "2012-2018" ~ 3),
         log_recorder_effort = log(recorder_effort)) |>
  # select relevant columns
  select(cell_ID, lc_time_period, x, y, taxonomic_group, species_change, 
         change_type, abs_change, total_spp_before, total_spp_after, recorder_effort,
         forest_no_change, forest_to_tws, has_forest_to_tws, delta_recorder_effort,
         tws_no_change, tws_to_forest, has_tws_to_forest, log_recorder_effort,
         urban_no_change, all_to_urban, has_all_to_urban, time_numeric)

# 3. PLANTS SPECIES CHANGES MODELS ---------------------------------------------

## 3.1. Forest -> TWS ----------------------------------------------------------

### 3.1.1. Model 1: GLS --------------------------------------------------------

# Define model
plants_sp_richness_FTWS_model1_gls <- gls(total_spp_after ~ forest_to_tws + forest_no_change +
                                            delta_recorder_effort + log_recorder_effort + lc_time_period,
                                          correlation = corExp(form = ~ x + y | time_numeric), 
                                          data = plants_all_changes,
                                          method = "REML")
# Save model
save(plants_sp_richness_FTWS_model1_gls,
     file = here("data", "models", "exploratory",
                 "plants_sp_richness_FTWS_model1_gls.RData"))

# Model diagnostic revealed major issues - trying a GAM

### 3.1.2. Model 2: GAM --------------------------------------------------------

# Convert lc_time_period to factor
plants_all_changes <- plants_all_changes |>
  mutate(lc_time_period = as.factor(lc_time_period))

# Run GAM
plants_sp_richness_FTWS_model2_GAM <- gam(species_change ~ forest_to_tws + forest_no_change +
                                            delta_recorder_effort + log_recorder_effort + 
                                            lc_time_period +
                                            s(x, y, by = lc_time_period, k = 120),
                                          data = plants_all_changes,
                                          family = gaussian(),
                                          method = "REML")

# Save model
save(plants_sp_richness_FTWS_model2_GAM,
     file = here("data", "models", "exploratory",
                 "plants_sp_richness_FTWS_model2_GAM.RData"))

# Check model
gam.check(plants_sp_richness_FTWS_model2_GAM)
# some serious issues with the model

### 3.1.3. Model 3: GAM with weights -------------------------------------------

# Calculate weights based on residual variance
  # group fitted values and calculate variance in each group
plants_all_changes <- plants_all_changes |>
  mutate(fitted_initial = fitted(plants_sp_richness_FTWS_model2_GAM),
         fitted_group = cut(fitted_initial, breaks = 10))  # create 10 groups

# Calculate variance by group
variance_by_group <- plants_all_changes |>
  mutate(resid = residuals(plants_sp_richness_FTWS_model2_GAM)) |>
  group_by(fitted_group) |>
  summarise(group_variance = var(resid, na.rm = TRUE))

# Merge back and create weights
plants_all_changes <- plants_all_changes |>
  left_join(variance_by_group, by = "fitted_group") |>
  mutate(weights = 1 / group_variance)  # inverse variance weighting

# Re-run model with weights
plants_sp_richness_FTWS_model3_GAM <- gam(species_change ~ forest_to_tws + forest_no_change + 
                                            delta_recorder_effort + log_recorder_effort + 
                                            lc_time_period +
                                            s(x, y, by = lc_time_period, k = 120),
                                          data = plants_all_changes,
                                          family = gaussian(),
                                          weights = weights,
                                          method = "REML")

# Save model
save(plants_sp_richness_FTWS_model3_GAM,
     file = here("data", "models", "exploratory",
                 "plants_sp_richness_FTWS_model3_GAM.RData"))

# Check model output
gam.check(plants_sp_richness_FTWS_model3_GAM)

# Diagnostic revealed model performance is acceptable - we will use this model
save(plants_sp_richness_FTWS_model3_GAM,
     file = here("data", "models", "final",
                 "plants_sp_richness_FTWS_model3_GAM.RData"))

## 3.2. TWS -> Forest ----------------------------------------------------------

### 3.2.1. Model 1: GAM --------------------------------------------------------

# Run GAM
plants_sp_richness_TWSF_model1_GAM <- gam(species_change ~ tws_to_forest + tws_no_change +
                                           delta_recorder_effort + log_recorder_effort + 
                                            lc_time_period +
                                            s(x, y, by = lc_time_period, k = 120),
                                          data = plants_all_changes,
                                          family = gaussian(),
                                          method = "REML")

# Save model
save(plants_sp_richness_TWSF_model1_GAM,
     file = here("data", "models", "exploratory",
                 "plants_sp_richness_TWSF_model1_GAM.RData"))

# Check model
gam.check(plants_sp_richness_TWSF_model1_GAM)
# plots are revealing some issues - will try the weighting approach

### 3.2.2. Model 2: GAM with weights -------------------------------------------

# Calculate weights based on variance groups
plants_all_changes <- plants_all_changes |>
  mutate(fitted_initial_twsf = fitted(plants_sp_richness_TWSF_model1_GAM),
         fitted_group_twsf = cut(fitted_initial_twsf, breaks = 10))

# Calculate variance by group
variance_by_group_twsf <- plants_all_changes |>
  mutate(resid_twsf = residuals(plants_sp_richness_TWSF_model1_GAM)) |>
  group_by(fitted_group_twsf) |>
  summarise(group_variance_twsf = var(resid_twsf, na.rm = TRUE))

# Merge and create weights
plants_all_changes <- plants_all_changes |>
  left_join(variance_by_group_twsf, by = "fitted_group_twsf") |>
  mutate(weights_twsf = 1 / group_variance_twsf)

# Refit with weights
plants_sp_richness_TWSF_model2_GAM <- gam(species_change ~ tws_to_forest + tws_no_change +
                                            delta_recorder_effort + log_recorder_effort + 
                                            lc_time_period +
                                            s(x, y, by = lc_time_period, k = 120),
                                          data = plants_all_changes,
                                          family = gaussian(),
                                          weights = weights_twsf,
                                          method = "REML")

# Save model output
save(plants_sp_richness_TWSF_model2_GAM,
     file = here("data", "models", "exploratory",
                 "plants_sp_richness_TWSF_model2_GAM.RData"))

# Check improved diagnostics
gam.check(plants_sp_richness_TWSF_model2_GAM)
# model has improved enough to justify use in the manuscript

# Save model output
save(plants_sp_richness_TWSF_model2_GAM,
     file = here("data", "models", "final",
                 "plants_sp_richness_TWSF_model2_GAM.RData"))

## 3.3. All -> Urban -----------------------------------------------------------

### 3.3.1. Model 1: GAM --------------------------------------------------------

# Run GAM
plants_sp_richness_AllUrban_model1_GAM <- gam(species_change ~ all_to_urban + urban_no_change +
                                            delta_recorder_effort + log_recorder_effort + 
                                            lc_time_period +
                                            s(x, y, by = lc_time_period, k = 120),
                                          data = plants_all_changes,
                                          family = gaussian(),
                                          method = "REML")

# Save model
save(plants_sp_richness_AllUrban_model1_GAM,
     file = here("data", "models", "exploratory",
                 "plants_sp_richness_AllUrban_model1_GAM.RData"))

# Check model
gam.check(plants_sp_richness_AllUrban_model1_GAM)
# plots are revealing some issues - will try the weighting approach

### 3.3.2. Model 2: GAM with weights -------------------------------------------

# Calculate weights based on variance groups
plants_all_changes <- plants_all_changes |>
  mutate(fitted_initial_allurban = fitted(plants_sp_richness_AllUrban_model1_GAM),
         fitted_group_allurban = cut(fitted_initial_allurban, breaks = 10))

# Calculate variance by group
variance_by_group_allurban <- plants_all_changes |>
  mutate(resid_allurban = residuals(plants_sp_richness_AllUrban_model1_GAM)) |>
  group_by(fitted_group_allurban) |>
  summarise(group_variance_allurban = var(resid_allurban, na.rm = TRUE))

# Merge and create weights
plants_all_changes <- plants_all_changes |>
  left_join(variance_by_group_allurban, by = "fitted_group_allurban") |>
  mutate(weights_allurban = 1 / group_variance_allurban)

# Refit with weights
plants_sp_richness_AllUrban_model2_GAM <- gam(species_change ~ all_to_urban + urban_no_change +
                                            delta_recorder_effort + log_recorder_effort + 
                                            lc_time_period +
                                            s(x, y, by = lc_time_period, k = 120),
                                          data = plants_all_changes,
                                          family = gaussian(),
                                          weights = weights_allurban,
                                          method = "REML")

# Save model output
save(plants_sp_richness_AllUrban_model2_GAM,
     file = here("data", "models", "exploratory",
                 "plants_sp_richness_AllUrban_model2_GAM.RData"))

# Check improved diagnostics
gam.check(plants_sp_richness_AllUrban_model2_GAM)
# model has improved enough to justify use in the manuscript

# Save model output
save(plants_sp_richness_AllUrban_model2_GAM,
     file = here("data", "models", "final",
                 "plants_sp_richness_AllUrban_model2_GAM.RData"))

# 4. BIRD SPECIES CHANGES MODELS -----------------------------------------------

## 4.1. Forest -> TWS ----------------------------------------------------------

### 4.1.1. Model 1: GAM --------------------------------------------------------

# Convert lc_time_period to factor
birds_all_changes <- birds_all_changes |>
  mutate(lc_time_period = as.factor(lc_time_period))

# Run GAM
birds_sp_change_FTWS_model1_GAM <- gam(species_change ~ forest_to_tws + forest_no_change +
                                            delta_recorder_effort + log_recorder_effort + 
                                            lc_time_period +
                                            s(x, y, by = lc_time_period, k = 120),
                                          data = birds_all_changes,
                                          family = gaussian(),
                                          method = "REML")

# Save model
save(birds_sp_change_FTWS_model1_GAM,
     file = here("data", "models", "exploratory",
                 "birds_sp_change_FTWS_model1_GAM.RData"))

# Check model
gam.check(birds_sp_change_FTWS_model1_GAM)
# some serious issues with the model

### 4.1.2. Model 2: GAM with weights -------------------------------------------

# Calculate weights based on residual variance
# group fitted values and calculate variance in each group
birds_all_changes <- birds_all_changes |>
  mutate(fitted_initial = fitted(birds_sp_change_FTWS_model1_GAM),
         fitted_group = cut(fitted_initial, breaks = 10))  # create 10 groups

# Calculate variance by group
variance_by_group <- birds_all_changes |>
  mutate(resid = residuals(birds_sp_change_FTWS_model1_GAM)) |>
  group_by(fitted_group) |>
  summarise(group_variance = var(resid, na.rm = TRUE))

# Merge back and create weights
birds_all_changes <- birds_all_changes |>
  left_join(variance_by_group, by = "fitted_group") |>
  mutate(weights = 1 / group_variance)  # inverse variance weighting

# Re-run model with weights
birds_sp_change_FTWS_model2_GAM <- gam(species_change ~ forest_to_tws + forest_no_change + 
                                            delta_recorder_effort + log_recorder_effort + 
                                            lc_time_period +
                                            s(x, y, by = lc_time_period, k = 120),
                                          data = birds_all_changes,
                                          family = gaussian(),
                                          weights = weights,
                                          method = "REML")

# Save model
save(birds_sp_change_FTWS_model2_GAM,
     file = here("data", "models", "exploratory",
                 "birds_sp_change_FTWS_model2_GAM.RData"))

# Check model output
gam.check(birds_sp_change_FTWS_model2_GAM)

# Diagnostic revealed model performance is acceptable - we will use this model
save(birds_sp_change_FTWS_model2_GAM,
     file = here("data", "models", "final",
                 "birds_sp_change_FTWS_model2_GAM.RData"))

## 4.2. TWS -> Forest ----------------------------------------------------------

### 4.2.1. Model 1: GAM --------------------------------------------------------

# Run GAM
birds_sp_change_TWSF_model1_GAM <- gam(species_change ~ tws_to_forest + tws_no_change +
                                            delta_recorder_effort + log_recorder_effort + 
                                            lc_time_period +
                                            s(x, y, by = lc_time_period, k = 120),
                                          data = birds_all_changes,
                                          family = gaussian(),
                                          method = "REML")

# Save model
save(birds_sp_change_TWSF_model1_GAM,
     file = here("data", "models", "exploratory",
                 "birds_sp_change_TWSF_model1_GAM.RData"))

# Check model
gam.check(birds_sp_change_TWSF_model1_GAM)
# plots are revealing some issues - will try the weighting approach

### 4.2.2. Model 2: GAM with weights -------------------------------------------

# Calculate weights based on variance groups
birds_all_changes <- birds_all_changes |>
  mutate(fitted_initial_twsf = fitted(birds_sp_change_TWSF_model1_GAM),
         fitted_group_twsf = cut(fitted_initial_twsf, breaks = 10))

# Calculate variance by group
variance_by_group_twsf <- birds_all_changes |>
  mutate(resid_twsf = residuals(birds_sp_change_TWSF_model1_GAM)) |>
  group_by(fitted_group_twsf) |>
  summarise(group_variance_twsf = var(resid_twsf, na.rm = TRUE))

# Merge and create weights
birds_all_changes <- birds_all_changes |>
  left_join(variance_by_group_twsf, by = "fitted_group_twsf") |>
  mutate(weights_twsf = 1 / group_variance_twsf)

# Refit with weights
birds_sp_change_TWSF_model2_GAM <- gam(species_change ~ tws_to_forest + tws_no_change +
                                            delta_recorder_effort + log_recorder_effort + 
                                            lc_time_period +
                                            s(x, y, by = lc_time_period, k = 120),
                                          data = birds_all_changes,
                                          family = gaussian(),
                                          weights = weights_twsf,
                                          method = "REML")

# Save model output
save(birds_sp_change_TWSF_model2_GAM,
     file = here("data", "models", "exploratory",
                 "birds_sp_change_TWSF_model2_GAM.RData"))

# Check improved diagnostics
gam.check(birds_sp_change_TWSF_model2_GAM)
# model has improved enough to justify use in the manuscript

# Save model output
save(birds_sp_change_TWSF_model2_GAM,
     file = here("data", "models", "final",
                 "birds_sp_change_TWSF_model2_GAM.RData"))

## 4.3. All -> Urban -----------------------------------------------------------

### 4.3.1. Model 1: GAM --------------------------------------------------------

# Run GAM
birds_sp_change_AllUrban_model1_GAM <- gam(species_change ~ all_to_urban + urban_no_change +
                                                delta_recorder_effort + log_recorder_effort + 
                                                lc_time_period +
                                                s(x, y, by = lc_time_period, k = 120),
                                              data = birds_all_changes,
                                              family = gaussian(),
                                              method = "REML")

# Save model
save(birds_sp_change_AllUrban_model1_GAM,
     file = here("data", "models", "exploratory",
                 "birds_sp_change_AllUrban_model1_GAM.RData"))

# Check model
gam.check(birds_sp_change_AllUrban_model1_GAM)
# plots are revealing some issues - will try the weighting approach

### 4.3.2. Model 2: GAM with weights -------------------------------------------

# Calculate weights based on variance groups
birds_all_changes <- birds_all_changes |>
  mutate(fitted_initial_allurban = fitted(birds_sp_change_AllUrban_model1_GAM),
         fitted_group_allurban = cut(fitted_initial_allurban, breaks = 10))

# Calculate variance by group
variance_by_group_allurban <- birds_all_changes |>
  mutate(resid_allurban = residuals(birds_sp_change_AllUrban_model1_GAM)) |>
  group_by(fitted_group_allurban) |>
  summarise(group_variance_allurban = var(resid_allurban, na.rm = TRUE))

# Merge and create weights
birds_all_changes <- birds_all_changes |>
  left_join(variance_by_group_allurban, by = "fitted_group_allurban") |>
  mutate(weights_allurban = 1 / group_variance_allurban)

# Refit with weights
birds_sp_change_AllUrban_model2_GAM <- gam(species_change ~ all_to_urban + urban_no_change +
                                                delta_recorder_effort + log_recorder_effort + 
                                                lc_time_period +
                                                s(x, y, by = lc_time_period, k = 120),
                                              data = birds_all_changes,
                                              family = gaussian(),
                                              weights = weights_allurban,
                                              method = "REML")

# Save model output
save(birds_sp_change_AllUrban_model2_GAM,
     file = here("data", "models", "exploratory",
                 "birds_sp_change_AllUrban_model2_GAM.RData"))

# Check improved diagnostics
gam.check(birds_sp_change_AllUrban_model2_GAM)
# model has improved enough to justify use in the manuscript

# Save model output
save(birds_sp_change_AllUrban_model2_GAM,
     file = here("data", "models", "final",
                 "birds_sp_change_AllUrban_model2_GAM.RData"))
