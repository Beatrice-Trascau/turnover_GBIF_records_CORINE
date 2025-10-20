##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 4.2_model_setup_turnover_tws_forest
# This script contains code which sets up the models exploring the impact of 
# Forest -> TWS land cover transition on temporal turnover for all occurrences,
# plant-only occurrences and bird-only occurrences
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# Turnover data
load(here("data", "derived_data", 
          "all_periods_turnover_all_land_cover_chanegs_15km.rda"))

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


# END OF SCRIPT ----------------------------------------------------------------