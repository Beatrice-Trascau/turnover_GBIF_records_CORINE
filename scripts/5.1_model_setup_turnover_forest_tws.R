##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 5.1_model_setup_bayesian_ordbeta_beta_jtu
# This script contains code which sets up the ordered beta models for temporal
# turnover (beta_jtu) across all three land-cover changes and both occurrence
# groups (vasculat plants and birds)
# Final model specification: ordered beta + dispersion submodel 
# (phi ~ predictors) + approximate spatial Gaussian process at k = 40
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

## 1.1. Load packages ---------------------------------------------------------- 

library(here)
source(here("scripts", "0_setup.R"))

# Bayesian-specific packages (not added to 0_setup.R to save load time in other scripts)
library(ordbetareg)
library(brms)
library(posterior)
library(patchwork)

## 1.2. Global model settings --------------------------------------------------

# Total number of 100m x 100m pixels in a 15km x 15km cell (150 x 150 = 22,500)
total_pixels_per_cell <- 22500

# GP basis dimensions (see section # 3.)
gp_k <- 40

# Diagnostic group
diagnostic_group <- "birds_FTWS"

# Reference period
time_ref <- "2000-2006"

# Time factor levels
time_levels <- c("2000-2006", "2006-2012", "2012-2018")

# Shared predictors
preds <- paste("lc_change_prop + lc_nochange_prop +",
                 "delta_recorder_effort + log_recorder_effort + lc_time_period")
gp_term <- paste0("gp(x, y, k = ", gp_k, ", c = 5/4)")

## 1.3. Model specifications ---------------------------------------------------

# One row per land-cover change x taxon combination
model_specs <- tibble::tribble(~id, ~group, ~transition,
                               "plants_FTWS", "plants", "FTWS",
                               "birds_FTWS",    "birds",   "FTWS",
                               "plants_TWSF",   "plants",  "TWSF",
                               "birds_TWSF",    "birds",   "TWSF",
                               "plants_urban",  "plants",  "urban",
                               "birds_urban",   "birds",   "urban")

## 1.4. Load data --------------------------------------------------------------

# Turnover data for plants
load(here("data", "derived_data", 
          "vascular_plants_turnover_all_land_cover_climate_15km.rda"))

# Turnover data for birds
load(here("data", "derived_data", 
          "bird_turnover_all_land_cover_climate_15km.rda"))

## 1.5. Helper function -------------------------------------------------------

# Build the modelling data for one group/transition. Selects the relevant land-
# cover columns, computes proportions and log recorder effort, and standardises
# the four continuous predictors. Also keeps the raw (un-standardised) predictor
# values (suffix _raw) so the prediction-figure axes can be back-transformed.
prep_turnover <- function(raw_df, transition){
  
  # pick the no-change and chnage columns for this transition
  cols <- switch(transition,
                 FTWS  = c(nochange = "Forest no change", change = "Forest to TWS"),
                 TWSF  = c(nochange = "TWS no change",    change = "TWS to Forest"),
                 urban = c(nochange = "Urban_no_change",  change = "all_to_urban"),
                 stop("Unknown transition: ", transition))
  nochange_col <- cols[["nochange"]]
  change_col   <- cols[["change"]]
  
  raw_df |>
    # pull the period-specific land-cover values into single columns
    mutate(
      lc_nochange = case_when(
        lc_time_period == "2000-2006" ~ .data[[paste0("2000-2006_", nochange_col)]],
        lc_time_period == "2006-2012" ~ .data[[paste0("2006-2012_", nochange_col)]],
        lc_time_period == "2012-2018" ~ .data[[paste0("2012-2018_", nochange_col)]],
        TRUE ~ NA_real_),
      lc_change = case_when(
        lc_time_period == "2000-2006" ~ .data[[paste0("2000-2006_", change_col)]],
        lc_time_period == "2006-2012" ~ .data[[paste0("2006-2012_", change_col)]],
        lc_time_period == "2012-2018" ~ .data[[paste0("2012-2018_", change_col)]],
        TRUE ~ NA_real_)) |>
    # keep only cells with coordinates
    filter(!is.na(x) & !is.na(y)) |>
    # derive predictors and their proportions
    mutate(lc_time_period = factor(lc_time_period, levels = time_levels),
           log_recorder_effort = log(recorder_effort),
           lc_nochange_prop    = lc_nochange / total_pixels_per_cell,
           lc_change_prop      = lc_change / total_pixels_per_cell,
           # raw copies for back-transformation in Section 5
           lc_change_prop_raw        = lc_change_prop,
           lc_nochange_prop_raw      = lc_nochange_prop,
           delta_recorder_effort_raw = delta_recorder_effort,
           log_recorder_effort_raw   = log_recorder_effort) |>
    # standardise the continuous predictors used in the models (1SD)
    mutate(across(c(lc_change_prop, lc_nochange_prop,
                    delta_recorder_effort, log_recorder_effort),
                  ~ as.numeric(scale(.)))) |>
    select(beta_jtu, x, y, lc_time_period,
           lc_change_prop, lc_nochange_prop, delta_recorder_effort, log_recorder_effort,
           ends_with("_raw"))
}

# 2. MODEL TRIALS --------------------------------------------------------------

# Two model options:
  # C1 = ordered beta, constant dispersion
  # C2 = ordered beta, dispersion submodel (phi ~ predictors)

## 2.1. Prepare dignostic group data -------------------------------------------

# Extract the modelling data for the diagnostic group
dat_d <- prep_turnover(
  if (model_specs$group[model_specs$id == diagnostic_group] == "plants")
    vascular_plants_turnover_with_climate else birds_turnover_with_climate,
  model_specs$transition[model_specs$id == diagnostic_group])

## 2.2. C1: ordered beta, constant dispersion ----------------------------------

# Fit model
C1 <- ordbetareg(formula = bf(as.formula(paste0("beta_jtu ~ ", preds, " + ", gp_term))),
                 data = dat_d, chains = 4, iter = 2000, cores = 4, seed = 1234,
                 control = list(adapt_delta = 0.99), backend = "cmdstanr")

# Write model to file
save(C1, 
     file = here("data", "models", "exploratory", 
                 "bayes_birds_FTWS_C1_ordbeta_const_spatial.RData"))

## 2.3. C2: ordered beta, dispersion submodel ----------------------------------

# Fit model
C2 <- ordbetareg(formula = bf(as.formula(paste0("beta_jtu ~ ", preds, " + ", gp_term)),
                              as.formula(paste0("phi ~ ", preds))),
                 phi_reg = "only", data = dat_d, chains  = 4, iter = 2000, cores = 4, seed = 1234,
                 control = list(adapt_delta = 0.99), backend = "cmdstanr")

# Write model to file
save(C1, 
     file = here("data", "models", "exploratory", 
                 "bayes_birds_FTWS_C2_ordbeta_disp_spatial.RData"))

## 2.4. Compare moodels --------------------------------------------------------

# Check the convergence
print(summarise_draws(C1)[, c("variable", "rhat", "ess_bulk", "ess_tail")])
print(summarise_draws(C2)[, c("variable", "rhat", "ess_bulk", "ess_tail")])

# Check the posterior predictive shape
print(pp_check(C1, ndraws = 100) + ggtitle("C1: ordered beta, constant phi"))
print(pp_check(C2, ndraws = 100) + ggtitle("C2: ordered beta, phi submodel"))

# Check loo
print(loo_compare(loo(C1), loo(C2)))


















# END OF SCRIPT ----------------------------------------------------------------