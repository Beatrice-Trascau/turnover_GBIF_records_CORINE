##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 6.1_model_setup_bayesian_ordebta_beta_jne
# This script contains code which sets up the ordered beta models for temporal
# nestedness (beta_jne) across all three land-cover changes and both occurrence
# groups (vascular plants and birds)
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
    select(beta_jne, x, y, lc_time_period,
           lc_change_prop, lc_nochange_prop, delta_recorder_effort, log_recorder_effort,
           ends_with("_raw"))
}

# 2. MODEL TRIALS --------------------------------------------------------------

# Two model options:
# C1 = ordered beta, constant dispersion
# C2 = ordered beta, dispersion submodel (phi ~ predictors)

## 2.1. Prepare diagnostic group data -------------------------------------------

# Extract the modelling data for the diagnostic group
dat_d <- prep_nestedness(
  if (model_specs$group[model_specs$id == diagnostic_group] == "plants")
    vascular_plants_turnover_with_climate else birds_turnover_with_climate,
  model_specs$transition[model_specs$id == diagnostic_group])

## 2.2. C1: ordered beta, constant dispersion ----------------------------------

# Fit model
C1 <- ordbetareg(formula = bf(as.formula(paste0("beta_jne ~ ", preds, " + ", gp_term))),
                 data = dat_d, chains = 4, iter = 2000, cores = 4, seed = 1234,
                 control = list(adapt_delta = 0.99), backend = "cmdstanr")

# Write model to file
saveRDS(C1, 
        here("data", "models", "exploratory", 
             paste0("bayes_nestedness_", diagnostic_group, "_C1_ordbeta_const_spatial.rds")))

## 2.3. C2: ordered beta, dispersion submodel ----------------------------------

# Fit model
C2 <- ordbetareg(formula = bf(as.formula(paste0("beta_jne ~ ", preds, " + ", gp_term)),
                              as.formula(paste0("phi ~ ", preds))),
                 phi_reg = "only", data = dat_d, chains  = 4, iter = 2000, cores = 4, seed = 1234,
                 control = list(adapt_delta = 0.99), backend = "cmdstanr")

# Write model to file
saveRDS(C2, 
        here("data", "models", "exploratory", 
             paste0("bayes_nestedness_", diagnostic_group, "_C2_ordbeta_disp_spatial.rds")))

## 2.4. Compare moodels --------------------------------------------------------

# Check the convergence
print(summarise_draws(C1)[, c("variable", "rhat", "ess_bulk", "ess_tail")])
print(summarise_draws(C2)[, c("variable", "rhat", "ess_bulk", "ess_tail")])

# Check the posterior predictive shape
print(pp_check(C1, ndraws = 100) + ggtitle("C1: ordered beta, constant phi"))
print(pp_check(C2, ndraws = 100) + ggtitle("C2: ordered beta, phi submodel")) 
# C2 much better at approximating where the hump in the data is

# Check the boundary: proportion of exact 0s (the zero-inflation in nestedness)
print(pp_check(C1, type = "stat", stat = function(y) mean(y == 0)) +
        ggtitle("C1: proportion of exact 0s"))
print(pp_check(C2, type = "stat", stat = function(y) mean(y == 0)) +
        ggtitle("C2: proportion of exact 0s"))

# Check loo
print(loo_compare(loo(C1), loo(C2)))

# 3. GP BASIS DIMENSION (k) SENSITIVITY AND PLATEAU ----------------------------

# Will fit the better model (C2) across a k ladder on the diagnostic group
# Then check where the spatial term (sdgp) plateaus and loo stops improving

## 3.1. Fit the k ladder -------------------------------------------------------

# Define a list of k values to try
k_ladder <- c(15, 30, 40, 50)

# Create an empty list in which to store the model outputs
k_fits <- list()

# Fit C2 for all 4 ks
for (K in k_ladder) {
  
  # fit C2 at this resolution
  message("Fitting k = ", K, " ...")
  fit_k <- ordbetareg(
    formula = bf(as.formula(paste0("beta_jne ~ ", preds,
                                   " + gp(x, y, k = ", K, ", c = 5/4)")),
                 as.formula(paste0("phi ~ ", preds))),
    phi_reg = "only",
    data    = dat_d,
    chains  = 4, iter = 2000, cores = 4, seed = 1234,
    control = list(adapt_delta = 0.99), backend = "cmdstanr")
  saveRDS(fit_k, here("data", "models", "exploratory",
                      paste0("bayes_nestedness_", diagnostic_group,
                             "_ordbeta_disp_spatial_k", K, ".rds")))
  k_fits[[as.character(K)]] <- fit_k
}

## 3.2. sdgp and focal coefficient across k ------------------------------------

# Get the spatial process SD at each k to get the plateau diagnostic
sdgp_vs_k <- bind_rows(lapply(names(k_fits), function(K) {
  ds <- summarise_draws(k_fits[[K]])
  r  <- ds[ds$variable == "sdgp_gpxy", ]
  tibble::tibble(k = as.integer(K), sdgp = r$mean, lwr = r$q5, upr = r$q95)
}))
print(sdgp_vs_k)

# Inspect land-cover change coefficient at each k 
focal_vs_k <- bind_rows(lapply(names(k_fits), function(K) {
  fx <- fixef(k_fits[[K]], probs = c(0.025, 0.975))
  tibble::tibble(k = as.integer(K),
                 lc_change_prop = fx["lc_change_prop", "Estimate"],
                 lwr = fx["lc_change_prop", "Q2.5"],
                 upr = fx["lc_change_prop", "Q97.5"])
}))
print(focal_vs_k)

# Create plot to check when sdgp vs k plateaus
p_plateau <- ggplot(sdgp_vs_k, aes(k, sdgp)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, fill = "#2196F3") +
  geom_line(color = "#2196F3") +
  geom_point(size = 2.5, color = "#2196F3") +
  labs(x = "GP basis dimension (k)", y = "sdgp (spatial process SD)") +
  theme_minimal(base_size = 12)
print(p_plateau)

# Save figure to file
ggsave(here("figures", paste0("SupplementaryFigure8_nestedness_gp_k_plateau_", diagnostic_group, ".png")),
       p_plateau, width = 7, height = 5, dpi = 300, bg = "white")

# Check predictive power across k
print(loo_compare(lapply(k_fits, loo)))

# 4. FIT FINAL MODELS ACROSS LAND-COVERS AND GROUPS ----------------------------

# The final model = ordered beta + dispersion submodel + spatial GP at k = 40)

## 4.1. Fit models in a loop ---------------------------------------------------

for (this_id in model_specs$id) {
  
  # add a progress message
  message("=== Fitting ", this_id, " ===")
  
  # define model specifications
  spec <- model_specs[model_specs$id == this_id, ]
  
  # modelling data for this group
  raw_df <- if (spec$group == "plants") vascular_plants_turnover_with_climate else birds_turnover_with_climate
  dat    <- prep_nestedness(raw_df, spec$transition)
  
  # fit the unified specification
  fit <- ordbetareg(
    formula = bf(as.formula(paste0("beta_jne ~ ", preds, " + ", gp_term)),
                 as.formula(paste0("phi ~ ", preds))),
    phi_reg = "only",
    data    = dat,
    chains  = 4, iter = 2000, cores = 4, seed = 1234,
    control = list(adapt_delta = 0.99), backend = "cmdstanr")
  
  saveRDS(fit, here("data", "models", "final",
                    paste0("bayes_nestedness_", this_id, "_ordbeta_disp_spatial.rds")))
}

## 4.2. Check convergence and posterior predictive shape -----------------------

# Loop over the final models 
for (this_id in model_specs$id) {
  
  fit <- readRDS(here("data", "models", "final",
                      paste0("bayes_nestedness_", this_id, "_ordbeta_disp_spatial.rds")))
  
  # convergence one-liner
  ds <- summarise_draws(fit)
  n_div <- sum(nuts_params(fit) |> filter(Parameter == "divergent__") |> pull(Value))
  message(sprintf("%-13s  max Rhat %.3f | min ESS %d | divergences %d",
                  this_id, max(ds$rhat, na.rm = TRUE),
                  round(min(ds$ess_bulk, na.rm = TRUE)), n_div))
  
  # posterior predictive checks
  print(pp_check(fit, ndraws = 100) + ggtitle(paste0(this_id, " - density")))
  print(pp_check(fit, type = "stat", stat = function(y) mean(y == 1)) +
          ggtitle(paste0(this_id, " - proportion of exact 1s")))
}

# 5. MODEL PREDICTION FIGURE ---------------------------------------------------

# Predicted turnover across each predictor's observed range (holds the other
# continuous predictors at their mean (0 standardised) and time at the reference
# period)
# Axes are back-transformed to original units using thr raw columns

## 5.1. Generate predictions ---------------------------------------------------

# List of the continuous predictors
cont_preds  <- c("lc_change_prop", "lc_nochange_prop",
                 "delta_recorder_effort", "log_recorder_effort")

# List of labels for the predictors
cont_labels <- c(lc_change_prop = "Land-cover change\n(proportion of cell)",
                 lc_nochange_prop = "Land-cover no change\n(proportion of cell)",
                 delta_recorder_effort = "\u0394 recorder effort",
                 log_recorder_effort = "log recorder effort")

# Create empty lists for the predictions
cont_rows <- list()
time_rows <- list()
point_rows <- list() 
presid_rows <- list() 

# Loop through predictors
for (this_id in model_specs$id) {
  
  spec <- model_specs[model_specs$id == this_id, ]
  raw_df <- if (spec$group == "plants") vascular_plants_turnover_with_climate else birds_turnover_with_climate
  dat <- prep_nestedness(raw_df, spec$transition)
  fit <- readRDS(here("data", "models", "final",
                      paste0("bayes_nestedness_", this_id, "_ordbeta_disp_spatial.rds")))
  
  # baseline row: all continuous predictors at mean (0), time at reference,
  # GP at the map centroid
  base <- tibble::tibble(lc_change_prop = 0, lc_nochange_prop = 0,
                         delta_recorder_effort = 0, log_recorder_effort = 0,
                         lc_time_period = factor(time_ref, levels = time_levels),
                         x = mean(dat$x), y = mean(dat$y))
  
  # full fitted values (all observed predictors incl GP) - for partial residuals
  yhat_full  <- colMeans(posterior_epred(fit, newdata = dat))
  resid_full <- dat$beta_jne - yhat_full
  
  # one curve per continuous predictor
  for (p in cont_preds) {
    raw <- dat[[paste0(p, "_raw")]]
    grid_raw <- seq(min(raw), max(raw), length.out = 100)
    nd <- base[rep(1, 100), ]
    nd[[p]] <- (grid_raw - mean(raw)) / sd(raw)   # standardise as in fitting
    ep <- posterior_epred(fit, newdata = nd)
    cont_rows[[length(cont_rows) + 1]] <- tibble::tibble(id = this_id, 
                                                         group = spec$group, 
                                                         transition = spec$transition,
                                                         predictor = p, 
                                                         x_value = grid_raw,
                                                         est = apply(ep, 2, median),
                                                         lwr = apply(ep, 2, quantile, 0.025),
                                                         upr = apply(ep, 2, quantile, 0.975))
    
    # overlay A: raw observed points (raw predictor vs observed nestedness)
    point_rows[[length(point_rows) + 1]] <- tibble::tibble(
      id = this_id, group = spec$group, transition = spec$transition,
      predictor = p, x_value = raw, value = dat$beta_jne)
    
    # overlay B: partial residuals (partial prediction at each observed value of
    # p, with this observation's full-model residual added back)
    nd_obs       <- base[rep(1, nrow(dat)), ]
    nd_obs[[p]]  <- dat[[p]]
    pred_partial <- colMeans(posterior_epred(fit, newdata = nd_obs))
    presid_rows[[length(presid_rows) + 1]] <- tibble::tibble(
      id = this_id, group = spec$group, transition = spec$transition,
      predictor = p, x_value = raw, value = pred_partial + resid_full)
  }
  
  # time-period predictions
  nd_t <- base[rep(1, length(time_levels)), ]
  nd_t$lc_time_period <- factor(time_levels, levels = time_levels)
  ep_t <- posterior_epred(fit, newdata = nd_t)
  time_rows[[length(time_rows) + 1]] <- tibble::tibble(id = this_id, 
                                                       group = spec$group, 
                                                       transition = spec$transition,
                                                       period = factor(time_levels, 
                                                                       levels = time_levels),
                                                       est = apply(ep_t, 2, median),
                                                       lwr = apply(ep_t, 2, quantile, 0.025),
                                                       upr = apply(ep_t, 2, quantile, 0.975))
}

# Combine continuous predictions and add display labels (taxon + transition)
cont_df <- bind_rows(cont_rows) |>
  mutate(taxon = dplyr::recode(group, plants = "Vascular plants", birds = "Birds"),
         transition = dplyr::recode(transition, FTWS = "Forest \u2192 TWS",
                                    TWSF = "TWS \u2192 Forest", urban = "All \u2192 Urban"),
         transition = factor(transition,
                             levels = c("Forest \u2192 TWS", "TWS \u2192 Forest", "All \u2192 Urban")),
         label = factor(cont_labels[predictor], levels = unname(cont_labels)))

# Combine time predictions and add display labels (taxon + transition)
time_df <- bind_rows(time_rows) |>
  mutate(taxon = dplyr::recode(group, plants = "Vascular plants", birds = "Birds"),
         transition = dplyr::recode(transition, FTWS = "Forest \u2192 TWS",
                                    TWSF = "TWS \u2192 Forest", urban = "All \u2192 Urban"),
         transition = factor(transition,
                             levels = c("Forest \u2192 TWS", "TWS \u2192 Forest", "All \u2192 Urban")),
         period = dplyr::recode(period, "2000-2006" = "2000\u20132006",
                                "2006-2012" = "2006\u20132012", "2012-2018" = "2012\u20132018"))

# Combine raw points and add display labels (taxon + transition)
points_df <- bind_rows(point_rows) |>
  mutate(taxon = dplyr::recode(group, plants = "Vascular plants", birds = "Birds"),
         transition = dplyr::recode(transition, FTWS = "Forest \u2192 TWS",
                                    TWSF = "TWS \u2192 Forest", urban = "All \u2192 Urban"),
         transition = factor(transition,
                             levels = c("Forest \u2192 TWS", "TWS \u2192 Forest", "All \u2192 Urban")),
         label = factor(cont_labels[predictor], levels = unname(cont_labels)))

# Combine partial residuals and add display labels (taxon + transition)
presid_df <- bind_rows(presid_rows) |>
  mutate(taxon = dplyr::recode(group, plants = "Vascular plants", birds = "Birds"),
         transition = dplyr::recode(transition, FTWS = "Forest \u2192 TWS",
                                    TWSF = "TWS \u2192 Forest", urban = "All \u2192 Urban"),
         transition = factor(transition,
                             levels = c("Forest \u2192 TWS", "TWS \u2192 Forest", "All \u2192 Urban")),
         label = factor(cont_labels[predictor], levels = unname(cont_labels)))

# Define taxon colours
taxon_cols <- c("Vascular plants" = "#1B7837", "Birds" = "#762A83")

## 5.2. Predictor figure for land-cover and recorder effort --------------------

# Transitions (rows) x predictor (cols), taxa in the same panel
p_cont <- ggplot(cont_df, aes(x_value, est, color = taxon, fill = taxon)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.7) +
  facet_grid(transition ~ label, scales = "free_x") +
  scale_color_manual(values = taxon_cols) +
  scale_fill_manual(values = taxon_cols) +
  labs(x = NULL, y = "Predicted nestedness (\u03b2_jne)", color = NULL, fill = NULL) +
  theme_bw(base_size = 14, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "top")

# Inspect figure
print(p_cont)

# Save figure to file
ggsave(here("figures", "Fig3a_predicted_nestedness_continuous.png"),
       p_cont, width = 11, height = 6.5, dpi = 300, bg = "white")

### 5.2.1. Overlay A: raw observed points --------------------------------------

# Plot the same figure with raw data points added
p_cont_raw <- ggplot(cont_df, aes(x_value, est, color = taxon, fill = taxon)) +
  geom_point(data = points_df, aes(x = x_value, y = value, color = taxon),
             inherit.aes = FALSE, alpha = 0.08, size = 0.4) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.7) +
  facet_grid(transition ~ label, scales = "free_x") +
  scale_color_manual(values = taxon_cols) +
  scale_fill_manual(values = taxon_cols) +
  labs(x = NULL, y = "Predicted nestedness (\u03b2_jne)", color = NULL, fill = NULL) +
  theme_bw(base_size = 14, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "top")

# Inspect figure
print(p_cont_raw)

### 5.2.2. Overlay B: partial residuals ----------------------------------------

# Plot the same figure as above with partial residuals: each observation's full-model residual is
# added to the partial-effect prediction, which strips out the other predictors
# and the spatial term so the points scatter around the curve by unexplained
# variation only.
p_cont_presid <- ggplot(cont_df, aes(x_value, est, color = taxon, fill = taxon)) +
  geom_point(data = presid_df, aes(x = x_value, y = value, color = taxon),
             inherit.aes = FALSE, alpha = 0.08, size = 0.4) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.7) +
  facet_grid(transition ~ label, scales = "free_x") +
  scale_color_manual(values = taxon_cols) +
  scale_fill_manual(values = taxon_cols) +
  labs(x = NULL, y = "Predicted nestedness (\u03b2_jne)", color = NULL, fill = NULL) +
  theme_bw(base_size = 14, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "top")

# Inspect figure
print(p_cont_presid)

## 5.3. Time-period predictions ------------------------------------------------

# Transition (columns), taxa as dodge points + intervals
p_time <- ggplot(time_df, aes(period, est, color = taxon)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  position = position_dodge(width = 0.4), size = 0.45) +
  facet_wrap(~ transition, nrow = 1) +
  scale_color_manual(values = taxon_cols) +
  labs(x = "Time period", y = "Predicted nestedness (\u03b2_jne)", color = NULL) +
  theme_bw(base_size = 14, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "top")

# Inspect figure
print(p_time)

# Save figure to file
ggsave(here("figures", "Fig3b_predicted_nestedness_time.png"),
       p_time, width = 9, height = 3.5, dpi = 300, bg = "white")

## 5.4. Combine figure (raw points) --------------------------------------------

# Combine the two panels
combined_fig_raw <- p_cont_raw / p_time +
  plot_layout(heights = c(3, 1.4)) +
  plot_annotation(tag_levels = "a")

# Save figure
ggsave(here("figures", "Fig3_predicted_nestedness_all_predictors.png"),
       combined_fig_raw, width = 11, height = 9.5, dpi = 300, bg = "white")
ggsave(here("figures", "Fig2_predicted_turnover_all_predictors.pdf"),
       combined_fig_raw, width = 11, height = 9.5, dpi = 300, bg = "white")

## 5.5. Combine figure (partial residuals) -------------------------------------

# Combine the two panels
combined_fig_residuals <- p_cont_presid / p_time +
  plot_layout(heights = c(3, 1.4)) +
  plot_annotation(tag_levels = "a")

# Save figure
ggsave(here("figures", "SupplementaryFigure8_predicted_nestedness_all_predictors_partial_residuals.png"),
       combined_fig_residuals, width = 11, height = 9.5, dpi = 300, bg = "white")


# END OF SCRIPT ----------------------------------------------------------------