##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 3.2_turnover_land_cover_exploration_and_models
# This script contains code which plots and models the relationship between
# turnover and land cover changes
##----------------------------------------------------------------------------##

# 1. PREPARE DATA --------------------------------------------------------------

## 1.1. Load data --------------------------------------------------------------
turnover <- readRDS(here("data", "derived_data", 
                         "jaccard_temporal_turnover_with_land_cover.rds"))

## 1.2. Prepare df for analyisi ------------------------------------------------

# Add land cover change columns
turnover1 <- turnover |>
  mutate(LC2000_2006_change = ifelse(LC2000 != LC2006, "Y", "N"),
         LC2006_2012_change = ifelse(LC2006 != LC2012, "Y", "N"),
         LC2012_2018_change = ifelse(LC2012 != LC2018, "Y", "N"))

# Convert to long format
turnover_long <- turnover1 |>
  pivot_longer(cols = c("LC2000_2006_change", "LC2006_2012_change", "LC2012_2018_change"),
               names_to = "change_period", values_to = "LC_change") |>
  # Match change_period to start_period
  mutate(start_period = case_when(
    change_period == "LC2000_2006_change" ~ "1997-2000",
    change_period == "LC2006_2012_change" ~ "2003-2006",
    change_period == "LC2012_2018_change" ~ "2009-2012",
    TRUE ~ as.character(start_period)))

# Convert numerical LC values to factorial
turnover_lc <- turnover_long |>
  mutate(LC2000 = case_when(LC2000 == 1 ~ "urban",
                            LC2000 == 80 ~ "complex_agri",
                            LC2000 == 103 ~ "agri_sig_veg",
                            LC2000 == 250 ~ "forests",
                            LC2000 == 380 ~ "moors_heath_grass",
                            LC2000 == 590 ~ "woodland_shrub",
                            LC2000 == 711 ~ "sparse_veg",
                            LC2000 == NA ~ "other"),
         LC2006 = case_when(LC2006 == 1 ~ "urban",
                            LC2006 == 80 ~ "complex_agri",
                            LC2006 == 103 ~ "agri_sig_veg",
                            LC2006 == 250 ~ "forests",
                            LC2006 == 380 ~ "moors_heath_grass",
                            LC2006 == 590 ~ "woodland_shrub",
                            LC2006 == 711 ~ "sparse_veg",
                            LC2006 == NA ~ "other"),
         LC2012 = case_when(LC2012 == 1 ~ "urban",
                            LC2012 == 80 ~ "complex_agri",
                            LC2012 == 103 ~ "agri_sig_veg",
                            LC2012 == 250 ~ "forests",
                            LC2012 == 380 ~ "moors_heath_grass",
                            LC2012 == 590 ~ "woodland_shrub",
                            LC2012 == 711 ~ "sparse_veg",
                            LC2012 == NA ~ "other"),
         LC2018 = case_when(LC2018 == 1 ~ "urban",
                            LC2018 == 80 ~ "complex_agri",
                            LC2018 == 103 ~ "agri_sig_veg",
                            LC2018 == 250 ~ "forests",
                            LC2018 == 380 ~ "moors_heath_grass",
                            LC2018 == 590 ~ "woodland_shrub",
                            LC2018 == 711 ~ "sparse_veg",
                            LC2018 == NA ~ "other"))

# Change values to run beta model with the formula
#  ( Y * (N - 1) + 0.5 ) / N
# Y = response variable, N = sample size
N <- nrow(turnover_lc)

turnover_model <- turnover_lc |>
  mutate(jaccard_dissimilarity_beta = ( jaccard_dissimilarity * (N - 1) + 0.5 ) / N)

# Set up model
model1.1_SSB <- glmmTMB(jaccard_dissimilarity_beta ~ LC_change * start_period + (1 | SSBID),
                        family = beta_family(),
                        data = turnover_model)

# Save model output to file to save time next time
save(model1.1_SSB, file = here::here("data", "models", "model1.1_SSB.RData"))

# Beta GAM
model1.2 <- gam(jaccard_dissimilarity_beta ~ LC_change * start_period, 
             family = betar(link = "logit"), 
             data = turnover_model)

# Save model output to file to save time next time
save(model1.2, file = here::here("data", "models", "model1.2.RData"))
