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

# Turnover data for plants
load(here("data", "derived_data", 
          "vascular_plants_turnover_all_land_cover_climate_15km.rda"))

# Turnover data for birds
load(here("data", "derived_data", 
          "bird_turnover_all_land_cover_climate_15km.rda"))

# 2. FOREST -> TWS TURNOVER (BETA_JTU) -----------------------------------------

# Total number of 100m x 100m pixels in a 15km x 15km cell
# 15,000m / 100m = 150 pixels per side
# 150 x 150 = 22,500 total pixels per cell
total_pixels_per_cell <- 22500


# END OF SCRIPT ----------------------------------------------------------------