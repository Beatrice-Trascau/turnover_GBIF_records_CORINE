##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 5.1_species_gain_loss_forest_tws
# This script calculates and visualises species gain and loss for birds and
# vascular plants in areas with forest to TWS land cover changes
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

## 1.1. Vascular plants data ---------------------------------------------------

# Load data
load(here("data", "derived_data", 
          "vascular_plants_all_periods_turnover_all_land_cover_chanegs_15km.rda"))

# Rename for easier handling
plants_turnover <- vascular_plants_all_periods_turnover_all_land_cover_chanegs_15km