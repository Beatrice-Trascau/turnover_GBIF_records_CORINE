##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 3.1_land_cover_calculations
# This script contains code which calculates the land cover transitions in the
# original rasters and the aggregations (1km, 5km, 15km)
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# CORINE Status Layers
norway_corine_status_modified_stack <- rast(here("data", "derived_data",
                                                 "norway_corine_status_modified_stack.tif"))

# 2. CALCULATE LC CHANGES FOR 100M RESOLUTION ----------------------------------

# Change name of layers 
names(norway_corine_status_modified_stack) <- c("lc2000", "lc2006", 
                                                "lc2012", "lc2018")

## 2.1. Forest -> TWS ----------------------------------------------------------

# 2000-2006
transitions_2000_2006 <- analyse_forest_transition(norway_corine_status_modified_stack$"lc2000", 
                                                   norway_corine_status_modified_stack$"lc2006")

# 2006-2012
transitions_2006_2012 <- analyse_forest_transition(norway_corine_status_modified_stack$"lc2006", 
                                                   norway_corine_status_modified_stack$"lc2012")

# 2012-2018
transitions_2012_2018 <- analyse_forest_transition(norway_corine_status_modified_stack$"lc2012", 
                                                   norway_corine_status_modified_stack$"lc2018")

# Combine into single raster
clc_100m_forest_tws <- c(transitions_2000_2006, transitions_2006_2012,
                         transitions_2012_2018)

# Set categories for layers (so that you know what each value represents)
categories_forest <- data.frame(
  value = 0:3,
  class = c("Non-forest", "Forest no change", 
            "Forest to TWS", "Other forest conversion"))

# Change categories
levels(clc_100m_forest_tws) <- categories_forest

# Write raster to file
terra::writeRaster(clc_100m_forest_tws, 
                   here("data", "derived_data", 
                        "clc_status_100m_forest_tws.tif"), 
                   overwrite = TRUE)

## 2.3. TWS -> Forest ----------------------------------------------------------

# 2000-2006
tws_2000_2006 <- analyse_tws_transition(norway_corine_status_modified_stack$"lc2000",
                                        norway_corine_status_modified_stack$"lc2006")

# 2006-2012
tws_2006_2012 <- analyse_tws_transition(norway_corine_status_modified_stack$"lc2006",
                                        norway_corine_status_modified_stack$"lc2012")

# 2012-2018
tws_2012_2018 <- analyse_tws_transition(norway_corine_status_modified_stack$"lc2012",
                                        norway_corine_status_modified_stack$"lc2018")

# Combine into single raster
clc_100m_tws_forest <- c(tws_2000_2006, tws_2006_2012, tws_2012_2018)

# Set categories for layers (so that you know what each value represents)
categories_tws <- data.frame(
  value = 0:3,
  class = c("Non-TWS", "TWS no change", 
            "TWS to Forest", "Other TWS conversion"))

# Change categories
levels(clc_100m_tws_forest) <- categories_tws

# Write raster to file
terra::writeRaster(clc_100m_tws_forest, 
                   here("data", "derived_data", 
                        "clc_status_100m_tws_forest.tif"), 
                   overwrite = TRUE)