##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 1.2_CORINE_aggregation
# This script contains code which calculates the land cover transitions in the
# original rasters and the aggregation (15km)
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

# Set categories for layers (so that you know what each value represents)
categories_forest <- data.frame(
  value = 0:3,
  class = c("Non_forest", "Forest_no_change", 
            "Forest_to_TWS", "Other_forest_conversion"))

# Combine into single raster
clc_100m_forest_tws <- c(transitions_2000_2006, transitions_2006_2012,
                         transitions_2012_2018)

# Change categories
for(i in 1:nlyr(clc_100m_forest_tws)) {
  levels(clc_100m_forest_tws[[i]]) <- categories_forest
}

# Set names for layers
names(clc_100m_forest_tws) <- c("2000-2006", "2006-2012", "2012-2018")

# Write raster to file
terra::writeRaster(clc_100m_forest_tws, 
                   here("data", "derived_data", 
                        "clc_status_100m_forest_tws.tif"), 
                   overwrite = TRUE)

## 2.2. TWS -> Forest ----------------------------------------------------------

# 2000-2006
tws_2000_2006 <- analyse_tws_transition(norway_corine_status_modified_stack$"lc2000",
                                        norway_corine_status_modified_stack$"lc2006")

# 2006-2012
tws_2006_2012 <- analyse_tws_transition(norway_corine_status_modified_stack$"lc2006",
                                        norway_corine_status_modified_stack$"lc2012")

# 2012-2018
tws_2012_2018 <- analyse_tws_transition(norway_corine_status_modified_stack$"lc2012",
                                        norway_corine_status_modified_stack$"lc2018")

# Set categories for layers (so that you know what each value represents)
categories_tws <- data.frame(
  value = 0:3,
  class = c("Non_TWS", "TWS_no_change", 
            "TWS_to_Forest", "Other_TWS_conversion"))

# Combine into single raster
clc_100m_tws_forest <- c(tws_2000_2006, tws_2006_2012, tws_2012_2018)

# Change categories
for(i in 1:nlyr(clc_100m_tws_forest)) {
  levels(clc_100m_tws_forest[[i]]) <- categories_tws
}

# Set names for layers
names(clc_100m_tws_forest) <- c("2000-2006", "2006-2012", "2012-2018")

# Write raster to file
terra::writeRaster(clc_100m_tws_forest, 
                   here("data", "derived_data", 
                        "clc_status_100m_tws_forest.tif"), 
                   overwrite = TRUE)

## 2.3. All -> Urban -----------------------------------------------------------

# Calculate cells where land covers (any of them) are converted to urban areas

# 2000-2006
urban_2000_2006 <- analyse_urban_conversion(norway_corine_status_modified_stack$"lc2000",
                                            norway_corine_status_modified_stack$"lc2006")

# 2006-2012
urban_2006_2012 <- analyse_urban_conversion(norway_corine_status_modified_stack$"lc2006",
                                            norway_corine_status_modified_stack$"lc2012")

# 2012-2018
urban_2012_2018 <- analyse_urban_conversion(norway_corine_status_modified_stack$"lc2012",
                                            norway_corine_status_modified_stack$"lc2018")

# Set categories for layers (so that you know what each value represents)
categories_urban <- data.frame(
  value = 0:7,
  class = c("Urban_no_change",
            "Forest_to_urban",
            "TWS_to_urban", 
            "Complex_agriculture_to_urban",
            "Agriculture_vegetation_to_urban",
            "Moors_heathland_grassland_to_urban",
            "Sparse_vegetation_to_urban",
            "No_urban_conversion"))

# Combine into single raster
clc_100m_all_urban <- c(urban_2000_2006, urban_2006_2012, urban_2012_2018)

# Change categories
for(i in 1:nlyr(clc_100m_all_urban)) {
  levels(clc_100m_all_urban[[i]]) <- categories_urban
}

# Set names for layers
names(clc_100m_all_urban) <- c("2000-2006", "2006-2012", "2012-2018")

# Write raster to file
terra::writeRaster(clc_100m_all_urban, 
                   here("data", "derived_data", 
                        "clc_status_100m_all_urban.tif"), 
                   overwrite = TRUE)

# 3. AGGREGATE CHANGES TO 15 KM ------------------------------------------------

## 3.1. Forest -> TWS ----------------------------------------------------------

# Aggregate to 15km
forest_tws_15km <- aggregate_transitions(clc_100m_forest_tws, 150)

# Save aggregated layer
writeRaster(forest_tws_15km,
            here("data", "derived_data", "clc_status_15km_forest_tws.tif"),
            overwrite = TRUE)

## 3.2. TWS -> Forest ----------------------------------------------------------

# Aggregate to 15
tws_forest_15km <- aggregate_transitions(clc_100m_tws_forest, 150)

# Save aggregated layer
writeRaster(tws_forest_15km,
            here("data", "derived_data", "clc_status_15km_tws_forest.tif"),
            overwrite = TRUE)

## 3.3. All -> Urban -----------------------------------------------------------

# Aggregate to 15
all_urban_15km <- aggregate_transitions(clc_100m_all_urban, 150)

# Check names of layer
names(all_urban_15km)

# Create rasters that sum up all transitions to urban for each time period
all_urban_2000_2006 <- sum(all_urban_15km[[2:7]])    # layers 2-7 for 2000-2006
all_urban_2006_2012 <- sum(all_urban_15km[[10:15]])  # layers 10-15 for 2006-2012
all_urban_2012_2018 <- sum(all_urban_15km[[18:23]])  # layers 18-23 for 2012-2018

# Combine into single layer
all_urban_15km_combined <- c(all_urban_15km[[1]], all_urban_2000_2006, all_urban_15km[[8]],
                             all_urban_15km[[9]], all_urban_2006_2012, all_urban_15km[[16]],
                             all_urban_15km[[17]], all_urban_2012_2018, all_urban_15km[[24]])


# Change names
names(all_urban_15km_combined) <- c("2000-2006_Urban_no_change",
                                    "2000-2006_all_to_urban",
                                    "2000-2006_No_urban_conversion",
                                    "2006-2012_Urban_no_change",
                                    "2006-2012_all_to_urban",
                                    "2006-2012_No_urban_conversion",
                                    "2012-2018_Urban_no_change",
                                    "2012-2018_all_to_urban",
                                    "2012-2018_No_urban_conversion")

# Save to file
writeRaster(all_urban_15km_combined,
            here("data", "derived_data", "clc_status_15km_all_urban_combined.tif"),
            overwrite = TRUE)

# END OF SCRIPT ----------------------------------------------------------------