##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 1.2_CORINE_aggregation
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

# Set categories for layers (so that you know what each value represents)
categories_forest <- data.frame(
  value = 0:3,
  class = c("Non-forest", "Forest no change", 
            "Forest to TWS", "Other forest conversion"))

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
  class = c("Non-TWS", "TWS no change", 
            "TWS to Forest", "Other TWS conversion"))

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
  class = c("Urban no change",
            "Forest to urban",
            "TWS to urban", 
            "Complex agriculture to urban",
            "Agriculture & vegetation to urban",
            "Moors, heathland & grassland to urban",
            "Sparse vegetation to urban",
            "No urban conversion"))

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

# 3. CALCULATE LC CHANGES FOR HIGHER RESOLUTIONS -------------------------------

## 3.1. Forest -> TWS ----------------------------------------------------------

# Aggregate to 1km, 5km, 15km
forest_tws_1km <- aggregate_transitions(clc_100m_forest_tws, 10)
forest_tws_5km <- aggregate_transitions(clc_100m_forest_tws, 50)
forest_tws_15km <- aggregate_transitions(clc_100m_forest_tws, 150)

# Save aggregated layers
writeRaster(forest_tws_1km,
            here("data", "derived_data", "clc_status_1km_forest_tws.tif"),
            overwrite = TRUE)
writeRaster(forest_tws_5km,
            here("data", "derived_data", "clc_status_5km_forest_tws.tif"),
            overwrite = TRUE)
writeRaster(forest_tws_15km,
            here("data", "derived_data", "clc_status_15km_forest_tws.tif"),
            overwrite = TRUE)

## 3.2. TWS -> Forest ----------------------------------------------------------

# Aggregate to 1km, 5km, 15km
tws_forest_1km <- aggregate_transitions(clc_100m_tws_forest, 10)
tws_forest_5km <- aggregate_transitions(clc_100m_tws_forest, 50)
tws_forest_15km <- aggregate_transitions(clc_100m_tws_forest, 150)

# Save aggregated layers
writeRaster(tws_forest_1km,
            here("data", "derived_data", "clc_status_1km_tws_forest.tif"),
            overwrite = TRUE)
writeRaster(tws_forest_5km,
            here("data", "derived_data", "clc_status_5km_tws_forest.tif"),
            overwrite = TRUE)
writeRaster(tws_forest_15km,
            here("data", "derived_data", "clc_status_15km_tws_forest.tif"),
            overwrite = TRUE)

## 3.3. All -> Urban -----------------------------------------------------------

# Aggregate to 1km, 5km, 15km
all_urban_1km <- aggregate_transitions(clc_100m_all_urban, 10)
all_urban_5km <- aggregate_transitions(clc_100m_all_urban, 50)
all_urban_15km <- aggregate_transitions(clc_100m_all_urban, 150)

# Save aggregated layers
writeRaster(all_urban_1km,
            here("data", "derived_data", "clc_status_1km_all_urban.tif"),
            overwrite = TRUE)
writeRaster(all_urban_5km,
            here("data", "derived_data", "clc_status_5km_all_urban.tif"),
            overwrite = TRUE)
writeRaster(all_urban_15km,
            here("data", "derived_data", "clc_status_15km_all_urban.tif"),
            overwrite = TRUE)

# Check names of layer
names(all_urban_1km)

# Create rasters that sum up all transitions to urban for each time period
# Define list of resolutions
urban_rasters <- list(
  "1km" = all_urban_1km,
  "5km" = all_urban_5km,
  "15km" = all_urban_15km)

# Loop through each resolution
for(res in resolutions) {
  # Load the appropriate raster
  all_urban <- urban_rasters[[res]]
  
  # Create combined urban transitions for each time period
  urban_2000_2006 <- sum(all_urban[[2:7]])    # layers 2-7 for 2000-2006
  urban_2006_2012 <- sum(all_urban[[10:15]])  # layers 10-15 for 2006-2012
  urban_2012_2018 <- sum(all_urban[[18:23]])  # layers 18-23 for 2012-2018
  
  # Combine into a single raster with 3 layers
  all_urban_combined <- c(urban_2000_2006, urban_2006_2012, urban_2012_2018)
  
  # Set names for the layers
  names(all_urban_combined) <- c("2000-2006_all_to_urban", 
                                 "2006-2012_all_to_urban", 
                                 "2012-2018_all_to_urban")
  
  # Save the new raster
  writeRaster(all_urban_combined, 
              filename = here("data", "derived_data", 
                              paste0("clc_status_", res, "_all_urban_combined.tif")),
              overwrite = TRUE)
  
  # Print progress
  print(paste("Completed processing", res, "resolution"))
}

# END OF SCRIPT ----------------------------------------------------------------