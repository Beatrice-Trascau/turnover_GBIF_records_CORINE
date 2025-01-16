##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 3.2_turnover_aggregated_land_cover
# This script contains code which calculates the temporal turnover of GBIF
# records in the aggregate CORINE land cover layers (1km, 5km and 15km)
##----------------------------------------------------------------------------##

source(here("scripts", "0_setup.R"))

# 1. LOAD DATA -----------------------------------------------------------------

# CORINE Aggregated Layers
norway_corine_status_modified_stack <- rast(here("data",
                                                 "derived_data",
                                                 "norway_corine_status_modified_stack.tif"))
corine_1km <- rast(here("data", "derived_data", "corine_1km.tif"))
corine_5km <- rast(here("data", "derived_data", "corine_5km.tif"))
corine_15km <- rast(here("data", "derived_data", "corine_15km.tif"))

# Cleaned occurrence records
occurrences_norway <- fread(here("data", "derived_data", 
                                 "cleaned_occurrences_july24.txt"))

# 2. PREPARE DATA FOR ANALYSIS -------------------------------------------------

## 2.1. Prepare occurrence records ---------------------------------------------

# Convert occurrence records df to sf object
occurrences_sf <- st_as_sf(occurrences_norway,
                           coords = c("decimalLongitude", "decimalLatitude"),
                           crs = 4326)
# Re-project 
occurrences_sf_reprojected <- st_transform(occurrences_sf, 
                                           st_crs(norway_corine_status_modified_stack[[1]]))

# Define periods 
clean_occurrences_sf <- occurrences_sf_reprojected |>
  mutate(period = case_when(
    year %in% c(1997:2000) ~ "1997-2000",
    year %in% c(2006:2009) ~ "2006-2009",
    year %in% c(2003:2006) ~ "2003-2006",
    year %in% c(2012:2015) ~ "2012-2015",
    year %in% c(2009:2012) ~ "2009-2012",
    year %in% c(2015:2018) ~ "2015-2018",
    TRUE ~ NA_character_)) |>
  filter(!is.na(period))

## 2.2. Prepare ID rasters from the aggregated raster --------------------------

# 1km
land_cover_id_1km <- corine_1km[[1]]
land_cover_id_1km[] <- 1:ncell(land_cover_id_1km)

# 5km
land_cover_id_5km <- corine_5km[[1]]
land_cover_id_5km[] <- 1:ncell(land_cover_id_5km)

# 15km
land_cover_id_15km <- corine_15km[[1]]
land_cover_id_15km[] <- 1:ncell(land_cover_id_15km)

# 3. REPROJECT AND EXTRACT CELL IDS FOR OCCURRENCES ----------------------------

## 3.1. Extract cell_id to each occurrence -------------------------------------

# Extract coordinates from clean occurrences
coords <- st_coordinates(clean_occurrences_sf)

# Extract land cover cell IDs for each occurrence 
cell_ids_1km <- terra::extract(land_cover_id_1km, coords)[, 1]
cell_ids_5km <- terra::extract(land_cover_id_5km, coords)[, 1]
cell_ids_15km <- terra::extract(land_cover_id_15km, coords)[, 1]

## 3.2. Add cell_ids to the occurrence df --------------------------------------

# 1km
clean_occurrences_1km <- clean_occurrences_sf |>
  mutate(cell = cell_ids_1km)

# 5km
clean_occurrences_5km <- clean_occurrences_sf |>
  mutate(cell = cell_ids_5km)

# 15km
clean_occurrences_15km <- clean_occurrences_sf |>
  mutate(cell = cell_ids_15km)

# 4. EXCLUDE CELLS WITH LESS THAN 3 SPECIES ------------------------------------

# Using custom function, see script 0_setup.R

# 1km
unique_occurrences_1km <- filter_cells_with_species_each_period(clean_occurrences_1km)

# 5km
unique_occurrences_5km <- filter_cells_with_species_each_period(clean_occurrences_5km)

# 15km
unique_occurrences_15km <- filter_cells_with_species_each_period(clean_occurrences_15km)

# 5. CALCULATE TEMPORAL TURNOVER -----------------------------------------------

# Using custom function, see script 0_setup.R

# Define period combinations 
period_combinations <- list(
  c("1997-2000", "2006-2009"),
  c("2003-2006", "2012-2015"),
  c("2009-2012", "2015-2018"))

# Calculate Jaccard's dissimilarity for each aggregation
temporal_turnover_1km <- calculate_turnover_for_resolutions(unique_occurrences_1km)
temporal_turnover_5km <- calculate_turnover_for_resolutions(unique_occurrences_5km)
temporal_turnover_15km <- calculate_turnover_for_resolutions(unique_occurrences_15km)

# 6. ADD CELL COORDINATES FOR AGGREGATED CELLS ---------------------------------

# 6.1. Create df with % changes for aggregated rasters -------------------------

# 1km
corine_1km_df <- as.data.frame(corine_1km, xy = TRUE) |>
  mutate(cell = row_number(),
         urban_change2000.2006 = (urban2006 - urban2000)/100,
         complex_agri_change2000.2006 = (complex_agri2006 - complex_agri2000)/100,
         agri_sig_veg_change2000.2006 = (agri_sig_veg2006 - agri_sig_veg2000)/100,
         forest_change2000.2006 = (forest2006 - forest2000)/100,
         moors_change2000.2006 = (moors2006 - moors2000)/100,
         woodland_change2000.2006 = (woodland2006 - woodland2000)/100,
         sparse_veg_change2000.2006 = (sparse_veg2006 - sparse_veg2000)/100,
         urban_change2006.2012 = (urban2012 - urban2006)/100,
         complex_agri_change2006.2012 = (complex_agri2012 - complex_agri2006)/100,
         agri_sig_veg_change2006.2012 = (agri_sig_veg2012 - agri_sig_veg2006)/100,
         forest_change2006.2012 = (forest2012 - forest2006)/100,
         moors_change2006.2012 = (moors2012 - moors2006)/100,
         woodland_change2006.2012 = (woodland2012 - woodland2006)/100,
         sparse_veg_change2006.2012 = (sparse_veg2012 - sparse_veg2006)/100,
         urban_change2012.2018 = (urban2018 - urban2012)/100,
         complex_agri_change2012.2018 = (complex_agri2018 - complex_agri2012)/100,
         agri_sig_veg_change2012.2018 = (agri_sig_veg2018 - agri_sig_veg2012)/100,
         forest_change2012.2018 = (forest2018 - forest2012)/100,
         moors_change2012.2018 = (moors2018 - moors2012)/100,
         woodland_change2012.2018 = (woodland2018 - woodland2012)/100,
         sparse_veg_change2012.2018 = (sparse_veg2018 - sparse_veg2012)/100) |>
  select(cell, starts_with("urban_change"), starts_with("complex_agri_change"),
         starts_with("agri_sig_veg_change"), starts_with("forest_change"),
         starts_with("moors_change"), starts_with("woodland_change"),
         starts_with("sparse_veg_change"))

# 5km
corine_5km_df <- as.data.frame(corine_5km, xy = TRUE) |>
  mutate(cell = row_number(),
         urban_change2000.2006 = (urban2006 - urban2000)/100,
         complex_agri_change2000.2006 = (complex_agri2006 - complex_agri2000)/100,
         agri_sig_veg_change2000.2006 = (agri_sig_veg2006 - agri_sig_veg2000)/100,
         forest_change2000.2006 = (forest2006 - forest2000)/100,
         moors_change2000.2006 = (moors2006 - moors2000)/100,
         woodland_change2000.2006 = (woodland2006 - woodland2000)/100,
         sparse_veg_change2000.2006 = (sparse_veg2006 - sparse_veg2000)/100,
         urban_change2006.2012 = (urban2012 - urban2006)/100,
         complex_agri_change2006.2012 = (complex_agri2012 - complex_agri2006)/100,
         agri_sig_veg_change2006.2012 = (agri_sig_veg2012 - agri_sig_veg2006)/100,
         forest_change2006.2012 = (forest2012 - forest2006)/100,
         moors_change2006.2012 = (moors2012 - moors2006)/100,
         woodland_change2006.2012 = (woodland2012 - woodland2006)/100,
         sparse_veg_change2006.2012 = (sparse_veg2012 - sparse_veg2006)/100,
         urban_change2012.2018 = (urban2018 - urban2012)/100,
         complex_agri_change2012.2018 = (complex_agri2018 - complex_agri2012)/100,
         agri_sig_veg_change2012.2018 = (agri_sig_veg2018 - agri_sig_veg2012)/100,
         forest_change2012.2018 = (forest2018 - forest2012)/100,
         moors_change2012.2018 = (moors2018 - moors2012)/100,
         woodland_change2012.2018 = (woodland2018 - woodland2012)/100,
         sparse_veg_change2012.2018 = (sparse_veg2018 - sparse_veg2012)/100) |>
  select(cell, starts_with("urban_change"), starts_with("complex_agri_change"),
         starts_with("agri_sig_veg_change"), starts_with("forest_change"),
         starts_with("moors_change"), starts_with("woodland_change"),
         starts_with("sparse_veg_change"))

# 15 km
corine_15km_df <- as.data.frame(corine_15km, xy = TRUE) |>
  mutate(cell = row_number(),
         urban_change2000.2006 = (urban2006 - urban2000)/100,
         complex_agri_change2000.2006 = (complex_agri2006 - complex_agri2000)/100,
         agri_sig_veg_change2000.2006 = (agri_sig_veg2006 - agri_sig_veg2000)/100,
         forest_change2000.2006 = (forest2006 - forest2000)/100,
         moors_change2000.2006 = (moors2006 - moors2000)/100,
         woodland_change2000.2006 = (woodland2006 - woodland2000)/100,
         sparse_veg_change2000.2006 = (sparse_veg2006 - sparse_veg2000)/100,
         urban_change2006.2012 = (urban2012 - urban2006)/100,
         complex_agri_change2006.2012 = (complex_agri2012 - complex_agri2006)/100,
         agri_sig_veg_change2006.2012 = (agri_sig_veg2012 - agri_sig_veg2006)/100,
         forest_change2006.2012 = (forest2012 - forest2006)/100,
         moors_change2006.2012 = (moors2012 - moors2006)/100,
         woodland_change2006.2012 = (woodland2012 - woodland2006)/100,
         sparse_veg_change2006.2012 = (sparse_veg2012 - sparse_veg2006)/100,
         urban_change2012.2018 = (urban2018 - urban2012)/100,
         complex_agri_change2012.2018 = (complex_agri2018 - complex_agri2012)/100,
         agri_sig_veg_change2012.2018 = (agri_sig_veg2018 - agri_sig_veg2012)/100,
         forest_change2012.2018 = (forest2018 - forest2012)/100,
         moors_change2012.2018 = (moors2018 - moors2012)/100,
         woodland_change2012.2018 = (woodland2018 - woodland2012)/100,
         sparse_veg_change2012.2018 = (sparse_veg2018 - sparse_veg2012)/100) |>
  select(cell, starts_with("urban_change"), starts_with("complex_agri_change"),
         starts_with("agri_sig_veg_change"), starts_with("forest_change"),
         starts_with("moors_change"), starts_with("woodland_change"),
         starts_with("sparse_veg_change"))

# 6.2. Add cell coordinates and LC values to temporal turnover dfs -------------
# Using custom function, see script 0_setup.R

# Add cell coordinates for aggregated CORINE cells
temporal_turnover_1km <- add_cell_coordinates_and_land_cover(temporal_turnover_1km,
                                                             land_cover_id_1km,
                                                             corine_1km_df)

temporal_turnover_5km <- add_cell_coordinates_and_land_cover(temporal_turnover_5km,
                                                             land_cover_id_5km,
                                                             corine_5km_df)

temporal_turnover_15km <- add_cell_coordinates_and_land_cover(temporal_turnover_15km,
                                                             land_cover_id_15km,
                                                             corine_15km_df)
# Write results to file
saveRDS(temporal_turnover_1km, here("data", "derived_data", 
                                    "jaccard_temporal_turnover_1km.rds"))
saveRDS(temporal_turnover_5km, here("data", "derived_data", 
                                    "jaccard_temporal_turnover_5km.rds"))
saveRDS(temporal_turnover_15km, here("data", "derived_data", 
                                     "jaccard_temporal_turnover_15km.rds"))

# END OF SCRIPT ----------------------------------------------------------------