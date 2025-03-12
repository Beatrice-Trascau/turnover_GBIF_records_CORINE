##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 3.1_turnover_land_cover
# This script contains code which calculates the temporal turnover of GBIF
# records in CORINE land cover pixels
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

## 1.1. 15km CORINE Status Layers ----------------------------------------------

# Forest -> TWS
forest_tws_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_forest_tws_masked.tif"))

# TWS -> Forest
tws_forest_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_tws_forest_masked.tif"))

# All -> Urban
all_urban_15km <- rast(here("data", "derived_data", 
                            "clc_status_15km_all_urban_combined_masked.tif"))


## 1.2. Cleaned GBIF Occurrences -----------------------------------------------

# Cleaned occurrence records
occurrences_norway <- fread(here("data", "derived_data", 
                                 "clean_occurrences_15km.txt"))

# 2. SPATIAL REFERENCE GRID ----------------------------------------------------

## 2.1. Identify non-NA cells across all raster stacks -------------------------

# Create a reference grid raster with same properties as CLC
reference_grid_15km <- forest_tws_15km[[1]]

# Create a cell validity mask with all values = FALSE
valid_cell_mask <- rep(FALSE, ncell(reference_grid_15km))

# Check all layers in Forest -> TWS raster
for(i in 1:nlyr(forest_tws_15km)){
  valid_cell_mask <- valid_cell_mask | !is.na(values(forest_tws_15km[[i]]))
}

# Check all layers in TWs -> Forest raster
for(i in 1:nlyr(tws_forest_15km)){
  valid_cell_mask <- valid_cell_mask | !is.na(values(tws_forest_15km[[i]]))
}

# Check all layers in TWs -> Forest raster
for(i in 1:nlyr(all_urban_15km)){
  valid_cell_mask <- valid_cell_mask | !is.na(values(all_urban_15km[[i]]))
}

# Get indices of cells with valid data in at least one layer
valid_cells <- which(valid_cell_mask)

## 2.2. Create reference grid with consecutive cell IDs ------------------------

# Create the grid with all cells = NA
reference_grid_15km[] <- NA

# Assign consecutive cell IDs
reference_grid_15km[valid_cells] <- 1:length(valid_cells)

# Convert to dataframe with coordinates and cell IDs
grid_15km_df <- as.data.frame(reference_grid_15km, xy = TRUE) |>
  rename(cell_id = '2000-2006_Non-forest') |>
  na.omit()

# Convert df to sf object
grid_15km_sf <- st_as_sf(grid_15km_df, coords = c("x", "y"),
                         crs = st_crs(forest_tws_15km))

## 2.3. Validate spatial coverage of the grid ----------------------------------

# Count cells that are NA across all layers in all raster stacks
all_na_count <- sum(apply(is.na(values(forest_tws_15km)), 1, all)&
                    apply(is.na(values(tws_forest_15km)), 1, all) &
                    apply(is.na(values(all_urban_15km)), 1, all))

# Validate spatial coverage
cat("Total raster cells in grid extent:", ncell(reference_grid_15km), "\n")
cat("Cells with valid data in at least one layer:", nrow(grid_15km_df), "\n")
cat("Cells without valid data in any layer:", all_na_count, "\n")
cat("Verification check:", ncell(reference_grid_15km) - nrow(grid_15km_df) == all_na_count, "\n")
cat("Proportion of grid with analyzable data:", 
    round(nrow(grid_15km_df)/ncell(reference_grid_15km)*100, 2), "%\n")

# 3. EXTRACT LAND COVER CHANGE VALUES FOR EACH CELL ----------------------------

## 3.1. Extract all land cover values ------------------------------------------

# Create df for the land cover values for each cell
lc_values <- data.frame(cell_id = grid_15km_df$cell_id)

# Extract land cover values for Forest -> TWS
for(layer_name in names(forest_tws_15km)){
  layer_values <- terra::extract(forest_tws_15km[[layer_name]], 
                          grid_15km_df[, c("x", "y")])
  lc_values[[layer_name]] <- layer_values[, 2]
}

# Extract land cover values for TWS -> Forest
for(layer_name in names(tws_forest_15km)){
  layer_values <- terra::extract(tws_forest_15km[[layer_name]], 
                                 grid_15km_df[, c("x", "y")])
  lc_values[[layer_name]] <- layer_values[, 2]
}

# Extract land cover values for All -> Urban
for(layer_name in names(all_urban_15km)){
  layer_values <- terra::extract(all_urban_15km[[layer_name]], 
                                 grid_15km_df[, c("x", "y")])
  lc_values[[layer_name]] <- layer_values[, 2]
}

# Check if all layers were exracted
ncol(lc_values)-1 #33
nlyr(forest_tws_15km) #12
nlyr(tws_forest_15km) #12
nlyr(all_urban_15km) #9

# 4. PREPARE OCCURRENCES FOR ANALYSIS ------------------------------------------

## 4.1. Convert occurrences to spatial data ------------------------------------

# Convert occurrence records to sf 
occurrences_sf <- st_as_sf(occurrences_norway,
                           coords = c("decimalLongitude", "decimalLatitude"),
                           crs = 4326)

# Re-project to match CORINE CRS
occurrences_sf_reprojected <- st_transform(occurrences_sf,
                                           st_crs(forest_tws_15km))


## 4.2. Assign before and after period for each time period --------------------

# Define time periods (2000-2006, 2006-2012, 2012-2018)
occurrences_with_periods <- occurrences_sf_reprojected |>
  mutate(period_2000_2006 = case_when(year %in% 1997:2000 ~ "1997-2000", #Before
                                      year %in% 2006:2009 ~ "2006-2009", # After
                                      TRUE ~ NA_character_),
         period_2006_2012 = case_when(year %in% 2003:2006 ~ "2003-2006", #Before
                                      year %in% 2012:2015 ~ "2012-2015", # After
                                      TRUE ~ NA_character_),
         period_2012_2018 = case_when(year %in% 2008:2012 ~ "2008-2012", #Before
                                      year %in% 2015:2018 ~ "2015-2018", # After
                                      TRUE ~ NA_character_)) |>
  # Keep only rows with at least one period assigned
  filter(!is.na(period_2000_2006) | !is.na(period_2006_2012) | !is.na(period_2012_2018))

## 4.3. Assign grid cell IDs to occurrences ------------------------------------

# Create a 1m buffer for grid cells to account for boundary precision issues

