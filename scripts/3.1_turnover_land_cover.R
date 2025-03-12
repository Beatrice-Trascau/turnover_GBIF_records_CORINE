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
valud_cells <- which(valid_cell_mask)

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
