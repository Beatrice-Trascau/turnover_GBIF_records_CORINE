##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 3.1_turnover_plants_land_cover
# This script contains code which calculates the temporal turnover of GBIF
# vascular plant records in CORINE land cover pixels
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

# 2. CREATE A SPATIAL REFERENCE GRID WITH CELL IDS -----------------------------

# Create a reference grid from the CLC rasters
reference_grid <- forest_tws_15km[[1]]

# Reset all values to NA
reference_grid[] <-NA

# Find all valid cells across all raster layers
valid_cells <- which(!is.na(values(forest_tws_15km[[1]])) |
                       !is.na(values(tws_forest_15km[[1]])) |
                       !is.na(values(all_urban_15km[[1]])))

# Assign sequential cell IDs to valid cells
reference_grid[valid_cells] <- 1:length(valid_cells)

# Rename layer to cell_id
names(reference_grid) <- "cell_id"

# Create dataframe with cell IDs and coordinates
grid_df <- as.data.frame(reference_grid, xy = TRUE) |>
  filter(!is.na(cell_id))

# 3. PREPARE OCCURRENCE DATA ---------------------------------------------------

## 3.1. Convert occurrences to sf ----------------------------------------------

# Filter out everything except for the vascular plants
vascular_plant_occurrences <- occurrences_norway |>
  filter(phylum == "Tracheophyta")

# Convert occurrences to sf and reproject to match CLC
vascular_plant_occurrences_sf <- st_as_sf(vascular_plant_occurrences,
                           coords = c("decimalLongitude", "decimalLatitude"),
                           crs = 4326) |>
  st_transform(crs(reference_grid))

## 3.2. Filter for the before and after periods for 2000-2006 ------------------

# Filter occurrences and label the periods
plants_sf_2000_2006 <- vascular_plant_occurrences_sf |>
  filter(year %in% c(1997:2000, 2006:2009)) |>
  mutate(time_period = case_when(year %in% 1997:2000 ~ "1997-2000",
                                 year %in% 2006:2009 ~ "2006-2009",
                                 TRUE ~ NA_character_))

## 3.3. Filter for the before and after periods for 2006-2012 ------------------

# Filter occurrences and label the periods
plants_sf_2006_2012 <- vascular_plant_occurrences_sf |>
  filter(year %in% c(2003:2006, 2012:2015)) |>
  mutate(time_period = case_when(year %in% 2003:2006 ~ "2003-2006",
                                 year %in% 2012:2015 ~ "2012-2015",
                                 TRUE ~ NA_character_))

## 3.4. Filter for the before and after periods for 2012-2018 ------------------

# Filter occurrences and label the periods
plants_sf_2012_2018 <- vascular_plant_occurrences_sf |>
  filter(year %in% c(2009:2012, 2015:2018)) |>
  mutate(time_period = case_when(year %in% 2009:2012 ~ "2009-2012",
                                 year %in% 2015:2018 ~ "2015-2018",
                                 TRUE ~ NA_character_))

# 4. EXTRACT CELL ID -----------------------------------------------------------

## 4.1. First period: 2000-2006 ------------------------------------------------

# Extract cell IDs from the reference grid
plants_sf_2000_2006$cell_ID <- terra::extract(reference_grid,
                                              st_coordinates(plants_sf_2000_2006))[, "cell_id"]

# Remove occurrences that fall outside valid grid cells
plants_sf_2000_2006_valid <- plants_sf_2000_2006 |>
  filter(!is.na(cell_ID))

## 4.2. Second period: 2006-2012 -----------------------------------------------

# Extract cell IDs from the reference grid
plants_sf_2006_2012$cell_ID <- terra::extract(reference_grid,
                                              st_coordinates(plants_sf_2006_2012))[, "cell_id"]

# Remove occurrences that fall outside valid grid cells
plants_sf_2006_2012_valid <- plants_sf_2006_2012 |>
  filter(!is.na(cell_ID))

## 4.3. Third period: 2012-2018 ------------------------------------------------

# Extract cell IDs from the reference grid
plants_sf_2012_2018$cell_ID <- terra::extract(reference_grid,
                                                   st_coordinates(plants_sf_2012_2018))[, "cell_id"]

# Remove occurrences that fall outside valid grid cells
plants_sf_2012_2018_valid <- plants_sf_2012_2018 |>
  filter(!is.na(cell_ID))
