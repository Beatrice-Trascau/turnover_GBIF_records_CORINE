##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 3.1_turnover_land_cover
# This script contains code which calculates the temporal turnover of GBIF
# records in CORINE land cover pixels
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

## 1.1. Download files if neccessary -------------------------------------------

# CORINE Status Layers
download_files("https://ntnu.box.com/shared/static/97g9x4839ij4lnlldji2wh8e0e2lm5bf.tif", 
              "data/derived_data/norway_corine_change_modified_stack.tif")

# SSB Grid
download_files("https://ntnu.box.com/shared/static/pjb2qr9aptugu7awlqmomx2o9d40sdhp.zip", 
              "data/raw_data/norway_corine_change_modified_stack.tif")

## 1.2. Read in data -----------------------------------------------------------

# CORINE Status Layers
norway_corine_status_modified_stack <- rast(here("data", "derived_data",
                                                 "norway_corine_status_modified_stack.tif"))

# SSB Grid
ssb_grids <- vect(here("data", "raw_data",
                       "SSB050KM", "ssb50km.shp"))

# Cleaned occurrence records
occurrences_norway <- fread(here("data", "derived_data", 
                                 "cleaned_occurrences_july24.txt"))

# 2. PREPARE DATA FOR ANALYSIS -------------------------------------------------

## 2.1. Raster with unique ID for each cell ------------------------------------

# Create an ID raster with the same properties as CORINE
land_cover_id <- norway_corine_status_modified_stack[[1]]

#Assign each cell a unique ID from 1 to ncell
land_cover_id[] <- 1:ncell(norway_corine_status_modified_stack[[1]])

## 2.2. Occurrence records as sf -----------------------------------------------

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

## 2.3. Assign species to land cover cells -------------------------------------

# Extract coordinates from the occurrence records
coords <- st_coordinates(clean_occurrences_sf)

# Extract land cover cell ID for each occurrence record
cell_ids <- terra::extract(land_cover_id, coords)[, 1]

# Ensure we have the correct number of cell IDs
stopifnot(length(cell_ids) == nrow(clean_occurrences_sf))

# Add cell_ids in occurrence df in new column
clean_occurrences_for_turnover <- clean_occurrences_sf |>
  mutate(cell = cell_ids)

# Find the cells that have more than 3 species
species_count <- clean_occurrences_for_turnover |>
  st_drop_geometry() |>
  group_by(cell, period) |>
  summarise(species_count = n_distinct(species), .groups = 'drop')

# Get the cells that have more than 3 species
cells_with_more_than_3_species <- species_count |>
  filter(species_count > 3) |>
  pull(cell)

# Filter original data for the specific cells and convert to long format
filtered_occurrences <- clean_occurrences_for_turnover |>
  filter(cell %in% cells_with_more_than_3_species) |>
  st_drop_geometry() |>
  select(period, species, cell) |>
  na.omit()

# Aggregate to ensure unique species records per cell and period
unique_occurrences <- filtered_occurrences |>
  distinct(cell, species, period, .keep_all = TRUE)

# 3. CALUCATE TEMPORAL TURNOVER FOR EACH PAIR OF CHANGE PERIODS ----------------

# Define pairs of periods
period_combinations <- list(
  c("1997-2000", "2006-2009"),
  c("2003-2006", "2012-2015"),
  c("2009-2012", "2015-2018")
)

# Apply (custom) function to calculate Jaccard's dissimilarity index between the
# two periods in a pair
jaccard_index_list <- lapply(period_combinations, function(periods) {
  calculate_jaccard_for_periods(unique_occurrences, periods[1], periods[2])})

# Combine all results into one df
jaccard_results_combined <- bind_rows(jaccard_results_list)

# Write df to file
saveRDS(jaccard_results_combined, here("data", "derived_data",
                                       "jaccard_temporal_turnover_with_land_cover.rds"))

