##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 3.2_turnover_aggregated_land_cover
# This script contains code which calculates the temporal turnover of GBIF
# records in the aggregate CORINE land cover layers (1km, 5km and 15km)
##----------------------------------------------------------------------------##

# # 1. LOAD DATA -----------------------------------------------------------------

## 1.1. Download files if neccessary -------------------------------------------

# SSB Grid
# download_files("https://ntnu.box.com/shared/static/pjb2qr9aptugu7awlqmomx2o9d40sdhp.zip", 
#                "data/raw_data/norway_corine_change_modified_stack.tif")

## 1.2. Read in data -----------------------------------------------------------

# CORINE Aggregated Layers
corine_1km <- rast(here("data", "derived_data", "corine_1km.tif"))
corine_5km <- rast(here("data", "derived_data", "corine_5km.tif"))
corine_15km <- rast(here("data", "derived_data", "corine_15km.tif"))

# SSB Grid
ssb_grids <- vect(here("data", "raw_data",
                       "SSB050KM", "ssb50km.shp"))

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

