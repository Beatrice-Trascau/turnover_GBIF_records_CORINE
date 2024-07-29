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
    year %in% c(1997:2000) ~ "before_1",
    year %in% c(2006:2009) ~ "after_1",
    year %in% c(2003:2006) ~ "before_2",
    year %in% c(2012:2015) ~ "after_2",
    year %in% c(2009:2012) ~ "before_3",
    year %in% c(2015:2018) ~ "after_3",
    TRUE ~ NA_character_
  )) %>% 
  filter(!is.na(period))

## 2.3. Assign species to land cover cells -------------------------------------

# Extract land cover cell ID for each occurrence record
clean_occurrences_for_turnover <- clean_occurrences_sf |>
  mutate(cell = terra::extract(land_cover_id, 
                               as.matrix(st_coordinates(clean_occurrences_sf))))

# Filter cells with more than 3 species
species_count <- clean_occurrences_for_turnover |>
  st_drop_geometry() |>
  group_by(cell, period) |>
  summarise(species_count = n_distinct(species))

# Get cells with more than 3 species 
cells_with_more_than_3_species <- species_count |>
  filter(species_count > 3) |>
  pull(cell)

# Filter original data for the specific cells
filtered_occurrences <- clean_occurrences_for_turnover |>
  filter(cell %in% cells_with_more_than_3_species) |>
  st_drop_geometry() |>
  select(period, species, cell) |>
  na.omit()


# Convert to long format
occurrences_long <- clean_occurrences_for_turnover |>
  st_drop_geometry() |>
  select(period, species, cell) %>%
  na.omit()

# 3. CALCULATE TURNOVER --------------------------------------------------------