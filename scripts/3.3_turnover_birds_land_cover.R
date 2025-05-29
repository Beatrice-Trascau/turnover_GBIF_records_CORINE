##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 3.3_turnover_birds_land_cover
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

## 1.3. SSB ID Grid ------------------------------------------------------------

# SSB ID Grid
ssb_grid <- st_read(here("data", "raw_data", "SSB050KM", "ssb50km.shp"))

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
bird_occurrences <- occurrences_norway |>
  filter(class == "Aves")

# Convert occurrences to sf and reproject to match CLC
bird_occurrences_sf <- st_as_sf(bird_occurrences,
                                          coords = c("decimalLongitude", "decimalLatitude"),
                                          crs = 4326) |>
  st_transform(crs(reference_grid))

## 3.2. Filter for the before and after periods for 2000-2006 ------------------

# Filter occurrences and label the periods
birds_sf_2000_2006 <- bird_occurrences_sf |>
  filter(year %in% c(1997:2000, 2006:2009)) |>
  mutate(time_period = case_when(year %in% 1997:2000 ~ "1997-2000",
                                 year %in% 2006:2009 ~ "2006-2009",
                                 TRUE ~ NA_character_))

## 3.3. Filter for the before and after periods for 2006-2012 ------------------

# Filter occurrences and label the periods
birds_sf_2006_2012 <- bird_occurrences_sf |>
  filter(year %in% c(2003:2006, 2012:2015)) |>
  mutate(time_period = case_when(year %in% 2003:2006 ~ "2003-2006",
                                 year %in% 2012:2015 ~ "2012-2015",
                                 TRUE ~ NA_character_))

## 3.4. Filter for the before and after periods for 2012-2018 ------------------

# Filter occurrences and label the periods
birds_sf_2012_2018 <- bird_occurrences_sf |>
  filter(year %in% c(2009:2012, 2015:2018)) |>
  mutate(time_period = case_when(year %in% 2009:2012 ~ "2009-2012",
                                 year %in% 2015:2018 ~ "2015-2018",
                                 TRUE ~ NA_character_))

# 4. EXTRACT CELL ID -----------------------------------------------------------

## 4.1. First period: 2000-2006 ------------------------------------------------

# Extract cell IDs from the reference grid
birds_sf_2000_2006$cell_ID <- terra::extract(reference_grid,
                                              st_coordinates(birds_sf_2000_2006))[, "cell_id"]

# Remove occurrences that fall outside valid grid cells
birds_sf_2000_2006_valid <- birds_sf_2000_2006 |>
  filter(!is.na(cell_ID))

## 4.2. Second period: 2006-2012 -----------------------------------------------

# Extract cell IDs from the reference grid
birds_sf_2006_2012$cell_ID <- terra::extract(reference_grid,
                                              st_coordinates(birds_sf_2006_2012))[, "cell_id"]

# Remove occurrences that fall outside valid grid cells
birds_sf_2006_2012_valid <- birds_sf_2006_2012 |>
  filter(!is.na(cell_ID))

## 4.3. Third period: 2012-2018 ------------------------------------------------

# Extract cell IDs from the reference grid
birds_sf_2012_2018$cell_ID <- terra::extract(reference_grid,
                                              st_coordinates(birds_sf_2012_2018))[, "cell_id"]

# Remove occurrences that fall outside valid grid cells
birds_sf_2012_2018_valid <- birds_sf_2012_2018 |>
  filter(!is.na(cell_ID))

# 5. GET SPECIES LISTS AND NUMBER OF OCCURRENCES FOR EACH CELL -----------------

## 5.1. First period: 2000-2006 ------------------------------------------------

# Create species list and total number of occurrences for each before and after period
cell_2000_2006_summary <- birds_sf_2000_2006_valid |>
  st_drop_geometry() |>
  group_by(cell_ID, time_period) |>
  summarize(species_list = list(unique(species)),
            n_species = length(unique(species)),
            n_occurrences = n(),
            .group = "drop")

# Reshape to wide format to separate before and after columns
cell_2000_2006_summary_wide <- cell_2000_2006_summary |>
  # convert to wide format
  pivot_wider(id_cols = cell_ID,
              names_from = time_period,
              values_from = c(species_list, n_species, n_occurrences)) |>
  # rename columns 
  rename(species_list_before = 'species_list_1997-2000',
         species_list_after = 'species_list_2006-2009',
         total_spp_before = 'n_species_1997-2000',
         total_spp_after = 'n_species_2006-2009',
         total_occ_before = 'n_occurrences_1997-2000',
         total_occ_after = 'n_occurrences_2006-2009') |>
  # replace NA values in occurrence columns with 0
  mutate(total_spp_before = ifelse(is.na(total_spp_before), 0, total_spp_before),
         total_spp_after = ifelse(is.na(total_spp_after), 0, total_spp_after),
         total_occ_before = ifelse(is.na(total_occ_before), 0, total_occ_before),
         total_occ_after = ifelse(is.na(total_occ_after), 0, total_occ_after)) |>
  # add time period, recorder effort and delta recorder effort columns
  mutate(lc_time_period = "2000-2006",
         recorder_effort = total_occ_before + total_occ_after,
         delta_recorder_effort = abs((total_occ_before - total_occ_after)/(total_occ_before + total_occ_after)))

# Replace NULL values in the species_list_columns
cell_2000_2006_summary_wide <- cell_2000_2006_summary_wide |>
  mutate(species_list_before = lapply(species_list_before, 
                                      function(x) if(is.null(x)) character(0) else x),
         species_list_after = lapply(species_list_after, 
                                     function(x) if(is.null(x)) character(0) else x))

## 5.2. Second period: 2006-2012 -----------------------------------------------

# Create species list and total number of occurrences for each before and after period
cell_2006_2012_summary <- birds_sf_2006_2012_valid |>
  st_drop_geometry() |>
  group_by(cell_ID, time_period) |>
  summarize(species_list = list(unique(species)),
            n_species = length(unique(species)),
            n_occurrences = n(),
            .group = "drop")

# Reshape to wide formate to separate before and after columns
cell_2006_2012_summary_wide <- cell_2006_2012_summary |>
  # convert to wide format
  pivot_wider(id_cols = cell_ID,
              names_from = time_period,
              values_from = c(species_list, n_species, n_occurrences)) |>
  # rename columns 
  rename(species_list_before = 'species_list_2003-2006',
         species_list_after = 'species_list_2012-2015',
         total_spp_before = 'n_species_2003-2006',
         total_spp_after = 'n_species_2012-2015',
         total_occ_before = 'n_occurrences_2003-2006',
         total_occ_after = 'n_occurrences_2012-2015') |>
  # replace NA values in occurrence columns with 0
  mutate(total_spp_before = ifelse(is.na(total_spp_before), 0, total_spp_before),
         total_spp_after = ifelse(is.na(total_spp_after), 0, total_spp_after),
         total_occ_before = ifelse(is.na(total_occ_before), 0, total_occ_before),
         total_occ_after = ifelse(is.na(total_occ_after), 0, total_occ_after)) |>
  # add time period, recorder effort and delta recorder effort columns
  mutate(lc_time_period = "2006-2012",
         recorder_effort = total_occ_before + total_occ_after,
         delta_recorder_effort = abs((total_occ_before - total_occ_after)/(total_occ_before + total_occ_after)))

# Replace NULL values in the species_list_columns
cell_2006_2012_summary_wide <- cell_2006_2012_summary_wide |>
  mutate(species_list_before = lapply(species_list_before, 
                                      function(x) if(is.null(x)) character(0) else x),
         species_list_after = lapply(species_list_after, 
                                     function(x) if(is.null(x)) character(0) else x))

## 5.3. Third period: 2012-2018 ------------------------------------------------

# Create species list and total number of occurrences for each before and after period
cell_2012_2018_summary <- birds_sf_2012_2018_valid |>
  st_drop_geometry() |>
  group_by(cell_ID, time_period) |>
  summarize(species_list = list(unique(species)),
            n_species = length(unique(species)),
            n_occurrences = n(),
            .group = "drop")

# Reshape to wide formate to separate before and after columns
cell_2012_2018_summary_wide <- cell_2012_2018_summary |>
  # convert to wide format
  pivot_wider(id_cols = cell_ID,
              names_from = time_period,
              values_from = c(species_list, n_species, n_occurrences)) |>
  # rename columns 
  rename(species_list_before = 'species_list_2009-2012',
         species_list_after = 'species_list_2015-2018',
         total_spp_before = 'n_species_2009-2012',
         total_spp_after = 'n_species_2015-2018',
         total_occ_before = 'n_occurrences_2009-2012',
         total_occ_after = 'n_occurrences_2015-2018') |>
  # replace NA values in occurrence columns with 0
  mutate(total_spp_before = ifelse(is.na(total_spp_before), 0, total_spp_before),
         total_spp_after = ifelse(is.na(total_spp_after), 0, total_spp_after),
         total_occ_before = ifelse(is.na(total_occ_before), 0, total_occ_before),
         total_occ_after = ifelse(is.na(total_occ_after), 0, total_occ_after)) |>
  # add time period, recorder effort and delta recorder effort columns
  mutate(lc_time_period = "2012-2018",
         recorder_effort = total_occ_before + total_occ_after,
         delta_recorder_effort = abs((total_occ_before - total_occ_after)/(total_occ_before + total_occ_after)))

# Replace NULL values in the species_list_columns
cell_2012_2018_summary_wide <- cell_2012_2018_summary_wide |>
  mutate(species_list_before = lapply(species_list_before, 
                                      function(x) if(is.null(x)) character(0) else x),
         species_list_after = lapply(species_list_after, 
                                     function(x) if(is.null(x)) character(0) else x))

# 6. GET LC VALUES -------------------------------------------------------------

## 6.1. First period: 2000-2006 ------------------------------------------------

# Check which layers are needed
names(forest_tws_15km) #[[2]] = "2000-2006_Forest no change", [[3]] = "2000-2006_Forest to TWS"
names(tws_forest_15km) #[[2]] = "2000-2006_TWS no change", [[3]] = "2000-2006_TWS to Forest"
names(all_urban_15km) #[[1]] = "2000-2006_Urban_no_change", [[2]] = "2000-2006_all_to_urban"

### 6.1.1. Forests -> TWS ------------------------------------------------------

# Create df with raster layers of interest
forest_tws_df <- as.data.frame(forest_tws_15km[[c(2,3)]], xy = TRUE)

# Extract cell IDs for each coordinate pair
forest_tws_df$cell_ID <- terra::extract(reference_grid,
                                        forest_tws_df[, c("x", "y")])[, "cell_id"]

# Remove rows with NA cell_ID
forest_tws_df <- forest_tws_df |>
  filter(!is.na(cell_ID))

# Merge these values with the occurrences df
cell_2000_2006_summary_wide_LC1 <- cell_2000_2006_summary_wide |>
  left_join(forest_tws_df, by = "cell_ID")

### 6.1.2. TWS -> Forests ------------------------------------------------------

# Create df with raster values
tws_forest_df <- as.data.frame(tws_forest_15km[[c(2,3)]], xy = TRUE)

# Extract cell IDs for each coordinate pair
tws_forest_df$cell_ID <- terra::extract(reference_grid,
                                        tws_forest_df[, c("x", "y")])[, "cell_id"]

# Remove rows with NA cell_ID
tws_forest_df <- tws_forest_df |>
  filter(!is.na(cell_ID))

# Merge these values with the occurrences df
cell_2000_2006_summary_wide_LC2 <- cell_2000_2006_summary_wide_LC1 |>
  left_join(tws_forest_df |>
              select(cell_ID, '2000-2006_TWS no change', '2000-2006_TWS to Forest'),
            by = "cell_ID")

### 6.1.3. All -> Urban --------------------------------------------------------

# Create df with raster values
all_urban_df <- as.data.frame(all_urban_15km[[c(1,2)]], xy = TRUE)

# Extract cell IDs for each coordinate pair
all_urban_df$cell_ID <- terra::extract(reference_grid,
                                       all_urban_df[, c("x", "y")])[, "cell_id"]

# Remove rows with NA cell_ID
all_urban_df <- all_urban_df |>
  filter(!is.na(cell_ID))

# Merge these values with the occurrences df
cell_2000_2006_summary_wide_LC3 <- cell_2000_2006_summary_wide_LC2 |>
  left_join(all_urban_df |>
              select(cell_ID, '2000-2006_Urban_no_change', '2000-2006_all_to_urban'),
            by = "cell_ID")

## 6.2. Second period: 2006-2012 -----------------------------------------------

# Check which layers are needed
names(forest_tws_15km) #[[6]] = "2006-2012_Forest no change", [[7]] = "2006-2012_Forest to TWS"
names(tws_forest_15km) #[[6]] = "2006-2012_TWS no change", [[7]] = "2006-2012_TWS to Forest"
names(all_urban_15km) #[[4]] = "2006-2012_Urban_no_change", [[5]] = "2006-2012_all_to_urban"

### 6.2.1. Forests -> TWS ------------------------------------------------------

# Create df with raster layers of interest
forest_tws_2006_2012_df <- as.data.frame(forest_tws_15km[[c(6,7)]], xy = TRUE)

# Extract cell IDs for each coordinate pair
forest_tws_2006_2012_df$cell_ID <- terra::extract(reference_grid,
                                                  forest_tws_2006_2012_df[, c("x", "y")])[, "cell_id"]

# Remove rows with NA cell_ID
forest_tws_2006_2012_df <- forest_tws_2006_2012_df |>
  filter(!is.na(cell_ID))

# Merge these values with the occurrences df
cell_2006_2012_summary_wide_LC1 <- cell_2006_2012_summary_wide |>
  left_join(forest_tws_2006_2012_df, by = "cell_ID")

### 6.2.2. TWS -> Forests ------------------------------------------------------

# Create df with raster values
tws_forest_2006_2012_df <- as.data.frame(tws_forest_15km[[c(6,7)]], xy = TRUE)

# Extract cell IDs for each coordinate pair
tws_forest_2006_2012_df$cell_ID <- terra::extract(reference_grid,
                                                  tws_forest_2006_2012_df[, c("x", "y")])[, "cell_id"]

# Remove rows with NA cell_ID
tws_forest_2006_2012_df <- tws_forest_2006_2012_df |>
  filter(!is.na(cell_ID))

# Merge these values with the occurrences df
cell_2006_2012_summary_wide_LC2 <- cell_2006_2012_summary_wide_LC1 |>
  left_join(tws_forest_2006_2012_df |>
              select(cell_ID, '2006-2012_TWS no change', '2006-2012_TWS to Forest'),
            by = "cell_ID")

### 6.2.3. All -> Urban --------------------------------------------------------

# Create df with raster values
all_urban_2006_2012_df <- as.data.frame(all_urban_15km[[c(4,5)]], xy = TRUE)

# Extract cell IDs for each coordinate pair
all_urban_2006_2012_df$cell_ID <- terra::extract(reference_grid,
                                                 all_urban_2006_2012_df[, c("x", "y")])[, "cell_id"]

# Remove rows with NA cell_ID
all_urban_2006_2012_df <- all_urban_2006_2012_df |>
  filter(!is.na(cell_ID))

# Merge these values with the occurrences df
cell_2006_2012_summary_wide_LC3 <- cell_2006_2012_summary_wide_LC2 |>
  left_join(all_urban_2006_2012_df |>
              select(cell_ID, '2006-2012_Urban_no_change', '2006-2012_all_to_urban'),
            by = "cell_ID")


## 6.3. Third period: 2012-2018 ------------------------------------------------

# Check which layers are needed
names(forest_tws_15km) #[[10]] = "2012-2018_Forest no change", [[11]] = "2012-2018_Forest to TWS"
names(tws_forest_15km) #[[10]] = "2012-2018_TWS no change", [[11]] = "2012-2018_TWS to Forest"
names(all_urban_15km) #[[7]] = "2012-2018_Urban_no_change", [[8]] = "2012-2018_all_to_urban"


### 6.3.1. Forests -> TWS ------------------------------------------------------

# Create df with raster layers of interest
forest_tws_2012_2018_df <- as.data.frame(forest_tws_15km[[c(10,11)]], xy = TRUE)

# Extract cell IDs for each coordinate pair
forest_tws_2012_2018_df$cell_ID <- terra::extract(reference_grid,
                                                  forest_tws_2012_2018_df[, c("x", "y")])[, "cell_id"]

# Remove rows with NA cell_ID
forest_tws_2012_2018_df <- forest_tws_2012_2018_df |>
  filter(!is.na(cell_ID))

# Merge these values with the occurrences df
cell_2012_2018_summary_wide_LC1 <- cell_2012_2018_summary_wide |>
  left_join(forest_tws_2012_2018_df, by = "cell_ID")

### 6.3.2. TWS -> Forests ------------------------------------------------------

# Create df with raster values
tws_forest_2012_2018_df <- as.data.frame(tws_forest_15km[[c(10,11)]], xy = TRUE)

# Extract cell IDs for each coordinate pair
tws_forest_2012_2018_df$cell_ID <- terra::extract(reference_grid,
                                                  tws_forest_2012_2018_df[, c("x", "y")])[, "cell_id"]

# Remove rows with NA cell_ID
tws_forest_2012_2018_df <- tws_forest_2012_2018_df |>
  filter(!is.na(cell_ID))

# Merge these values with the occurrences df
cell_2012_2018_summary_wide_LC2 <- cell_2012_2018_summary_wide_LC1 |>
  left_join(tws_forest_2012_2018_df |>
              select(cell_ID, '2012-2018_TWS no change', '2012-2018_TWS to Forest'),
            by = "cell_ID")

### 6.3.3. All -> Urban --------------------------------------------------------

# Create df with raster values
all_urban_2012_2018_df <- as.data.frame(all_urban_15km[[c(7,8)]], xy = TRUE)

# Extract cell IDs for each coordinate pair
all_urban_2012_2018_df$cell_ID <- terra::extract(reference_grid,
                                                 all_urban_2012_2018_df[, c("x", "y")])[, "cell_id"]

# Remove rows with NA cell_ID
all_urban_2012_2018_df <- all_urban_2012_2018_df |>
  filter(!is.na(cell_ID))

# Merge these values with the occurrences df
cell_2012_2018_summary_wide_LC3 <- cell_2012_2018_summary_wide_LC2 |>
  left_join(all_urban_2012_2018_df |>
              select(cell_ID, '2012-2018_Urban_no_change', '2012-2018_all_to_urban'),
            by = "cell_ID")

# 7. CALCULATE JACCARD'S DISSIMILARITY INDEX -----------------------------------

## 7.1. First period: 2000-2006 ------------------------------------------------

# Filter cells with < 3 species in each time step
cell_2000_2006_filtered <- cell_2000_2006_summary_wide_LC3 |>
  filter(total_spp_before >= 3 & total_spp_after >= 3)

# Calculate Jaccard's Dissimilarity Index
turnover_2000_2006_lc <- cell_2000_2006_filtered |>
  rowwise() |>
  # calculate shared (intersection) and total (union) species
  mutate(intersection_size = length(intersect(species_list_before, species_list_after)),
         union_size = length(union(species_list_before, species_list_after)),
         JDI = 1 - (intersection_size / union_size)) |>
  ungroup()

## 7.2. Second period: 2006-2012 -----------------------------------------------

# Filter cells with < 3 species in each time step
cell_2006_2012_filtered <- cell_2006_2012_summary_wide_LC3 |>
  filter(total_spp_before >= 3 & total_spp_after >= 3)

# Calculate Jaccard's Dissimilarity Index
turnover_2006_2012_lc <- cell_2006_2012_filtered |>
  rowwise() |>
  # calculate shared (intersection) and total (union) species
  mutate(intersection_size = length(intersect(species_list_before, species_list_after)),
         union_size = length(union(species_list_before, species_list_after)),
         JDI = 1 - (intersection_size / union_size)) |>
  ungroup()

## 7.3. Third period: 2012-2018 ------------------------------------------------

# Filter cells with < 3 species in each time step
cell_2012_2018_filtered <- cell_2012_2018_summary_wide_LC3 |>
  filter(total_spp_before >= 3 & total_spp_after >= 3)

# Calculate Jaccard's Dissimilarity Index
turnover_2012_2018_lc <- cell_2012_2018_filtered |>
  rowwise() |>
  # calculate shared (intersection) and total (union) species
  mutate(intersection_size = length(intersect(species_list_before, species_list_after)),
         union_size = length(union(species_list_before, species_list_after)),
         JDI = 1 - (intersection_size / union_size)) |>
  ungroup()

## 7.4 Combine all into single df ----------------------------------------------

# Combine all periods into a single dataframe 
birds_all_periods_turnover_all_land_cover_chanegs_15km <- bind_rows(turnover_2000_2006_lc,
                                                                    turnover_2006_2012_lc,
                                                                    turnover_2012_2018_lc)

## 7.5. Add SSB ID Grid --------------------------------------------------------

# Check column names in SSB ID grid
names(ssb_grid)

# Transform grid to match reference grid CRS
ssb_grid_match_crs <- st_transform(ssb_grid, crs(reference_grid))

# Convert the grid df to sf object
grid_centroid <- st_as_sf(grid_df,
                          coords = c("x", "y"),
                          crs = crs(reference_grid))

# Intersect the ssb grid and centroids of cells
ssb_intersect <- st_intersection(grid_centroid, ssb_grid_match_crs)

# Check geometries
any(!st_is_valid(ssb_grid)) #FALSE

# Create lookup table
ssb_lookup <- data.frame(cell_ID = ssb_intersect$cell_id,
                         ssb_id  = ssb_intersect$SSBID)

# Check intersection was done correctly - only one SSB ID per cell ID
cell_id_counts <- table(ssb_intersect$cell_id)
cells_with_multiple_ssb <- cell_id_counts[cell_id_counts > 1]

if(length(cells_with_multiple_ssb) > 0) {
  cat("Warning: Some 15km cells intersect multiple SSB cells\n")
  print(cells_with_multiple_ssb)
} else {
  cat("Good: Each 15km cell maps to exactly one SSB cell\n")
}

# Add SSB ID to final dataframe
birds_all_periods_turnover_all_land_cover_chanegs_15km <- birds_all_periods_turnover_all_land_cover_chanegs_15km |>
  left_join(ssb_lookup, by = "cell_ID")

# Save df to file
save(birds_all_periods_turnover_all_land_cover_chanegs_15km,
     file = here("data", "derived_data", 
                 "birds_all_periods_turnover_all_land_cover_chanegs_15km.rda"))

# END OF SCRIPT ----------------------------------------------------------------