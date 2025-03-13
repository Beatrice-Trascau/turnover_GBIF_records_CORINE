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

# Convert occurrences to sf and reproject to match CLC
occurrences_sf <- st_as_sf(occurrences_norway,
                           coords = c("decimalLongitude", "decimalLatitude"),
                           crs = 4326) |>
  st_transform(crs(reference_grid))
 
## 3.2. Filter for the before and after periods for 2000-2006 ------------------

# Filter occurrences and label the periods
occurrences_sf_2000_2006 <- occurrences_sf |>
  filter(year %in% c(1997:2000, 2006:2009)) |>
  mutate(time_period = case_when(year %in% 1997:2000 ~ "1997-2000",
                                 year %in% 2006:2009 ~ "2006-2009",
                                 TRUE ~ NA_character_))

# Check filtering results
cat("Period 1 (2000-2006) occurrences:\n")
cat("Total records:", nrow(occurrences_sf_2000_2006), "\n")
cat("Before period (1997-2000):",
    sum(occurrences_sf_2000_2006$time_period == "1997-2000"), "records\n")
cat("After period (2006-2009):", 
    sum(occurrences_sf_2000_2006$time_period == "2006-2009"), "records\n")

## 3.3. Filter for the before and after periods for 2006-2012 ------------------

# Filter occurrences and label the periods
occurrences_sf_2006_2012 <- occurrences_sf |>
  filter(year %in% c(2003:2006, 2012:2015)) |>
  mutate(time_period = case_when(year %in% 2003:2006 ~ "2003-2006",
                                 year %in% 2012:2015 ~ "2012-2015",
                                 TRUE ~ NA_character_))

# Check filtering results
cat("Period 2 (2006-2012) occurrences:\n")
cat("Total records:", nrow(occurrences_sf_2006_2012), "\n")
cat("Before period (2003-2006):",
    sum(occurrences_sf_2006_2012$time_period == "2003-2006"), "records\n")
cat("After period (2012-2015):", 
    sum(occurrences_sf_2006_2012$time_period == "2012-2015"), "records\n")

# 4. EXTRACT CELL ID -----------------------------------------------------------

## 4.1. First period: 2000-2006 ------------------------------------------------

# Extract cell IDs from the reference grid
occurrences_sf_2000_2006$cell_ID <- terra::extract(reference_grid,
                                                   st_coordinates(occurrences_sf_2000_2006))[, "cell_id"]

# Remove occurrences that fall outside valid grid cells
occurrences_sf_2000_2006_valid <- occurrences_sf_2000_2006 |>
  filter(!is.na(cell_ID))

# Calculate proportion of records retained with valid cells
rentention_rate_2000_2006 <- nrow(occurrences_sf_2000_2006_valid) / nrow(occurrences_sf_2000_2006) * 100

# Summarize spatial assignment results
cat("Cell ID extraction results:\n")
cat("Total occurrences processed:", nrow(occurrences_sf_2000_2006), "\n")
cat("Occurrences assigned to valid grid cells:", nrow(occurrences_sf_2000_2006_valid), 
    sprintf("(%.2f%%)\n", rentention_rate_2000_2006))
cat("Number of unique cells with occurrences:", 
    length(unique(occurrences_sf_2000_2006_valid$cell_ID)), "\n")

## 4.2. Second period: 2006-2012 -----------------------------------------------

# Extract cell IDs from the reference grid
occurrences_sf_2006_2012$cell_ID <- terra::extract(reference_grid,
                                                   st_coordinates(occurrences_sf_2006_2012))[, "cell_id"]

# Remove occurrences that fall outside valid grid cells
occurrences_sf_2006_2012_valid <- occurrences_sf_2006_2012 |>
  filter(!is.na(cell_ID))

# Calculate proportion of records retained with valid cells
rentention_rate_2006_2012 <- nrow(occurrences_sf_2006_2012_valid) / nrow(occurrences_sf_2006_2012) * 100

# Summarize spatial assignment results
cat("Cell ID extraction results:\n")
cat("Total occurrences processed:", nrow(occurrences_sf_2006_2012), "\n")
cat("Occurrences assigned to valid grid cells:", nrow(occurrences_sf_2006_2012_valid), 
    sprintf("(%.2f%%)\n", rentention_rate_2006_2012))
cat("Number of unique cells with occurrences:", 
    length(unique(occurrences_sf_2006_2012_valid$cell_ID)), "\n")


# 5. GET SPECIES LISTS AND NUMBER OF OCCURRENCES FOR EACH CELL -----------------

## 5.1. First period: 2000-2006 ------------------------------------------------

# Create species list and total number of occurrences for each before and after period
cell_2000_2006_summary <- occurrences_sf_2000_2006_valid |>
  st_drop_geometry() |>
  group_by(cell_ID, time_period) |>
  summarize(species_list = list(unique(species)),
            n_species = length(unique(species)),
            n_occurrences = n(),
            .group = "drop")

# Reshape to wide formate to separate before and after columns
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
cell_2006_2012_summary <- occurrences_sf_2006_2012_valid |>
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

# Save df 
save(turnover_2000_2006_lc,
     file = here("data", "derived_data", "turnover_2000_2006_all_data.rda"))

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

# Save df 
save(turnover_2006_2012_lc,
     file = here("data", "derived_data", "turnover_2006_2012_all_data.rda"))
