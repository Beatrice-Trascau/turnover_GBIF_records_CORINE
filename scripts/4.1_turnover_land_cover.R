##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 4.1_turnover_land_cover
# This script contains code which calculates the temporal turnover of GBIF
# records in CORINE land cover pixels
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

library(here)
source(here("scripts", "0_setup.R"))

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

# Keep only records with occurrenceStatus = PRESENT
occurrences_norway <- occurrences_norway |>
  filter(occurrenceStatus == "PRESENT")

## 1.3. SSB ID Grid ------------------------------------------------------------

# SSB ID Grid - 50km
ssb_grid_50km <- st_read(here("data", "raw_data", "SSB050KM", "ssb50km.shp"))

# SSB ID Grid - 250km
ssb_grid <- st_read(here("data", "raw_data", "SSB050KM", "ssb50km.shp"))

## 1.4. CHELSA processed data --------------------------------------------------

# Load temperature stack
chelsa_tas_stack <- rast(here("data", "derived_data", "chelsa", 
                              "chelsa_tas_period_means_norway.tif"))

# Load precipitation stack
chelsa_pr_stack <- rast(here("data", "derived_data", "chelsa", 
                             "chelsa_pr_period_means_norway.tif"))

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

## 3.4. Filter for the before and after periods for 2012-2018 ------------------

# Filter occurrences and label the periods
occurrences_sf_2012_2018 <- occurrences_sf |>
  filter(year %in% c(2009:2012, 2018:2021)) |>
  mutate(time_period = case_when(year %in% 2009:2012 ~ "2009-2012",
                                 year %in% 2018:2021 ~ "2018-2021",
                                 TRUE ~ NA_character_))

# Check filtering results
cat("Period 3 (2012-2018) occurrences:\n")
cat("Total records:", nrow(occurrences_sf_2012_2018), "\n")
cat("Before period (2009-2012):",
    sum(occurrences_sf_2012_2018$time_period == "2009-2012"), "records\n")
cat("After period (2018-2021):", 
    sum(occurrences_sf_2012_2018$time_period == "2018-2021"), "records\n")

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

## 4.3. Third period: 2012-2018 ------------------------------------------------

# Extract cell IDs from the reference grid
occurrences_sf_2012_2018$cell_ID <- terra::extract(reference_grid,
                                                   st_coordinates(occurrences_sf_2012_2018))[, "cell_id"]

# Remove occurrences that fall outside valid grid cells
occurrences_sf_2012_2018_valid <- occurrences_sf_2012_2018 |>
  filter(!is.na(cell_ID))

# Calculate proportion of records retained with valid cells
rentention_rate_2012_2018 <- nrow(occurrences_sf_2012_2018_valid) / nrow(occurrences_sf_2012_2018) * 100

# Summarize spatial assignment results
cat("Cell ID extraction results:\n")
cat("Total occurrences processed:", nrow(occurrences_sf_2012_2018), "\n")
cat("Occurrences assigned to valid grid cells:", nrow(occurrences_sf_2012_2018_valid), 
    sprintf("(%.2f%%)\n", rentention_rate_2012_2018))
cat("Number of unique cells with occurrences:", 
    length(unique(occurrences_sf_2012_2018_valid$cell_ID)), "\n")

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
         delta_recorder_effort = (total_occ_after - total_occ_before)/(total_occ_before + total_occ_after),
         abs_delta_recorder_effort = abs((total_occ_after - total_occ_before)/(total_occ_before + total_occ_after)))

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
         delta_recorder_effort = (total_occ_after - total_occ_before)/(total_occ_before + total_occ_after),
         abs_delta_recorder_effort = abs((total_occ_after - total_occ_before)/(total_occ_before + total_occ_after)))

# Replace NULL values in the species_list_columns
cell_2006_2012_summary_wide <- cell_2006_2012_summary_wide |>
  mutate(species_list_before = lapply(species_list_before, 
                                      function(x) if(is.null(x)) character(0) else x),
         species_list_after = lapply(species_list_after, 
                                     function(x) if(is.null(x)) character(0) else x))

## 5.3. Third period: 2012-2018 ------------------------------------------------

# Create species list and total number of occurrences for each before and after period
cell_2012_2018_summary <- occurrences_sf_2012_2018_valid |>
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
         species_list_after = 'species_list_2018-2021',
         total_spp_before = 'n_species_2009-2012',
         total_spp_after = 'n_species_2018-2021',
         total_occ_before = 'n_occurrences_2009-2012',
         total_occ_after = 'n_occurrences_2018-2021') |>
  # replace NA values in occurrence columns with 0
  mutate(total_spp_before = ifelse(is.na(total_spp_before), 0, total_spp_before),
         total_spp_after = ifelse(is.na(total_spp_after), 0, total_spp_after),
         total_occ_before = ifelse(is.na(total_occ_before), 0, total_occ_before),
         total_occ_after = ifelse(is.na(total_occ_after), 0, total_occ_after)) |>
  # add time period, recorder effort and delta recorder effort columns
  mutate(lc_time_period = "2012-2018",
         recorder_effort = total_occ_before + total_occ_after,
         delta_recorder_effort = (total_occ_after - total_occ_before)/(total_occ_before + total_occ_after),
         abs_delta_recorder_effort = abs((total_occ_after - total_occ_before)/(total_occ_before + total_occ_after)))

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

# 7. CALCULATE TEMPORAL TURNOVER -----------------------------------------------

## 7.0. Helper function to calculate betapart metrics --------------------------

# Function to calculate temporal turnover using betapart::beta.temp()
  # convers species lists to presence-absence matrices and calculates turnover
  # and nestedness components
calculate_betapart_temporal <- function(df, cell_id_col = "cell_ID"){
  
  # filter cells with >= 3 species in each time step
  df_filtered <- df |>
    filter(total_spp_before >= 3 & total_spp_after >= 3)
  
  # check how many cells were removed
  cat("  Cells retained after filtering (>= 3 species):", nrow(df_filtered), "\n")
  
  # get unique species across all cells in both time periods
  all_species <- df_filtered |>
    rowwise() |>
    summarise(species = list(c(species_list_before, species_list_after))) |>
    pull(species) |>
    unlist() |>
    unique()
  
  # check how many unique species there are
  cat("  Total unique species:", length(all_species), "\n")
  
  # create presence-absence matrices for before and after periods
  # each row is a cell, each column is a species (1 = present, 0 = absent)
  
  # sort all_species to ensure consistent ordering
  all_species <- sort(all_species)
  
  # before period matrix
  pa_before <- matrix(0, nrow = nrow(df_filtered), ncol = length(all_species))
  colnames(pa_before) <- all_species
  
  for(i in 1:nrow(df_filtered)) {
    species_present <- df_filtered$species_list_before[[i]]
    pa_before[i, all_species %in% species_present] <- 1
  }
  
  # after period matrix
  pa_after <- matrix(0, nrow = nrow(df_filtered), ncol = length(all_species))
  colnames(pa_after) <- all_species
  
  for(i in 1:nrow(df_filtered)) {
    species_present <- df_filtered$species_list_after[[i]]
    pa_after[i, all_species %in% species_present] <- 1
  }
  
  # verify matrices before passing to betapart
  cat("  Matrix dimensions - before:", paste(dim(pa_before), collapse=" x "), "\n")
  cat("  Matrix dimensions - after:", paste(dim(pa_after), collapse=" x "), "\n")
  cat("  Column names match:", identical(colnames(pa_before), colnames(pa_after)), "\n")
  cat("  All values 0 or 1 - before:", all(pa_before %in% c(0, 1)), "\n")
  cat("  All values 0 or 1 - after:", all(pa_after %in% c(0, 1)), "\n")
  
  # convert to dataframe for betapart
  # use check.names=FALSE to preserve original species names without modification
  pa_before <- as.data.frame(pa_before, check.names = FALSE)
  pa_after <- as.data.frame(pa_after, check.names = FALSE)
  
  # final check after conversion
  cat("  After df conversion - names match:", identical(names(pa_before), names(pa_after)), "\n")
  
  # calculate temporal beta diversity using betapart
  beta_temp_result <- betapart::beta.temp(pa_before, pa_after, 
                                          index.family = "jaccard")
  
  # add betapart results to the filtered dataframe
  # beta.temp returns a dataframe with one row per site (cell)
  df_filtered$beta_jtu <- beta_temp_result$beta.jtu    # Turnover component
  df_filtered$beta_jne <- beta_temp_result$beta.jne    # Nestedness component
  df_filtered$beta_jac <- beta_temp_result$beta.jac    # Total dissimilarity (jtu + jne)
  
  # Also calculate manual JDI for comparison/validation
  df_filtered <- df_filtered |>
    rowwise() |>
    mutate(
      intersection_size = length(intersect(species_list_before, species_list_after)),
      union_size = length(union(species_list_before, species_list_after)),
      JDI_manual = 1 - (intersection_size / union_size)
    ) |>
    ungroup()
  
  # update if function is finished
  cat("  ✓ Betapart calculations complete\n")
  
  return(df_filtered)
}

## 7.1. First period: 2000-2006 ------------------------------------------------

# Calculate temporal turnover for first period
turnover_2000_2006_lc <- calculate_betapart_temporal(cell_2000_2006_summary_wide_LC3)

# Save df 
save(turnover_2000_2006_lc,
     file = here("data", "derived_data", "turnover_2000_2006_all_data.rda"))

## 7.2. Second period: 2006-2012 -----------------------------------------------

# Calculate temporal turnover for second period
turnover_2006_2012_lc <- calculate_betapart_temporal(cell_2006_2012_summary_wide_LC3)

# Save df 
save(turnover_2006_2012_lc,
     file = here("data", "derived_data", "turnover_2006_2012_all_data.rda"))

## 7.3. Third period: 2012-2018 ------------------------------------------------

# Calculate temporal turnover for third period
turnover_2012_2018_lc <- calculate_betapart_temporal(cell_2012_2018_summary_wide_LC3)

# Save df 
save(turnover_2012_2018_lc,
     file = here("data", "derived_data", "turnover_2012_2018_all_data.rda"))

## 7.4. Validation: compare betapart & the manual JDI calculation --------------

# For each period, check if beta.jac matches manual JDI
for(period_name in c("2000-2006", "2006-2012", "2012-2018")) {
  
  if(period_name == "2000-2006") df <- turnover_2000_2006_lc
  if(period_name == "2006-2012") df <- turnover_2006_2012_lc
  if(period_name == "2012-2018") df <- turnover_2012_2018_lc
  
  cat("\nPeriod:", period_name, "\n")
  cat("  Correlation (beta.jac vs JDI_manual):", 
      round(cor(df$beta_jac, df$JDI_manual), 6), "\n")
  cat("  Mean absolute difference:", 
      round(mean(abs(df$beta_jac - df$JDI_manual)), 8), "\n")
  
  # Show first few values for inspection
  cat("  Sample comparison (first 5 cells):\n")
  print(df |> 
          select(cell_ID, beta_jac, JDI_manual, beta_jtu, beta_jne) |> 
          head(5))
}

## 7.5 Combine all into single df ----------------------------------------------

# Combine all periods into a single dataframe 
all_periods_turnover_all_land_cover_chanegs_15km <- bind_rows(turnover_2000_2006_lc,
                                                              turnover_2006_2012_lc,
                                                              turnover_2012_2018_lc)

## 7.6. Add SSB ID Grid --------------------------------------------------------

# Check column names in SSB ID grid
names(ssb_grid)

# Set crs of ssb grid to match that of the reference grid
if(is.na(st_crs(ssb_grid))) {
  st_crs(ssb_grid) <- crs(reference_grid)
  cat("CRS was missing - set to match reference grid CRS\n")
}

# Transform ssb grid to match reference grid CRS
ssb_grid_match_crs <- st_transform(ssb_grid, crs(reference_grid))

# Convert the ssb grid df to sf object
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
all_periods_turnover_all_land_cover_chanegs_15km <- all_periods_turnover_all_land_cover_chanegs_15km |>
  left_join(ssb_lookup, by = "cell_ID")

# Save df to file
# save(all_periods_turnover_all_land_cover_chanegs_15km,
#      file = here("data", "derived_data", 
#                  "all_periods_turnover_all_land_cover_chanegs_15km.rda"))

# 8. ADD CLIMATE DATA ----------------------------------------------------------

## 8.1. Prep climate data ------------------------------------------------------

# Check names ot the layers
print(names(chelsa_tas_stack))
print(names(chelsa_pr_stack))

## Define climate period mapping for the land cover periods
climate_period_mapping <- list("2000-2006" = list(before_climate = "1997-2000",
                                                  after_climate = "2006-2009",
                                                  before_precip = "1997-2000",
                                                  after_precip = "2006-2009"),
                               "2006-2012" = list(before_climate = "2003-2006",
                                                  after_climate = "2012-2015",
                                                  before_precip = "2003-2006", 
                                                  after_precip = "2012-2015"),
                               "2012-2018" = list(before_climate = "2009-2012",
                                                  after_climate = "2018-2021",
                                                  before_precip = "2009-2012",
                                                  after_precip = "2018-2021"))

## 8.2. Extract climate data for each period -----------------------------------

# Create list to store climate data for each period
climate_data_list <- list()

# Run through each period and extract climate data
for(lc_period in names(climate_period_mapping)){
  
  # track which period is being processed
  cat("\nProcessing climate data for LC period:", lc_period, "\n")
  
  # get climate layer names for the period
  climate_mapping <- climate_period_mapping[[lc_period]]
  
  # extract temperature layers (before and after)
  temp_before_layer <- climate_mapping$before_climate
  temp_after_layer <- climate_mapping$after_climate
  
  # extract precipitation layers (before and after)
  precip_before_layer <- climate_mapping$before_precip
  precip_after_layer <- climate_mapping$after_precip
  
  # check if the layers exist - give a warning if not
  if(!temp_before_layer %in% names(chelsa_tas_stack)) {
    cat("  Warning: Temperature layer", temp_before_layer, "not found\n")
    next
  }
  if(!temp_after_layer %in% names(chelsa_tas_stack)) {
    cat("  Warning: Temperature layer", temp_after_layer, "not found\n")
    next
  }
  if(!precip_before_layer %in% names(chelsa_pr_stack)) {
    cat("  Warning: Precipitation layer", precip_before_layer, "not found\n")
    next
  }
  if(!precip_after_layer %in% names(chelsa_pr_stack)) {
    cat("  Warning: Precipitation layer", precip_after_layer, "not found\n")
    next
  }
  
  # extract the specific layers
  temp_before <- chelsa_tas_stack[[temp_before_layer]]
  temp_after <- chelsa_tas_stack[[temp_after_layer]]
  precip_before <- chelsa_pr_stack[[precip_before_layer]]
  precip_after <- chelsa_pr_stack[[precip_after_layer]]
  
  # convert to dataframes with coordinates
  temp_before_df <- as.data.frame(temp_before, xy = TRUE)
  temp_after_df <- as.data.frame(temp_after, xy = TRUE)
  precip_before_df <- as.data.frame(precip_before, xy = TRUE)
  precip_after_df <- as.data.frame(precip_after, xy = TRUE)
  
  # extract cell IDs for each coordinate pair
  temp_before_df$cell_ID <- terra::extract(reference_grid, 
                                           temp_before_df[, c("x", "y")])[, "cell_id"]
  temp_after_df$cell_ID <- terra::extract(reference_grid,
                                          temp_after_df[, c("x", "y")])[, "cell_id"]
  precip_before_df$cell_ID <- terra::extract(reference_grid,
                                             precip_before_df[, c("x", "y")])[, "cell_id"]
  precip_after_df$cell_ID <- terra::extract(reference_grid,
                                            precip_after_df[, c("x", "y")])[, "cell_id"]
  
  # remove rows with NA cell_ID
  temp_before_df <- temp_before_df |> filter(!is.na(cell_ID))
  temp_after_df <- temp_after_df |> filter(!is.na(cell_ID))
  precip_before_df <- precip_before_df |> filter(!is.na(cell_ID))
  precip_after_df <- precip_after_df |> filter(!is.na(cell_ID))
  
  # rename columns so that they will make sense in the dataframe
  names(temp_before_df)[3] <- "temp_before"
  names(temp_after_df)[3] <- "temp_after"
  names(precip_before_df)[3] <- "precip_before"
  names(precip_after_df)[3] <- "precip_after"
  
  # merge all climate data for the period
  period_climate_data <- temp_before_df |>
    select(cell_ID, temp_before) |>
    left_join(temp_after_df |> select(cell_ID, temp_after), by = "cell_ID") |>
    left_join(precip_before_df |> select(cell_ID, precip_before), by = "cell_ID") |>
    left_join(precip_after_df |> select(cell_ID, precip_after), by = "cell_ID") |>
    mutate(lc_time_period = lc_period,
           # calculate climate change metrics
           temp_change = temp_after - temp_before,
           precip_change = precip_after - precip_before,
           precip_change_percent = ((precip_after - precip_before) / precip_before) * 100,
           # calculate mean values
           temp_mean = (temp_before + temp_after) / 2,
           precip_mean = (precip_before + precip_after) / 2)
  
  # store processed climate data in list
  climate_data_list[[lc_period]] <- period_climate_data
  
  # track which period was processed
  cat("  Processed", nrow(period_climate_data), "cells for period", lc_period, "\n")
  
}

## 8.3. Combine all climate data -----------------------------------------------

# Combine all climate periods into a signle df
all_climate_data <- bind_rows(climate_data_list)

# Check dimensions of the df
nrow(all_climate_data) # total cells = 
length(unique(all_climate_data$lc_time_period)) # periods =

# Check for any missing values
climate_summary <- all_climate_data |>
  group_by(lc_time_period) |>
  summarise(n_cells = n(),
            temp_before_na = sum(is.na(temp_before)),
            temp_after_na = sum(is.na(temp_after)),
            precip_before_na = sum(is.na(precip_before)),
            precip_after_na = sum(is.na(precip_after)),
            .groups = 'drop')

# Check the summary
summary(climate_summary)

## 8.4. Merge climate data with turnover data ----------------------------------

# Check structure of the turnover df
nrow(all_periods_turnover_all_land_cover_chanegs_15km)
unique(all_periods_turnover_all_land_cover_chanegs_15km$lc_time_period)

# Merge climate data with the existing turnover dataframe
all_periods_turnover_with_climate <- all_periods_turnover_all_land_cover_chanegs_15km |>
  left_join(all_climate_data, by = c("cell_ID", "lc_time_period"))

# Check it all went ok
nrow(all_periods_turnover_with_climate)
cat("Climate variables added:", sum(c("temp_before", "temp_after", "precip_before", 
                                      "precip_after", "temp_change", "precip_change") %in% 
                                      names(all_periods_turnover_with_climate)), "/ 6\n")

## 8.5. Quality checks ---------------------------------------------------------

# Check for missing climate data
climate_na_summary <- all_periods_turnover_with_climate |>
  summarise(temp_before_na = sum(is.na(temp_before)),
            temp_after_na = sum(is.na(temp_after)),
            precip_before_na = sum(is.na(precip_before)),
            precip_after_na = sum(is.na(precip_after)),
            total_rows = n())

print(climate_na_summary)

# Check value ranges
cat("\nClimate variable ranges:\n")
cat("Temperature before (°C):", round(range(all_periods_turnover_with_climate$temp_before, na.rm = TRUE), 2), "\n")
cat("Temperature after (°C):", round(range(all_periods_turnover_with_climate$temp_after, na.rm = TRUE), 2), "\n")
cat("Temperature change (°C):", round(range(all_periods_turnover_with_climate$temp_change, na.rm = TRUE), 2), "\n")
cat("Precipitation before (mm/yr):", round(range(all_periods_turnover_with_climate$precip_before, na.rm = TRUE), 0), "\n")
cat("Precipitation after (mm/yr):", round(range(all_periods_turnover_with_climate$precip_after, na.rm = TRUE), 0), "\n")
cat("Precipitation change (mm/yr):", round(range(all_periods_turnover_with_climate$precip_change, na.rm = TRUE), 0), "\n")

# Save dataframe of turnove + the new climate data
save(all_periods_turnover_with_climate,
     file = here("data", "derived_data", 
                 "all_periods_turnover_all_land_cover_climate_15km.rda"))

# One last check!
cat("\n✅ Climate integration complete!\n")
cat("Enhanced dataset saved with", ncol(all_periods_turnover_with_climate), "variables\n")
cat("New climate variables added:\n")
climate_vars <- c("temp_before", "temp_after", "temp_change", "temp_mean",
                  "precip_before", "precip_after", "precip_change", 
                  "precip_change_percent", "precip_mean")
for(var in climate_vars) {
  if(var %in% names(all_periods_turnover_with_climate)) {
    cat(" ✓", var, "\n")
  }
}

# END OF SCRIPT ----------------------------------------------------------------