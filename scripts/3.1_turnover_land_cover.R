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

# Extract cell IDs from CLC
occurrences_with_cell <- occurrences_with_periods |>
  mutate(cell_id = terra::extract(reference_grid_15km, 
                                  st_coordinates(occurrences_with_periods))[, 1])

# Check how many occurrences were successfully assigned to grid cells
n_with_cell <- sum(!is.na(occurrences_with_cell$cell_id))
n_total <- nrow(occurrences_with_cell)
cat("Occurrences assigned to grid cells:", n_with_cell, 
    "(", round(n_with_cell/n_total*100, 2), "% of total)\n")

# Filter out occurrences without cell assignment
occurrences_with_cell <- occurrences_with_cell |>
  filter(!is.na(cell_id))

# Verify spatial distribution of assigned cells
if(n_with_cell > 0) {
  # Count occurrences per cell as a quality check
  cell_counts <- occurrences_with_cell |>
    st_drop_geometry() |>
    count(cell_id) |>
    rename(occurrence_count = n)
  
  # Output summary statistics of occurrence distribution
  cat("Cell assignment summary:\n")
  cat("Number of cells with occurrences:", nrow(cell_counts), "\n")
  cat("Min occurrences per cell:", min(cell_counts$occurrence_count), "\n")
  cat("Max occurrences per cell:", max(cell_counts$occurrence_count), "\n")
  cat("Median occurrences per cell:", median(cell_counts$occurrence_count), "\n")
}

# 5. CALCULATE SPECIES TURNOVER AND RECORDER EFFORT ----------------------------
## 5.1. Create long format dataset for turnover calculation --------------------
# Transform to long format for easier analysis of temporal patterns
periods_long <- occurrences_with_cell |>
  # Select only needed columns for efficient processing
  select(species, cell_id, period_2000_2006, period_2006_2012, period_2012_2018) |>
  # Transform to long format
  pivot_longer(
    cols = starts_with("period_"),
    names_to = "lc_change_period",
    values_to = "time_period"
  ) |>
  # Remove NA values
  filter(!is.na(time_period)) |>
  # Extract the years from the lc_change_period for clearer identification
  mutate(
    lc_change_period = gsub("period_", "", lc_change_period)
  )

## 5.2. Calculate Jaccard's dissimilarity index for each cell and period -------
# List to store results for each period
turnover_results <- list()

# Define the period pairs based on our temporal analysis framework
period_pairs <- list(
  list(lc_period = "2000-2006", before = "1997-2000", after = "2006-2009"),
  list(lc_period = "2006-2012", before = "2003-2006", after = "2012-2015"),
  list(lc_period = "2012-2018", before = "2008-2012", after = "2015-2018")
)

# Calculate turnover for each land cover change period
for (pair in period_pairs) {
  # Filter data for this period
  period_data <- periods_long |>
    filter(lc_change_period == pair$lc_period & 
             time_period %in% c(pair$before, pair$after))
  
  # Group by cell and calculate species lists and counts
  cell_turnover <- period_data |>
    st_drop_geometry() |>  # Drop geometry for faster processing
    group_by(cell_id, time_period) |>
    summarize(
      species_list = list(unique(species)),
      n_species = length(unique(species)),
      n_occurrences = n(),
      .groups = "drop"
    ) |>
    # Pivot to get before and after columns
    # Use names_transform to replace dashes with dots in column names
    # This prevents R from interpreting dashes as subtraction operators
    pivot_wider(
      id_cols = cell_id,
      names_from = time_period,
      names_transform = list(time_period = ~gsub("-", ".", .)),
      values_from = c(species_list, n_species, n_occurrences)
    )
  
  # Create column names based on the specific periods
  # Replace dashes with dots to avoid R syntax issues
  before_species_col <- paste0("species_list_", gsub("-", ".", pair$before))
  after_species_col <- paste0("species_list_", gsub("-", ".", pair$after))
  before_count_col <- paste0("n_occurrences_", gsub("-", ".", pair$before))
  after_count_col <- paste0("n_occurrences_", gsub("-", ".", pair$after))
  before_species_count_col <- paste0("n_species_", gsub("-", ".", pair$before))
  after_species_count_col <- paste0("n_species_", gsub("-", ".", pair$after))
  
  # Calculate Jaccard's dissimilarity and effort metrics
  turnover_period <- cell_turnover |>
    rowwise() |>
    mutate(
      # Only calculate Jaccard's dissimilarity when sufficient data exists
      jaccard_dissimilarity = if (
        exists(before_species_col) && 
        exists(after_species_col) && 
        length(!!sym(before_species_col)) > 0 &&
        length(!!sym(after_species_col)) > 0
      ) {
        # Calculate Jaccard's dissimilarity
        # Jaccard = 1 - (intersection size / union size)
        1 - length(intersect(
          !!sym(before_species_col),
          !!sym(after_species_col)
        )) / length(union(
          !!sym(before_species_col), 
          !!sym(after_species_col)
        ))
      } else {
        NA_real_  # Can't calculate Jaccard's with missing data
      },
      
      # Calculate standardized recorder effort change
      # This measures the relative change in sampling intensity
      delta_recorder_effort = if (
        exists(before_count_col) && 
        exists(after_count_col) &&
        !is.na(!!sym(before_count_col)) &&
        !is.na(!!sym(after_count_col)) &&
        (!!sym(before_count_col) + !!sym(after_count_col)) > 0
      ) {
        (!!sym(before_count_col) - !!sym(after_count_col)) / 
          (!!sym(before_count_col) + !!sym(after_count_col))
      } else {
        NA_real_
      },
      
      # Calculate total recorder effort (sampling intensity)
      recorder_effort = if (
        exists(before_count_col) && 
        exists(after_count_col) &&
        !is.na(!!sym(before_count_col)) &&
        !is.na(!!sym(after_count_col))
      ) {
        !!sym(before_count_col) + !!sym(after_count_col)
      } else {
        NA_real_
      },
      
      # Record the time period for identification
      time_period = pair$lc_period
    ) |>
    ungroup()
  
  # Store results for this period
  turnover_results[[pair$lc_period]] <- turnover_period |>
    select(cell_id, jaccard_dissimilarity, delta_recorder_effort, 
           recorder_effort, time_period,
           !!sym(before_species_count_col), !!sym(after_species_count_col))
}

# Combine all periods into a single dataframe for analysis
all_turnover <- bind_rows(turnover_results)

# 6. PREPARE DATA FOR GLM ------------------------------------------------------

# Join turnover results with land cover transition values
occurrences_for_model <- all_turnover |>
  left_join(lc_values, by = "cell_id") |>
  rowwise() |>
  mutate(land_cover_change = case_when(time_period == "2000-2006" ~ sum(c_across(matches("2000-2006.*to.*")), 
                                                                       na.rm = TRUE), 
                                      time_period == "2006-2012" ~ sum(c_across(matches("2006-2012.*to.*")), 
                                                                       na.rm = TRUE),
                                      time_period == "2012-2018" ~ sum(c_across(matches("2012-2018.*to.*")), 
                                                                       na.rm = TRUE),
                                      TRUE ~ NA_real_),
         land_cover_no_change = case_when(time_period == "2000-2006" ~ sum(c_across(matches("2000-2006.*no_change")), 
                                                                           na.rm = TRUE), 
                                          time_period == "2006-2012" ~ sum(c_across(matches("2006-2012.*no_change")), 
                                                                           na.rm = TRUE),
                                          time_period == "2012-2018" ~ sum(c_across(matches("2012-2018.*no_change")), 
                                                                           na.rm = TRUE),
                                          TRUE ~ NA_real_),
         forest_to_tws = case_when(
           time_period == "2000-2006" ~ sum(c_across(matches("2000-2006.*Forest_to_TWS")), 
                                            na.rm = TRUE),
           time_period == "2006-2012" ~ sum(c_across(matches("2006-2012.*Forest_to_TWS")), 
                                            na.rm = TRUE),
           time_period == "2012-2018" ~ sum(c_across(matches("2012-2018.*Forest_to_TWS")), 
                                            na.rm = TRUE),
           TRUE ~ NA_real_),
         tws_to_forest = case_when(
           time_period == "2000-2006" ~ sum(c_across(matches("2000-2006.*TWS_to_Forest")), 
                                            na.rm = TRUE),
           time_period == "2006-2012" ~ sum(c_across(matches("2006-2012.*TWS_to_Forest")), 
                                            na.rm = TRUE),
           time_period == "2012-2018" ~ sum(c_across(matches("2012-2018.*TWS_to_Forest")), 
                                            na.rm = TRUE),
           TRUE ~ NA_real_),
         to_urban = case_when(
           time_period == "2000-2006" ~ sum(c_across(matches("2000-2006.*to_urban")), 
                                            na.rm = TRUE),
           time_period == "2006-2012" ~ sum(c_across(matches("2006-2012.*to_urban")), 
                                            na.rm = TRUE),
           time_period == "2012-2018" ~ sum(c_across(matches("2012-2018.*to_urban")), 
                                            na.rm = TRUE),
           TRUE ~ NA_real_)) |>
  ungroup()

# Filter data to ensure quality
occurrences_turnover_15km <- occurrences_for_model |>
  # Remove rows with NA
  filter(!is.na(jaccard_dissimilarity),
         !is.na(land_cover_change),
         !is.na(land_cover_no_change),
         !is.na(delta_recorder_effort),
         !is.na(recorder_effort)) |>
  # Remove cells without a minimum recorder effort
  filter(n_species_1997.2000 >= 3 | n_species_2003.2006 >= 3 | n_species_2008.2012 >= 3,
         n_species_2006.2009 >= 3 | n_species_2015.2018 >= 3 | n_species_2012.2015 >= 3)

# Save to file
saveRDS(occurrences_turnover_15km,
        here("data", "derived_data", "occurrences_turnover_15km.rds"))