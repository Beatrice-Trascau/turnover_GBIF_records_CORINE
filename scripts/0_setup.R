##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS 
# 0_setup
# This script contains code which loads/installs necessary packages and defines
# functions used in the analysis
##----------------------------------------------------------------------------##

# 1. FUNCTION TO LOAD/INSTALL PACKAGES NEEDED FOR ANALYIS ----------------------

# Define function
install_load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

# Define list of packages
package_vec <- c("here", "terra", "sf", "geodata", "mapview",
                 "tidyverse", "dplyr", "ggplot2", "ggalluvial",
                 "networkD3", "gt", "cowplot", "data.table",
                 "tidyterra", "patchwork", "styler", "scales",
                 "plotly", "lme4", "DHARMa", "glmmTMB", "mgcv",
                 "tidyterra", "ggspatial", "htmlwidgets",
                 "htmltools", "patchwork", "webshot2",
                 "rgbif", "CoordinateCleaner", "codyn",
                 "gratia", "lattice", "car") # specify packages

# Execute the function
sapply(package_vec, install_load_package)

# 2. FUNCTION TO ONLY DOWNLOAD FILES THAT ARE NOT ALREADY IN THE FOLDERS -------

download_files <- function(urls, filenames, dir = here("data", "raw_data")) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  for (i in seq_along(urls)) {
    file_path <- file.path(dir, filenames[i])
    if (!file.exists(file_path)) {
      download.file(urls[i], file_path)
    }
  }
}

# 3. FUNCTION TO READ RASTERS --------------------------------------------------

read_rasters <- function(filenames, dir = here("data/raw_data")) {
  rasters <- lapply(filenames, function(x) {
    file_path <- file.path(dir, x)
    if (!file.exists(file_path)) {
      stop(paste("File does not exist:", file_path))
    }
    rast(file_path)
  })
  return(do.call(c, rasters))
}

# 4. FUNCTION TO CROP AND MASK RASTERS TO NORWAY -------------------------------

modify_class_values <- function(raster_stack, class_modifications) {
  modified_stack <- raster_stack
  for (mod in class_modifications) {
    modified_stack <- app(modified_stack, fun = function(x) {
      x[x %in% mod$from] <- mod$to
      return(x)
    })
  }
  return(modified_stack)
}

# 5. FUNCTION TO MODIFY CLASSES IN RASTERS -------------------------------------

modify_class_values <- function(raster_stack, class_modifications) {
  modified_stack <- raster_stack
  for (mod in class_modifications) {
    modified_stack <- app(modified_stack, fun = function(x) {
      x[x %in% mod$from] <- mod$to
      return(x)
    })
  }
  return(modified_stack)
}

# 6. FUNCTION TO COUNT THE NUMBER OF SMALL PIXELS WITH CERTAIN VALUES ----------

calculate_counts <- function(x) {
  # remove NA values before counting
  x <- na.omit(x)
  
  # count occurrences of each land cover category value
  counts <- table(factor(x, levels = c(1, 80, 103, 250, 380, 590, 711)))
  
  # return the counts as a numeric vector
  return(as.numeric(counts))
}

# 7. FUNCTION TO EXTRACT LAND COVER FOR GIVEN PERIOD ---------------------------

extract_land_cover <- function(cell_ids, period) {
  if (period == "1997-2000") {
    land_cover <- terra::extract(norway_corine_status_modified_stack[[1]], cell_ids)[, 1]
  } else if (period == "2003-2006") {
    land_cover <- terra::extract(norway_corine_status_modified_stack[[2]], cell_ids)[, 1]
  } else if (period == "2009-2012") {
    land_cover <- terra::extract(norway_corine_status_modified_stack[[3]], cell_ids)[, 1]
  } else {
    land_cover <- NA
  }
  return(land_cover)
}

# 8. FUNCTION TO CALCULATE JACCARD DISSIMILARITY INDEX FOR PERIODS -------------

calculate_jaccard_for_periods <- function(df, start_period, end_period) {
  start_data <- df |>
    filter(period == start_period)
  end_data <- df |>
    filter(period == end_period)
  
  combined_data <- start_data |>
    full_join(end_data, by = "cell", suffix = c("_start", "_end"))
  
  # Extract land cover for start period
  combined_data <- combined_data |>
    mutate(land_cover_start = extract_land_cover(cell, start_period))
  
  jaccard_results <- combined_data |>
    group_by(cell) |>
    summarize(
      species_start = list(unique(species_start)),
      species_end = list(unique(species_end)),
      jaccard_dissimilarity = 1 - length(intersect(species_start[[1]], 
                                                   species_end[[1]])) / length(union(species_start[[1]], 
                                                                                     species_end[[1]])),
      land_cover_start = first(land_cover_start)) |>
    mutate(start_period = start_period, end_period = end_period)
  
  return(jaccard_results)
}

# 8. FUNCTION TO FILTER CELLS WITH MORE THAN 3 SPECIES -------------------------

# Filter cells with more than 3 species in every period (1km, 5km, 15km)
filter_cells_with_species_each_period <- function(occurrences_df) {
  # Count species per cell and period
  species_count <- occurrences_df |>
    st_drop_geometry() |>
    group_by(cell, period) |>
    summarise(species_count = n_distinct(species), .groups = 'drop')
  
  # Identify cells that have more than 3 species in every period
  cells_with_more_than_3_species_each_period <- species_count |>
    group_by(cell) |>
    filter(all(species_count > 3)) |>
    pull(cell) |>
    unique()
  
  # Filter occurrences for cells with more than 3 species in every period
  filtered_occurrences <- occurrences_df |>
    filter(cell %in% cells_with_more_than_3_species_each_period) |>
    st_drop_geometry() |>
    select(period, species, cell) |>
    na.omit()
  
  # Aggregate to ensure unique species records per cell and period
  unique_occurrences <- filtered_occurrences |>
    distinct(cell, species, period, .keep_all = TRUE)
  
  return(unique_occurrences)
}

# 9. FUNCTION TO CALCULATE JACCARDS' INDEX FOR EACH AGGREGATION AND PEROIOD ----

calculate_turnover_for_resolutions <- function(unique_occurrences_df) {
  # defined function to calculate Jaccard's dissimilarity index for each specified
  # period combination for the given data frame of unique occurrences
  jaccard_index_list <- lapply(period_combinations, function(periods) {
    # apply function to each pair of periods in the list "period_combinations"
    calculate_jaccard_for_periods(unique_occurrences_df, periods[1], periods[2])
    # use custom function to calculate the index between the two periods
  })
  return(bind_rows(jaccard_index_list))
}

# 10. FUNCTION TO ADD CELL COORDS AND LAND COVER CATEGORY ----------------------

# Function to add cell coordinates and land cover change values for each resolution
add_cell_coordinates_and_land_cover <- function(temporal_turnover_df, land_cover_id_raster, land_cover_df) {
  # Extract xy coordinates
  xy_coords <- terra::xyFromCell(land_cover_id_raster, 1:ncell(land_cover_id_raster))
  
  # Create df of the xy coordinates
  centroids_df <- data.frame(cell = 1:ncell(land_cover_id_raster), 
                             x = xy_coords[,1], y = xy_coords[,2])
  
  # Merge coordinates with temporal turnover data
  temporal_turnover_coords <- left_join(temporal_turnover_df, centroids_df, by = "cell")
  
  # Merge land cover change data with temporal turnover data
  temporal_turnover_with_lc <- left_join(temporal_turnover_coords, land_cover_df, by = "cell")
  
  return(temporal_turnover_with_lc)
}

# 11. FUNCTION TO CALCULATE FOREST -> TWS TRANSITIONS IN RASTERS ---------------

# Function to analyze forest transitions between two time periods
analyse_forest_transition <- function(rast_t1, rast_t2) {
  # Create transition raster
  # 0 = non-forest in t1
  # 1 = forest remained forest
  # 2 = forest converted to TWS
  # 3 = other forest conversion
  
  # Create "dummy" raster with matching spatial characteristics
  transition <- rast_t1
  
  # Give it 0 values to show non-forested areas
  transition[] <- 0
  
  # Identify forest cells in initial layer
  forest_t1 <- rast_t1 == 250
  
  # For forest cells in initial layer, categorize changes:
  transition[forest_t1] <- case_when(
    # Forest remained forest
    rast_t2[forest_t1] == 250 ~ 1,
    # Forest converted to shrubland
    rast_t2[forest_t1] == 590 ~ 2,
    # Forest converted to something else
    TRUE ~ 3
  )
  
  return(transition)
}


# 12. FUNCTION TO CALCULATE TWS -> FORESTS TRANSITIONS IN RASTERS --------------

# Function to analyze TWS transitions between two time periods
analyse_tws_transition <- function(rast_t1, rast_t2) {
  # Create transition raster
  # 0 = non-TWS in t1
  # 1 = TWS remained TWS
  # 2 = TWS converted to forest
  # 3 = other TWS conversion
  
  # Create "dummy" raster with matching spatial characteristics
  transition <- rast_t1
  
  # Give it 0 values to show non-forested areas
  transition[] <- 0
  
  # Identify TWS cells in initial layer
  tws_t1 <- rast_t1 == 590
  
  # For TWS cells in initial layer, categorize changes:
  transition[tws_t1] <- case_when(
    # TWS remained TWS
    rast_t2[tws_t1] == 590 ~ 1,
    # TWS converted to Forest
    rast_t2[tws_t1] == 250 ~ 2,
    # TWS converted to something else
    TRUE ~ 3
  )
  
  return(transition)
}

# 13. FUNCTION TO CALCULATE ALLL -> URBAN TRANSITIONS IN RASTERS ---------------

analyse_urban_conversion <- function(rast_t1, rast_t2) {
  # Create transition raster
  # 0 = already urban in t1
  # 1 = forest to urban
  # 2 = shrubland to urban
  # 3 = complex agriculture to urban
  # 4 = agriculture & vegetation to urban
  # 5 = moors, heathland & grassland to urban
  # 6 = sparse vegetation to urban
  # 7 = no conversion to urban
  
  # Create "dummy" raster with matching spatial characteristics
  transition <- rast_t1
  
  # Give it values = 7 = no conversion to urban
  transition[] <- 7
  
  # Mark existing urban areas as 0
  urban_t1 <- rast_t1 == 1
  transition[urban_t1] <- 0
  
  # Identify cells that became urban in t2
  urban_t2 <- rast_t2 == 1
  
  # For cells that became urban, categorize their original land cover
  transition[urban_t2 & !urban_t1] <- case_when(
    # Forest to urban
    rast_t1[urban_t2 & !urban_t1] == 250 ~ 1,
    # Shrubland to urban
    rast_t1[urban_t2 & !urban_t1] == 590 ~ 2,
    # Complex agriculture to urban
    rast_t1[urban_t2 & !urban_t1] == 80 ~ 3,
    # Agriculture & vegetation to urban
    rast_t1[urban_t2 & !urban_t1] == 103 ~ 4,
    # Moors, heathland & grassland to urban
    rast_t1[urban_t2 & !urban_t1] == 380 ~ 5,
    # Sparse vegetation to urban
    rast_t1[urban_t2 & !urban_t1] == 711 ~ 6,
    # This TRUE case should never be reached as we've covered all classes
    TRUE ~ 7
  )
  
  return(transition)
}

# END OF SCRIPT ----------------------------------------------------------------