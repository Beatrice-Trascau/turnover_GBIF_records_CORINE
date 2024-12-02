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

# Filter cells with more than 3 species (1km, 5km, 15km)
filter_cells_with_species <- function(occurrences_df) {
  species_count <- occurrences_df |>
    st_drop_geometry() |>
    group_by(cell, period) |>
    summarise(species_count = n_distinct(species), .groups = 'drop')
  
  cells_with_more_than_3_species <- species_count |>
    filter(species_count > 3) |>
    pull(cell)
  
  filtered_occurrences <- occurrences_df |>
    filter(cell %in% cells_with_more_than_3_species) |>
    st_drop_geometry() |>
    select(period, species, cell) |>
    na.omit()
  
  # Aggregate to ensure unique species records per cell and period
  unique_occurrences <- filtered_occurrences |>
    distinct(cell, species, period, .keep_all = TRUE)
  
  return(unique_occurrences)
}

# END OF SCRIPT ----------------------------------------------------------------