##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS 
# 0_setup
# This script contains code which loads/installs necessary packages and defines
# functions used in the analysis
##----------------------------------------------------------------------------##

# 1. LOAD/INSTALL PACKAGES NEEDED FOR ANALYIS ----------------------------------

# Function to check to install/load packages

# Define function
install_load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

# Define list of packages
package_vec <- c("here", "terra", "sf", "geodata", "mapview",
                 "tidyverse", "dplyr", "ggplot2", "gt", "cowplot", 
                 "data.table","tidyterra", "patchwork", "styler", 
                 "scales","plotly", "lme4", "DHARMa", "glmmTMB", 
                 "mgcv", "ggspatial", "htmlwidgets","htmltools",  
                 "webshot2", "rgbif", "CoordinateCleaner", "codyn",
                 "gratia", "lattice", "car", "kableExtra",
                 "betareg", "spdep", "corrplot", "leaflet",
                 "viridis", "DT", "broom", "nlme", "ordbetareg")

# Execute the function
sapply(package_vec, install_load_package)

# 2. DOWNLOAD FILES ------------------------------------------------------------

# Function to check if files are in directory and then download them if they
# aren't

download_files <- function(urls, filenames, dir = here("data", "raw_data")) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  for (i in seq_along(urls)) {
    file_path <- file.path(dir, filenames[i])
    if (!file.exists(file_path)) {
      download.file(urls[i], file_path)
    }
  }
}

# 3. READ RASTERS --------------------------------------------------------------

# Function to read rasters (from the raw_data file)

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

# 4. CROP AND MASK RASTERS TO NORWAY -------------------------------------------

# Function to crop and mask rasters to Norway

crop_mask_to_norway <- function(raster_stack, norway_shape) {
  return(crop(raster_stack, norway_shape, mask = TRUE))
}

# 5. MODIFY CLASSES IN RASTERS -------------------------------------------------

# Function to modify class values in the rasters

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

# 6. CALCULATE FOREST -> TWS TRANSITIONS ---------------------------------------

# Function to calculate Forest -> TWS transitions between two time periods
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


# 7. CALCULATE TWS -> FORESTS TRANSITIONS --------------------------------------

# Function to calculate TWS -> Forest transitions between two time periods
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

# 8. CALCULATE ALLL -> URBAN TRANSITIONS ---------------------------------------

# Function to calculate All -> Urban transitions between two time periods
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

# 9. AGGREGATE RASTERS TO LARGER GRIDS -----------------------------------------

# This function aggregates the transition rasters (created with the functions 
#   above) to larger grid sizes

aggregate_transitions <- function(transition_raster, factor) {
  # Create empty list to store aggregated count rasters
  agg_counts <- list()
  
  # Get the levels to know what we're counting
  cats <- levels(transition_raster)[[1]]
  
  # For each layer in the transition raster
  for(layer in 1:nlyr(transition_raster)) {
    # Create empty list for this time period
    layer_counts <- list()
    
    # For each transition category (except "no change" and "no conversion")
    for(val in cats$value) {
      # Create binary raster for this transition
      binary <- transition_raster[[layer]] == val
      
      # Aggregate and sum the occurrences
      count_rast <- terra::aggregate(binary, fact=factor, fun="sum")
      
      # Name the raster based on the category
      names(count_rast) <- paste0(names(transition_raster)[layer], "_",
                                  cats$class[cats$value == val])
      
      # Add to list
      layer_counts[[length(layer_counts) + 1]] <- count_rast
    }
    
    # Combine all counts for this time period
    agg_counts[[layer]] <- rast(layer_counts)
  }
  
  # Combine all time periods
  return(rast(agg_counts))
}

# 10. CHECK IF CELLS ARE NA ----------------------------------------------------

# Function to check if certain cells in a raster were already NA before the masking
#   this function is used to check that masking of cells with >50% of their area
#   outside of the boundary was done correctly
check_na_values <- function(raster_layer, coords) {
  # Extract values directly
  extracted <- terra::extract(raster_layer, as.matrix(coords[, c("x", "y")]))
  
  # Print the structure to understand the output
  print(str(extracted))
  
  # Check if extraction worked and what format it returned
  if (is.null(extracted)) {
    return(NA)
  } else if (is.data.frame(extracted)) {
    # If it's a data frame, find the value column (usually the second column)
    if (ncol(extracted) >= 2) {
      return(sum(is.na(extracted[,2])))
    } else if (ncol(extracted) == 1) {
      return(sum(is.na(extracted[,1])))
    }
  } else if (is.vector(extracted)) {
    # If it's a vector, count NAs directly
    return(sum(is.na(extracted)))
  }
  
  # Fallback
  return(NA)
}

# 14. CALCULATE JACCARD'S DISSIMILARITY INDEX ----------------------------------

# Function to calculate Jaccard dissimilarity for a period pair
calculate_jaccard_for_periods <- function(df, before_period, after_period, change_period) {
  # Get data for before and after periods
  before_data <- df |>
    filter(period == before_period)
  after_data <- df |>
    filter(period == after_period)
  
  # Combine data
  combined_data <- before_data |>
    full_join(after_data, by = "cell", suffix = c("_before", "_after"))
  
  # Calculate Jaccard dissimilarity and species counts
  jaccard_results <- combined_data |>
    group_by(cell) |>
    summarize(
      species_before = list(unique(species_before)),
      species_after = list(unique(species_after)),
      n_species_before = length(unique(species_before)),
      n_species_after = length(unique(species_after)),
      jaccard_dissimilarity = 1 - length(intersect(species_before[[1]], 
                                                   species_after[[1]])) / 
        length(union(species_before[[1]], 
                     species_after[[1]])),
      .groups = 'drop'
    ) |>
    mutate(
      before_period = before_period,
      after_period = after_period,
      change_period = change_period
    )
  
  return(jaccard_results)
}


# END OF SCRIPT ----------------------------------------------------------------