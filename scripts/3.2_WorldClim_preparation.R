##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS 
# 3.2_WorldClim_preparation
# This script contains code which reprojects, cuts and masks the WorldClim
# layers to Norway
##----------------------------------------------------------------------------##

# 1. LOAD NORWAY SHAPEFILE -----------------------------------------------------

# Read in Norway shapefile thatmatches CORINE projection
norway_corine_projection <- vect(here("data", "derived_data",
                                      "reprojected_norway_shapefile",
                                      "norway_corine_projection.shp"))

# 2. PROCESS EACH WORLDCLIM FOLDER ---------------------------------------------

# Define path to WorldClim folder
worldclim_dir <- here("data", "raw_data", "worldclim")

# Get all WorldClim folders
worldclim_folders <- list.dirs(worldclim_dir, full.names = FALSE, 
                               recursive = FALSE)

# Create lists to store the rasters for each variable
tmin_annual_rasters <- list()
tmax_annual_rasters <- list()
prec_annual_rasters <- list()

# Process each worldclim folder
for (i in seq_along(worldclim_folders)) {
  
  # get the folder
  folder <- worldclim_folders[i]
  # get the path to the folder
  folder_path <- file.path(worldclim_dir, folder)
  
  # add progress bar
  cat("\nProcessing folder", i, "of", length(worldclim_folders), ":", folder, "\n")
  
  # extract variable type from folder name
  variable_type <- case_when(str_detect(folder, "tmin") ~ "tmin",
                             str_detect(folder, "tmax") ~ "tmax",
                             str_detect(folder, "prec") ~ "prec",
                             TRUE ~ "unknown")
  
  # extract period from folder name
  period <- str_extract(folder, "\\d{4}-\\d{4}")
  
  # get all TIFF files in the folder
  tif_files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)
  
  # extract years from the file names
  file_info <- tibble(filepath = tif_files,
                      filename = basename(tif_files)) |>
    mutate(year_month = str_extract(filename, "\\d{4}-\\d{2}(?=\\.tif$)"),
           year = as.numeric(str_sub(year_month, 1, 4)),
           month = as.numeric(str_sub(year_month, 6, 7))) |>
    arrange(year, month)
  
  # get the unique years in the folder
  unique_years <- unique(file_info$year)
  
  # Process each year within the folder
  for (year in unique_years) {
    
    # add a progress tracker
    cat("    Processing year", year, "...\n")
    
    # get the 12 monthly files for the specific year
    year_files <- file_info |>
      filter(year == !!year) |>
      pull(filepath)
    
    # add a warning if it doesn't find 12 files per year
    if (length(year_files) != 12) {
      cat("      Warning: Found", length(year_files), "months for year", year, 
          "(expected 12)\n")
    }
    
    # read monthly rasters for the specific year
    monthly_stack <- rast(year_files)
    
    # reproject to CORINE CRS
    monthly_stack_reprojected <- project(monthly_stack, crs(norway_corine_projection))
    
    # crop and mask to Norway
    monthly_stack_norway <- crop(monthly_stack_reprojected, norway_corine_projection,
                                 mask = TRUE)
    
    # calculate annual summary for the year
    if(variable_type == "prec"){
      # sum across all months to get annual precipitation
      annual_value <- sum(monthly_stack_norway)
      names(annual_value) <- paste0(variable_type, "_", year, "_annual_total")
    } else {
      # calculate mean minimum and maximum temperature
      annual_value <- mean(monthly_stack_norway)
      names(annual_value) <- paste0(variable_type, "_", year, "_annual_mean")
    }
    
    # store in the lists
    if(variable_type == "tmin"){
      tmin_annual_rasters[[as.character(year)]] <- annual_value
    } else if(variable_type == "tmax"){
      tmax_annual_rasters[[as.character(year)]] <- annual_value
    } else if(variable_type == "prec"){
      prec_annual_rasters[[as.character(year)]] <- annual_value
    }
    
  }
  
  # progress tracker
  cat("  Completed processing for", folder, "\n")
  
}

# 3. CREATE ONE RASTER STACK FOR EACH VARIABLE ---------------------------------

# Stack all tmin years (make sure that you have more than 0 years processed)
if(length(tmin_annual_rasters) > 0){
  tmin_stack <- rast(tmin_annual_rasters)
  # order layers by years
  tmin_stack <- tmin_stack[[order(names(tmin_stack))]]
}

# Stack all tmax years
if(length(tmax_annual_rasters) > 0){
  tmax_stack <- rast(tmax_annual_rasters)
  # order layers by years
  tmax_stack <- tmax_stack[[order(names(tmax_stack))]]
}

# Stack all prec years
if(length(prec_annual_rasters) > 0){
  prec_stack <- rast(prec_annual_rasters)
  # order layers by years
  prec_stack <- prec_stack[[order(names(prec_stack))]]
}

# 4. SAVE PROCESSED WORLDCLIM VARIABLES ----------------------------------------

# Save minimum teperature stack
writeRaster(tmin_stack, 
            here("data", "derived_data", "worldclim_tmin_annual_norway.tif"),
            overwrite = TRUE)

# Save maximum temperature stack
writeRaster(tmax_stack, 
            here("data", "derived_data", "worldclim_tmax_annual_norway.tif"),
            overwrite = TRUE)

# Save annual precipitation stack
writeRaster(prec_stack, 
            here("data", "derived_data", "worldclim_prec_annual_norway.tif"),
            overwrite = TRUE)

# END OF SCRIPT ----------------------------------------------------------------