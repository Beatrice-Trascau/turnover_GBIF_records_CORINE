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
            here("data", "derived_data", "worldclim",
                 "worldclim_tmin_annual_norway.tif"),
            overwrite = TRUE)

# Save maximum temperature stack
writeRaster(tmax_stack, 
            here("data", "derived_data", "worldclim",
                 "worldclim_tmax_annual_norway.tif"),
            overwrite = TRUE)

# Save annual precipitation stack
writeRaster(prec_stack, 
            here("data", "derived_data", "worldclim",
                 "worldclim_prec_annual_norway.tif"),
            overwrite = TRUE)

# 5. VALIDATE CREATED RASTERS --------------------------------------------------

## 5.1. Check basic raster properties are correct ------------------------------

# Check the number of layers in each stack
nlyr(tmin_stack) #35
nlyr(tmax_stack) #35
nlyr(prec_stack) #35 - All 3 are as expected!

## 5.2. Check layer names and years --------------------------------------------

# Extract years from layer names - with function
extract_years <- function(stack_names) {
  str_extract(stack_names, "\\d{4}") |> as.numeric()
}

tmin_years <- extract_years(names(tmin_stack))
tmax_years <- extract_years(names(tmax_stack))
prec_years <- extract_years(names(prec_stack))

# Check which years there are in the layers
paste(sort(tmin_years), collapse = ", ") 
paste(sort(tmax_years), collapse = ", ") # 1990-2024 for all 3 stacks - correct!
paste(sort(prec_years), collapse = ", ") # For further analyis: need to remove <1997 & >2021

## 5.3. Check spatial properties -----------------------------------------------

# Load CORINE to use as example
corine <- rast(here("data", "derived_data", "clc_status_15km_forest_tws_masked.tif"))

# Check extents
as.vector(ext(tmin_stack))
as.vector(ext(tmax_stack))
as.vector(ext(prec_stack)) # All stacks have the same extent
as.vector(ext(corine[[1]])) # CLC has slightly different values but nothing drastic

# Check resolutions
res(tmin_stack)
res(tmax_stack)
res(prec_stack) # All have the same values (3711.318, 3711.318)
res(corine[[1]]) # 15000, 15000

# Check CRS - all looks well!
cat("  tmin matches CORINE:", compareGeom(tmin_stack, corine[[1]], crs = TRUE, ext = FALSE, rowcol = FALSE, res = FALSE), "\n")
cat("  tmax matches CORINE:", compareGeom(tmax_stack, corine[[1]], crs = TRUE, ext = FALSE, rowcol = FALSE, res = FALSE), "\n")
cat("  prec matches CORINE:", compareGeom(prec_stack, corine[[1]], crs = TRUE, ext = FALSE, rowcol = FALSE, res = FALSE), "\n")

## 5.4. Check if values are logical --------------------------------------------

# Check range for tmin values
tmin_range <- minmax(tmin_stack)

# Check range for tmin values
tmax_range <- minmax(tmax_stack)

# Check range for prec values
prec_range <- minmax(prec_stack)

# Check for NA values
for (i in 1:nlyr(tmin_stack)) {
  na_count <- global(is.na(tmin_stack[[i]]), sum)
  if (na_count[1,1] > 0) {
    cat("  tmin layer", i, "(", names(tmin_stack)[i], "):", na_count[1,1], "NA cells\n")
  }
}
cat("  Total NA cells in tmin:", sum(global(is.na(tmin_stack), sum)[,1]), "\n")
cat("  Total NA cells in tmax:", sum(global(is.na(tmax_stack), sum)[,1]), "\n")
cat("  Total NA cells in prec:", sum(global(is.na(prec_stack), sum)[,1]), "\n")
# there seem to be a lot of NAs

## 5.5. Logical checks ---------------------------------------------------------

# Check if tmin < tmax for all years
issues_found <- FALSE
for (year in tmin_years) {
  # Find layers for this year
  tmin_layer <- grep(as.character(year), names(tmin_stack), value = FALSE)[1]
  tmax_layer <- grep(as.character(year), names(tmax_stack), value = FALSE)[1]
  
  if (!is.na(tmin_layer) & !is.na(tmax_layer)) {
    # Check if any cell has tmin >= tmax
    diff_raster <- tmax_stack[[tmax_layer]] - tmin_stack[[tmin_layer]]
    problem_cells <- global(diff_raster <= 0, sum, na.rm = TRUE)[1,1]
    
    if (problem_cells > 0) {
      cat(" Year", year, ":", problem_cells, "cells where tmin >= tmax\n")
      issues_found <- TRUE
    }
  }
} # Year 2004 : 3 cells where tmin >= tmax ???
if (!issues_found) {
  cat("  ✓ All cells have tmax > tmin for all years\n")
}

# Check for unreasonable values
unrealistic_tmin <- global(tmin_stack < -50 | tmin_stack > 40, sum, na.rm = TRUE)
unrealistic_tmax <- global(tmax_stack < -50 | tmax_stack > 40, sum, na.rm = TRUE) # both look ok

# Check that precipitation is above 0 and lower than 10 000 mm
unrealistic_prec <- global(prec_stack < 0 | prec_stack > 10000, sum, na.rm = TRUE) # 0

## 5.6. Visual check -----------------------------------------------------------

# Plot first year of each variable for visual inspection
par(mfrow = c(2, 2))

# Plot tmin
plot(tmin_stack[[1]], main = paste("Min Temp", names(tmin_stack)[1]))
plot(norway_corine_projection, add = TRUE)

# Plot tmax
plot(tmax_stack[[1]], main = paste("Max Temp", names(tmax_stack)[1]))
plot(norway_corine_projection, add = TRUE)

# Plot prec
plot(prec_stack[[1]], main = paste("Precipitation", names(prec_stack)[1]))
plot(norway_corine_projection, add = TRUE)

# Plot difference (tmax - tmin) for first year
temp_diff <- tmax_stack[[1]] - tmin_stack[[1]]
plot(temp_diff, main = "Temperature Range (tmax - tmin)")
plot(norway_corine_projection, add = TRUE)

# END OF SCRIPT ----------------------------------------------------------------