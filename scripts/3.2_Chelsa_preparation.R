##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS 
# 3.2_WorldClim_preparation
# This script contains code which reprojects, cuts and masks the WorldClim
# layers to Norway
##----------------------------------------------------------------------------##

# 1. LOAD DATA AND DIRECTORIES -------------------------------------------------

# Read in Norway shapefile thatmatches CORINE projection
norway_corine_projection <- vect(here("data", "derived_data",
                                      "reprojected_norway_shapefile",
                                      "norway_corine_projection.shp"))

# Define paths
chelsa_dir <- here("data", "raw_data", "chelsa")
chelsa_processed_dir <- here("data", "derived_data", "chelsa")

# Load CLC for later resampling of CHELSA
corine_reference <- rast(here("data", "derived_data", "clc_status_15km_forest_tws_masked.tif"))[[1]]

# Define months for processing
months <- sprintf("%02d", 1:12)  # 01-12 with leading zeros

# 2. PROCESS CHELSA VARIABLES --------------------------------------------------

# Define time peirods for analysis
time_periods <- list("1997-2000" = 1997:2000,
                     "2003-2006" = 2003:2006,
                     "2006-2009" = 2006:2009,
                     "2009-2012" = 2009:2012,
                     "2012-2015" = 2012:2015,
                     "2018-2021" = 2018:2021)

## 2.1. Process temperature data -----------------------------------------------

# Get tracking of processing progress
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PROCESSING TEMPERATURE DATA\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Create list to store period means
tas_period_means <- list()

# Process each time period
for(period_name in names(time_periods)){
  
  # tracking
  cat("\nProcessing temperature for period:", period_name, "\n")
  
  # create variables needed
  period_years <- time_periods[[period_name]]
  monthly_rasters <- list()
  
  # load all monthly rasters for this periods (all months from all years)
  for(year in period_years){
    # track which year is being processed
    cat("  Processing year", year, "...")
    months_found <- 0
    
    # process each month
    for(month in months){
      # give filename and path
      filename <- paste0("CHELSA_tas_", month, "_", year, "_V.2.1.tif")
      filepath <- file.path(chelsa_dir, "tas", as.character(year), filename)
      
      # read raster if file exists
      if(file.exists(filepath)){
        #read raster
        r <- rast(filepath)
        
        # reproject to CORINE CRS
        r_proj <- project(r, crs(norway_corine_projection), method = "bilinear")
        
        # crop and mask to Norway
        r_norway <- crop(r_proj, norway_corine_projection, mask = TRUE)
        
        # convert from Kelvin to Celsius
        r_celsius <- (r_norway / 10) - 273.15
        
        # add to all monthly raster list
        monthly_rasters[[length(monthly_rasters) + 1]] <- r_celsius
        months_found <- months_found + 1
      } else {
        # get warning message
        cat("\n    Warning: Missing file", filename)
      }
    }
    
    # get progress update
    cat(" ✓ (", months_found, "months)\n")
  }
  
  # calculate mean across all months in the period (without averaging per year first)
  if(length(monthly_rasters) > 0){
    # create a stack
    period_stack <- rast(monthly_rasters)
    
    # calculate the mean
    period_mean <- mean(period_stack, na.rm = TRUE)
    
    # change the name of the period mean
    names(period_mean) <- paste0("tas_mean_", gsub("-", "_", period_name))
    
    # resample to match CORINE resolution and extent
    period_mean_resampled <- resample(period_mean, corine_reference, method = "average")
    
    # add to the one raster stack that will be saved
    tas_period_means[[period_name]] <- period_mean_resampled
    
    # get a progress update
    cat("  Completed period", period_name, "- processed", length(monthly_rasters), "months across", length(period_years), "years\n")
  } else {
    cat("  No data found for period", period_name, "\n")
  }
}

# Stack all temperature period means
if(length(tas_period_means) > 0) {
  tas_stack <- rast(tas_period_means)
  cat("\n Temperature stack created with", nlyr(tas_stack), "layers\n")
} else {
  stop("No temperature data processed successfully!")
}

## 2.2. Process precipitaion data ----------------------------------------------

# Start a tracker
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PROCESSING PRECIPITATION DATA\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Create list to store period means
pr_period_means <- list()

# Process each time period
for(period_name in names(time_periods)){
  
  # tracking
  cat("\nProcessing precipitation for period:", period_name, "\n")
  
  # create variables needed
  period_years <- time_periods[[period_name]]
  annual_totals <- list()
  
  # calculate total annual precipiation, then mean across years in the time period
  for(year in period_years){
    
    # tracking
    cat("  Processing year", year, "...")
    monthly_rasters <- list()
    
    # process each month in the specific year
    for(month in months){
      
      # give filename and path
      filename <- paste0("CHELSA_pr_", month, "_", year, "_V.2.1.tif")
      filepath <- file.path(chelsa_dir, "pr", as.character(year), filename)
      
      # read raster
      if(file.exists(filepath)){
        # read raster
        r <- rast(filepath)
        
        # reproject to CORINE CRS
        r_proj <- project(r, crs(norway_corine_projection), method = "bilinear")
        
        # crop and mask to Norway
        r_norway <- crop(r_proj, norway_corine_projection, mask = TRUE)
        
        # convert from kg/m²/month to mm/month
        monthly_rasters[[length(monthly_rasters) + 1]] <- r_norway
      } else {
        # get a warning
        cat("\n    Warning: Missing file", filename)
      }
    }
    
    # calculate annual total precipitation
    if(length(monthly_rasters) == 12) {
      annual_total <- sum(rast(monthly_rasters), na.rm = TRUE)
      annual_totals[[length(annual_totals) + 1]] <- annual_total
      cat(" (", length(monthly_rasters), "months)\n")
    } else {
      cat(" (only", length(monthly_rasters), "months found)\n")
    }
  }
  
  # calculate mean annual precipitaiton for the period
  if(length(annual_totals) > 0) {
    # calculate period mean
    period_mean_annual <- mean(rast(annual_totals), na.rm = TRUE)
    names(period_mean_annual) <- paste0("pr_mean_annual_", gsub("-", "_", period_name))
    
    # resample to match CORINE resolution and extent
    period_mean_resampled <- resample(period_mean_annual, corine_reference, method = "average")
    
    # add to the one raster stack that will be saved
    pr_period_means[[period_name]] <- period_mean_resampled
    
    # get a progress update
    cat("  Completed period", period_name, "- processed", length(annual_totals), "years\n")
  } else {
    cat("  No complete years found for period", period_name, "\n")
  }
}

# Stack all precipitation period means
if(length(pr_period_means) > 0) {
  pr_stack <- rast(pr_period_means)
  cat("\nPrecipitation stack created with", nlyr(pr_stack), "layers\n")
} else {
  stop("No precipitation data processed successfully!")
}

## 2.3. Save processed rasters -------------------------------------------------

# Save temperature stack
tas_output <- file.path(chelsa_processed_dir, "chelsa_tas_period_means_norway.tif")
writeRaster(tas_stack, tas_output, overwrite = TRUE)

# Save precipitation stack
pr_output <- file.path(chelsa_processed_dir, "chelsa_pr_period_means_norway.tif")
writeRaster(pr_stack, pr_output, overwrite = TRUE)

# 3. VALIDATION ----------------------------------------------------------------

## 3.1. Summary validation -----------------------------------------------------

# Check that all layers are present
cat("Temperature layers:", nlyr(tas_stack), "/ Expected:", length(time_periods), "\n")
cat("Precipitation layers:", nlyr(pr_stack), "/ Expected:", length(time_periods), "\n")

# Check spatial alignment with CORINE
cat("CRS match (temperature):", compareGeom(tas_stack, corine_reference, crs = TRUE, ext = FALSE, rowcol = FALSE, res = FALSE), "\n")
cat("CRS match (precipitation):", compareGeom(pr_stack, corine_reference, crs = TRUE, ext = FALSE, rowcol = FALSE, res = FALSE), "\n")
cat("Resolution match (temperature):", all(abs(res(tas_stack) - res(corine_reference)) < 0.1), "\n")
cat("Resolution match (precipitation):", all(abs(res(pr_stack) - res(corine_reference)) < 0.1), "\n")
cat("Extent match (temperature):", compareGeom(tas_stack, corine_reference, crs = FALSE, ext = TRUE, rowcol = FALSE, res = FALSE), "\n")
cat("Extent match (precipitation):", compareGeom(pr_stack, corine_reference, crs = FALSE, ext = TRUE, rowcol = FALSE, res = FALSE), "\n")

# Check temperature value ranges
temp_ranges <- minmax(tas_stack)
for(i in 1:ncol(temp_ranges)) {
  cat(" ", names(tas_stack)[i], ":", round(temp_ranges[1,i], 2), "to", round(temp_ranges[2,i], 2), "°C\n")
}

# Check precipitation ranges
prec_ranges <- minmax(pr_stack)
for(i in 1:ncol(prec_ranges)) {
  cat(" ", names(pr_stack)[i], ":", round(prec_ranges[1,i], 0), "to", round(prec_ranges[2,i], 0), "mm/year\n")
}

# Check that the values make sense/ are reasonable
temp_issues <- global(tas_stack < -50 | tas_stack > 40, sum, na.rm = TRUE)
prec_issues <- global(pr_stack < 0 | pr_stack > 10000, sum, na.rm = TRUE)
cat("Temperature values outside -50°C to 40°C:", sum(temp_issues[,1]), "cells\n")
cat("Precipitation values outside 0 to 10000 mm/year:", sum(prec_issues[,1]), "cells\n")

# Check for NA values
temp_nas <- global(is.na(tas_stack), sum, na.rm = TRUE)
prec_nas <- global(is.na(pr_stack), sum, na.rm = TRUE)
cat("NA cells in temperature data:", sum(temp_nas[,1]), "\n")
cat("NA cells in precipitation data:", sum(prec_nas[,1]), "\n")

## 3.2. Visualisation ----------------------------------------------------------

# Plot temperature across time periods
par(mfrow = c(2, 3))
for(i in 1:nlyr(tas_stack)) {
  plot(tas_stack[[i]], main = paste("Temperature", names(tas_stack)[i]))
  plot(norway_corine_projection, add = TRUE)
}

# Plot precipitation across time periods
par(mfrow = c(2, 3))
for(i in 1:nlyr(pr_stack)) {
  plot(pr_stack[[i]], main = paste("Precipitation", names(pr_stack)[i]))
  plot(norway_corine_projection, add = TRUE)
}

# END OF SCRIPT ----------------------------------------------------------------

