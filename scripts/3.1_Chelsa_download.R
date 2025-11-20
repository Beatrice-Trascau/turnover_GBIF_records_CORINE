##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS 
# 3.1_Chelsa_download
# This script contains code which downloads Chelsa Precipitation and Daily
# Mean Near-Surface Air Temperature (tas)
##----------------------------------------------------------------------------##

# 1. DOWNLOAD WORLDCLIM DATA ---------------------------------------------------

# Define directory for the data to be downloaded in
chelsa_dir <- here("data", "raw_data", "chelsa")

# Define years and months
years <- 1997:2021
months <- sprintf("%02d", 1:12)  # 01-12 with leading zeros

# Define bucket URL for CHELSA v2.1
bucket_url <- "https://os.unil.cloud.switch.ch/chelsa02/"

# Define base URL (combine the bucket URL with chelsa/global/monthly; naming convention taken from webiste)
base_url <- paste0(bucket_url, "chelsa/global/monthly/")

# Define function to download CHELSA files
download_chelsa_file <- function(variable, year, month, dest_dir) {
  # construct filename and URL based on the CHELSA structure
  filename <- paste0("CHELSA_", variable, "_", month, "_", year, "_V.2.1.tif")
  url <- paste0(base_url, variable, "/", year, "/", filename)
  
  # define the destination path
  dest_path <- file.path(dest_dir, variable, paste0(year))
  
  # check if directory already exists and create it if it does not
  if(!dir.exists(dest_path)) dir.create(dest_path, recursive = TRUE)
  dest_file <- file.path(dest_path, filename)
  
  # download if file doesn't exist
  if(!file.exists(dest_file)){
    # add progress indication
    cat("Downloading:", filename, "\n")
    
    # set longer time-out for downloads
    old_timeout <- getOption("timeout")
    options(timeout = 900)
    
    # download file with error catching method
    tryCatch({
      # use mode wb for binary files - to make sure the data will not get corrupted
      download.file(url, dest_file, mode = "wb", quiet = FALSE)
      
      # check if the file was downloaded successfuly and that it is not empty
      if(file.exists(dest_file) && file.size(dest_file) > 1000) {  # At least 1KB
        file_size_mb <- round(file.size(dest_file) / (1024^2), 2)
        cat("  ✓ Successfully downloaded:", filename, "(", file_size_mb, "MB)\n")
        result <- TRUE
      } else {
        cat("  ✗ Download failed (empty or too small file):", filename, "\n")
        if(file.exists(dest_file)) file.remove(dest_file)
        result <- FALSE
      }
      
      # Reset timeout
      options(timeout = old_timeout)
      
      return(result)
      
    }, error = function(e){
      # reset timeout
      options(timeout = old_timeout)
      
      # add error message
      cat("  ✗ Error downloading", filename, ":", e$message, "\n")
      
      # remove partial download if it exists
      if(file.exists(dest_file)){
        # add tracking text
        cat("  Removing partial download...\n")
        # remove file
        file.remove(dest_file)
      }
      return(FALSE)
    })
  } else {
    # get file size in mb
    file_size_mb <- round(file.size(dest_file) / (1024^2), 2)
    # get text alert of progress
    cat("  File already exists:", filename, "(", file_size_mb, "MB)\n")
    return(TRUE)
  }
}

# Initialise counters
total_files <- length(years) * 12 * 2
downloaded_files <- 0
failed_files <- 0
start_time <- Sys.time()

# Run download function
for(var in c("tas", "pr")){
  # get text indicating which variable is downloaded
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("Downloading", switch(var, "tas" = "Temperature", "pr" = "Precipitation"), "data (", var, ")\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # download each year for the variable
  for(year in years){
    # get text indicating which year is downloaded
    cat("Year", year, ":\n")
    year_start <- Sys.time()
    
    # download each month of each year
    for(month in months){
      # download chelsa file and store in success vector
      success <- download_chelsa_file(var, year, month, chelsa_dir)
      
      # check if the download happened successfully 
      if(success){
        downloaded_files <- downloaded_files + 1
      } else {
        failed_files <- failed_files + 1
        # get a warning if they failed to download
        cat("  Warning: Failed to download", var, "for", year, "-", month, "\n")
      }
      
      # add small delay to let Franklin breathe a little bit
      Sys.sleep(0.5)
    }
    
    # get stats about when download ended & how long it took
    year_end <- Sys.time()
    year_duration <- round(difftime(year_end, year_start, units = "mins"), 2)
    cat("  Year", year, "completed in", year_duration, "minutes\n\n")
  }
}

# Get total duration of download
end_time <- Sys.time()
total_duration <- round(difftime(end_time, start_time, units = "hours"), 2)

# Get a summary of successes/failures
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("DOWNLOAD SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Total files attempted:", total_files, "\n")
cat("Successfully downloaded:", downloaded_files, "\n")
cat("Failed downloads:", failed_files, "\n")
cat("Success rate:", round((downloaded_files/total_files)*100, 1), "%\n")
cat("Total download time:", total_duration, "hours\n")

if(failed_files > 0) {
  cat("\n⚠ Note: You can re-run this script to retry failed downloads.\n")
  cat("The script will skip files that were already downloaded successfully.\n")
} else {
  cat("\n✓ All downloads completed successfully!\n")
}

cat("\nFiles saved to:", chelsa_dir, "\n")

# END OF SCRIPT ----------------------------------------------------------------