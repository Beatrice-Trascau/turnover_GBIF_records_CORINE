##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS 
# 3.1_Chelsa_download
# This script contains code which downloads Chelsa Precipitation and Daily
# Mean Near-Surface Air Temperature (tas)
##----------------------------------------------------------------------------##

# 1. DOWNLOAD WORLDCLIM DATA ---------------------------------------------------

# Define directory for the data to be downloaded in
worldclim_dir <- here("data", "raw_data", "chelsa")

# Define years and months
years <- 1997:2021
months <- sprintf("%02d", 1:12)  # 01-12 with leading zeros

# Define base URL for CHELSA v2.1
base_url <- "https://os.unil.cloud.switch.ch/chelsa/chelsa_V2/GLOBAL/monthly/"

# Function to download CHELSA file
download_chelsa_file <- function(variable, year, month, dest_dir) {
  # construct filename based on CHELSA naming convention
  if(variable == "tas") {
    filename <- paste0("CHELSA_", variable, "_", month, "_", year, "_V.2.1.tif")
    url <- paste0(base_url, variable, "/", filename)
  } else if(variable == "pr") {
    filename <- paste0("CHELSA_", variable, "_", month, "_", year, "_V.2.1.tif")
    url <- paste0(base_url, variable, "/", filename)
  }
  
  # define destination path
  dest_path <- file.path(dest_dir, variable, paste0(year))
  if(!dir.exists(dest_path)) dir.create(dest_path, recursive = TRUE)
  
  dest_file <- file.path(dest_path, filename)
  
  # download if file doesn't exist
  if(!file.exists(dest_file)) {
    cat("Downloading:", filename, "\n")
    tryCatch({
      download.file(url, dest_file, mode = "wb", quiet = TRUE)
      return(TRUE)
    }, error = function(e) {
      cat("  Error downloading", filename, ":", e$message, "\n")
      return(FALSE)
    })
  } else {
    cat("  File already exists:", filename, "\n")
    return(TRUE)
  }
}

# Download all files
cat("Starting CHELSA data download...\n")
cat("This will download", length(years) * 12 * 2, "files\n\n")

for(var in c("tas", "pr")) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("Downloading", var, "data\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  for(year in years) {
    cat("Year", year, ":\n")
    for(month in months) {
      success <- download_chelsa_file(var, year, month, chelsa_dir)
      if(!success) {
        cat("  Warning: Failed to download", var, "for", year, "-", month, "\n")
      }
    }
  }
}

# END OF SCRIPT ----------------------------------------------------------------