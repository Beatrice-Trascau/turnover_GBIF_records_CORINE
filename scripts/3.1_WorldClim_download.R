##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS 
# 3.1_WorldClim_download
# This script contains code which downloads WorldClim minimum temperature,
# maximum temperature and precipitation
##----------------------------------------------------------------------------##

# 1. DOWNLOAD WORLDCLIM DATA ---------------------------------------------------

# Define directory for the data
worldclim_directory <- here("data", "raw_data", "worldclim")

# Historical WorldClim maximum temperature, minimum temperaute and precipitation were downloaded from https://worldclim.org/data/monthlywth.html on 08.11.2025 and uploaded to GoogleDrive
# Download folders for each period and variable

# Authenticate on Google Drive
drive_auth()

# Find Worlclim folder on Google Drive
worldclim_parent <- drive_find(pattern = "WorldClim", type = "folder") |>
  filter(name == "WorldClim") |>
  slice(1) # just in case it finds more than one folder - it shouldn't
cat("Using folder:", worldclim_parent$name, "\n")
cat("Folder ID:", worldclim_parent$id, "\n\n")

# Get all the folders within the WorldClim folder
worldclim_folders <- drive_ls(as_id(worldclim_parent$id), type = "folder")

# Check how many and which folders were found
cat("Found", nrow(worldclim_folders), "folders to download:\n")
print(worldclim_folders |> select(name))

# Download each folder - use function defined in the 0_setup.R script
for(i in seq_len(nrow(worldclim_folders))){
  folder_info <- worldclim_folders[i, ]
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("Processing folder", i, "of", nrow(worldclim_folders), ":", 
      folder_info$name, "\n")
  
  download_drive_folder(folder_info$id, 
                        folder_info$name, 
                        worldclim_dir)
}

# 2. CROP AND MASK CLIMATIC VARIABLES TO NORWAY --------------------------------

# Load Norway shapefile
norway_corine_projection <- vect(here("data", "derived_data", 
                                      "reprojected_norway_shapefile", 
                                      "norway_corine_projection.shp"))

# Check CRS
st_crs(norway_corine_projection)$input







