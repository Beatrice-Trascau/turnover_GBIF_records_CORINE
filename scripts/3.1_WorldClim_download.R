##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS 
# 3.1_WorldClim_download
# This script contains code which downloads WorldClim minimum temperature,
# maximum temperature and precipitation
##----------------------------------------------------------------------------##

# 1. DOWNLOAD WORLDCLIM DATA ---------------------------------------------------

# Historical WorldClim maximum temperature, minimum temperaute and precipitation were downloaded from https://worldclim.org/data/monthlywth.html on 08.11.2025 and uploaded to GoogleDrive
# Donwload folders for eah period and variable




# 1. LOAD NORWAY SHAPEFILE -----------------------------------------------------

# Read in the reprojected shapefile
norway_corine_projection <- vect(here("data", "derived_data", 
                                      "reprojected_norway_shapefile", 
                                      "norway_corine_projection.shp"))

# Convert to sf object
norway_sf <- st_as_sf(norway_corine_projection)

# Check CRS
st_crs(norway_sf)$input

# Check geometry type
st_geometry_type(norway_sf, by_geometry = FALSE)

# 3. DOWNLOAD CHELSA BIOCLIMATIC VARIABLES -------------------------------------

## 3.1. Set up download parameters ---------------------------------------------

# Define directory for the download
chelsa_raw_directory <- here("data", "raw_data", "chelsa")

# Define directory for the processed CHELSA fiels
chelsa_derived_directory <- here("data", "derived_data", "chelsa")

# Define bioclimatic variables to download
bioclimatic_variables <- paste0("bio", 1:19)