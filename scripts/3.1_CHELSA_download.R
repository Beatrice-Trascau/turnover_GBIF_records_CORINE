##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS 
# 3.1_CHELSA_download
# This script contains code which donwlodads CHELSA bioclimatic variables for
# use in downstream analyses
##----------------------------------------------------------------------------##

# 1. LOAD NORWAY SHAPEFILE -----------------------------------------------------

# Read in the reprojected shapefile
norway_corine_projection <- vect(here("data", "derived_data", 
                                      "reprojected_norway_shapefile", 
                                      "norway_corine_projection.shp"))

# Check CRS
crs(norway_corine_projection, proj = TRUE)

# 3. DOWNLOAD CHELSA BIOCLIMATIC VARIABLES -------------------------------------

## 3.1. Set up download parameters ---------------------------------------------

# Define directory for the download
chelsa_raw_directory <- here("data", "raw_data", "chelsa")

# Define directory for the processed CHELSA fiels
chelsa_derived_directory <- here("data", "derived_data", "chelsa")

# Define bioclimatic variables to download
bioclimatic_variables <- paste0("bio", 1:19)