##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 3.1_turnover_land_cover
# This script contains code which calculates the temporal turnover of GBIF
# records in CORINE land cover pixels
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

## 1.1. Download files if neccessary -------------------------------------------

# CORINE Status Layers
download_files("https://ntnu.box.com/shared/static/97g9x4839ij4lnlldji2wh8e0e2lm5bf.tif", 
              "data/derived_data/norway_corine_change_modified_stack.tif")

# SSB Grid
download_files("https://ntnu.box.com/shared/static/pjb2qr9aptugu7awlqmomx2o9d40sdhp.zip", 
              "data/raw_data/norway_corine_change_modified_stack.tif")

## 1.2. Read in data -----------------------------------------------------------

# CORINE Status Layers
norway_corine_status_modified_stack <- rast(here("data", "derived_data",
                                                 "norway_corine_status_modified_stack.tif"))

# SSB Grid
ssb_grids <- vect(here("data", "raw_data",
                       "SSB050KM", "ssb50km.shp"))

# Cleaned occurrence records
occurrences_norway <- fread(here("data", "derived_data", 
                                 "cleaned_occurrences_july24.txt"))