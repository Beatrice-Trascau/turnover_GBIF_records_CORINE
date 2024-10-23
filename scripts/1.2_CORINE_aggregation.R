##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 1.2_CORINE_aggregation
# This script contains code which aggregates CORINE land cover pixels to higher
# pixel sizes
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

## 1.1. Download files if neccessary -------------------------------------------

# CORINE Status Layers
download_files("https://ntnu.box.com/shared/static/97g9x4839ij4lnlldji2wh8e0e2lm5bf.tif", 
               "data/derived_data/norway_corine_change_modified_stack.tif")

## 1.2. Read in data -----------------------------------------------------------

# CORINE Status Layers
norway_corine_status_modified_stack <- rast(here("data", "derived_data",
                                                 "norway_corine_status_modified_stack.tif"))

# 2. AGGREGATE CORINE TO HIGHER CELL SIZES -------------------------------------

# Aggregate raster to 1km x 1km
corine_1km <- aggregate(norway_corine_status_modified_stack, 
                        fact = 10, fun = calculate_counts, cores = 2)

# Aggregate to 5km x 5km
corine_5km <- aggregate(norway_corine_status_modified_stack, 
                        fact = 50, fun = calculate_counts, cores = 2)

# Aggregate to 15km x 15km
corine_15km <- aggregate(norway_corine_status_modified_stack, 
                        fact = 150, fun = calculate_counts, cores = 2)

# 