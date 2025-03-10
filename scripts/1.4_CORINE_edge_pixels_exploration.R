##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 1.4_CORINE_edge_pixels_exploration
# This script contains code which identifies the pixels in the 15km rasters
# which fall partially outside of the Norway boundary
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# One of the 15km rasters
forest_tws_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_forest_tws.tif"))

# Norway shapefile
norway <- vex(here("data", "derived_data", "reprojected_norway_shapefile",
                   "norway_corine_projection.shp"))