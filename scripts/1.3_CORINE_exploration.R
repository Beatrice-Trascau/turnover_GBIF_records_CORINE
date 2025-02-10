##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 1.3_CORINE_exploration
# This script contains code which explores land cover transitions for each
# grain size considered
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# 100m resolution
forest_tws_100m <- rast(here("data", "derived_data", 
                             "clc_status_100m_forest_tws.tif"))
tws_forest_100m <- rast(here("data", "derived_data", 
                             "clc_status_100m_tws_forest.tif"))
all_urban_100m <- rast(here("data", "derived_data", 
                            "clc_status_100m_all_urban.tif"))

# 1km resolution
forest_tws_1km <- rast(here("data", "derived_data", 
                            "clc_status_1km_forest_tws.tif"))
tws_forest_1km <- rast(here("data", "derived_data", 
                            "clc_status_1km_tws_forest.tif"))
all_urban_1km <- rast(here("data", "derived_data", 
                           "clc_status_1km_all_urban.tif"))

# 5km resolution
forest_tws_5km <- rast(here("data", "derived_data", 
                            "clc_status_5km_forest_tws.tif"))
tws_forest_5km <- rast(here("data", "derived_data", 
                            "clc_status_5km_tws_forest.tif"))
all_urban_5km <- rast(here("data", "derived_data", 
                           "clc_status_5km_all_urban.tif"))

# 15km resolution
forest_tws_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_forest_tws.tif"))
tws_forest_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_tws_forest.tif"))
all_urban_15km <- rast(here("data", "derived_data", 
                            "clc_status_15km_all_urban.tif"))

# Load Norway shapefile
norway <- vect(here("data", "derived_data", "reprojected_norway_shapefile", 
                    "norway_corine_projection.shp"))