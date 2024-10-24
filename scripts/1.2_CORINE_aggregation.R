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

## 2.1. Aggregate rasters ------------------------------------------------------

# Aggregate raster to 1km x 1km
corine_1km <- aggregate(norway_corine_status_modified_stack, 
                        fact = 10, fun = calculate_counts, cores = 2)

# Aggregate to 5km x 5km
corine_5km <- aggregate(norway_corine_status_modified_stack, 
                        fact = 50, fun = calculate_counts, cores = 2)

# Aggregate to 15km x 15km
corine_15km <- aggregate(norway_corine_status_modified_stack, 
                        fact = 150, fun = calculate_counts, cores = 2)

## 2.2. Change layer names -----------------------------------------------------

# Define the list of rasters
corine_list <- list(corine_1km, corine_5km, corine_15km)

# Define the new names for the layers
new_names <- c("urban2000", "complex_agri2000", "agri_sig_veg2000", "forest2000",
               "moors2000", "woodland2000", "sparse_veg2000", "urban2006", 
               "complex_agri2006", "agri_sig_veg2006", "forest2006",
               "moors2006", "woodland2006", "sparse_veg2006", "urban2012", 
               "complex_agri2012", "agri_sig_veg2012", "forest2012",
               "moors2012", "woodland2012", "sparse_veg2012", "urban2018", 
               "complex_agri2018", "agri_sig_veg2018", "forest2018",
               "moors2018", "woodland2018", "sparse_veg2018")

# Loop through the list and rename the layers
for (i in seq_along(corine_list)) {
  names(corine_list[[i]]) <- new_names
}

# Re-assign the modified rasters back to the original variables
corine_1km <- corine_list[[1]]
corine_5km <- corine_list[[2]]
corine_15km <- corine_list[[3]]

# Write new rasters to file
terra::writeRaster(corine_1km, here("data", "derived_data", "corine_1km.tif"))
terra::writeRaster(corine_5km, here("data", "derived_data", "corine_5km.tif"))
terra::writeRaster(corine_15km, here("data", "derived_data", "corine_15km.tif"))

# 3. CREATE DF OF % CHANGES FOR AGGREGATIONS -----------------------------------

## 3.1. 1km --------------------------------------------------------------------
# Add change values to original df 
corine_1km_df <- as.data.frame(corine_1km, xy = TRUE) |>
  mutate(urban_change2000.2006 = (urban2006 - urban2000)/100,
         complex_agri_change2000.2006 = (complex_agri2006 - complex_agri2000)/100,
         agri_sig_veg_change2000.2006 = (agri_sig_veg2006 - agri_sig_veg2000)/100,
         forest_change2000.2006 = (forest2006 - forest2000)/100,
         moors_change2000.2006 = (moors2006 - moors2000)/100,
         woodland_change2000.2006 = (woodland2006 - woodland2000)/100,
         sparse_veg_change2000.2006 = (sparse_veg2006 - sparse_veg2000)/100,
         urban_change2006.2012 = (urban2012 - urban2006)/100,
         complex_agri_change2006.2012 = (complex_agri2012 - complex_agri2006)/100,
         agri_sig_veg_change2006.2012 = (agri_sig_veg2012 - agri_sig_veg2006)/100,
         forest_change2006.2012 = (forest2012 - forest2006)/100,
         moors_change2006.2012 = (moors2012 - moors2006)/100,
         woodland_change2006.2012 = (woodland2012 - woodland2006)/100,
         sparse_veg_change2006.2012 = (sparse_veg2012 - sparse_veg2006)/100,
         urban_change2012.2018 = (urban2018 - urban2012)/100,
         complex_agri_change2012.2018 = (complex_agri2018 - complex_agri2012)/100,
         agri_sig_veg_change2012.2018 = (agri_sig_veg2018 - agri_sig_veg2012)/100,
         forest_change2012.2018 = (forest2018 - forest2012)/100,
         moors_change2012.2018 = (moors2018 - moors2012)/100,
         woodland_change2012.2018 = (woodland2018 - woodland2012)/100,
         sparse_veg_change2012.2018 = (sparse_veg2018 - sparse_veg2012)/100)




















