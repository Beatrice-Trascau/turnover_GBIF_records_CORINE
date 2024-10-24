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



# 3. 
corine_1km_df <- as.data.frame(corine_1km, xy = TRUE)
