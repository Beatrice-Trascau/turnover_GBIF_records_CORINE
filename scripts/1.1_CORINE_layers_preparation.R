##----------------------------------------------------------------------------##
# PAPER 1: CORINE LAND COVER CHANGES AND GBIF BIODIVERSITY RECORDS
# 1.1_corine_layers_preparation
# This script contains code which loads and prepares the CORINE land cover 
# CHANGE and STATUS layers for further analysis
##----------------------------------------------------------------------------##

# 1. DEFINE URLS, NAMES AND NEW VALUES FOR RASTERS -----------------------------

# URLs and Filenames for the CORINE Change Layers
change_urls <- c(
  "https://ntnu.box.com/shared/static/vduuevecunldbrc7jb60jarjts99c4db.tif",
  "https://ntnu.box.com/shared/static/nmn2kguk9ipx0u4a2a8yfvcozdf6g9ij.tif",
  "https://ntnu.box.com/shared/static/pah7ig013inqeepg3gwvfan9w00anitp.tif",
  "https://ntnu.box.com/shared/static/g9grkxsvv20sz48rkbig8f9tb8gennfy.tif",
  "https://ntnu.box.com/shared/static/v51lua6b9fph0k7bmsbbc1g20tkdjqh9.tif",
  "https://ntnu.box.com/shared/static/x7ck0jnagfoxvjxvxf99l9lhknky5xlt.tif")

change_filenames <- c(
  "U2006_CHA0006_00_V2020_20u1.tif",
  "U2006_CHA0006_06_V2020_20u1.tif",
  "U2012_CHA0612_06_V2020_20u1.tif",
  "U2012_CHA0612_12_V2020_20u1.tif",
  "U2018_CHA1218_12_V2020_20u1.tif",
  "U2018_CHA1218_18_V2020_20u1.tif")

# URLs and Filenames for the CORINE Status Layers
status_urls <- c(
  "https://ntnu.box.com/shared/static/ffmbbb89aikwg64tg9ei30c8fnf7chl2.tif",
  "https://ntnu.box.com/shared/static/2x6g9jaov5rex3u0xt3hq9mmy91d63ew.tif",
  "https://ntnu.box.com/shared/static/ut1pcbnj7xgfwv3ptahu5c3krdy24l7d.tif",
  "https://ntnu.box.com/shared/static/iub514rfjnkopg3nu4nc18j4axq5jfon.tif")

status_filenames <- c(
  "U2006_CLC2000_V2020_20u1.tif",
  "U2012_CLC2006_V2020_20u1.tif",
  "U2018_CLC2012_V2020_20u1.tif",
  "U2018_CLC2018_V2020_20u1.tif")

# Changes for classes in the rasters
class_modifications <- list(
  list(from = 1:11, to = 1),
  list(from = c(12, 18, 20), to = 80),
  list(from = 21, to = 103),
  list(from = c(23, 24, 25), to = 250),
  list(from = c(26, 27), to = 380),
  list(from = 29, to = 590),
  list(from = 32, to = 711),
  list(from = c(30, 31, 33, 34, 35, 36, 39, 40, 41, 43, 44, 127, 128), to = NA))

# 2. READ IN CORINE LAYERS -----------------------------------------------------

## 2.1. Download layers (if needed) ----
download_files(change_urls, change_filenames)
download_files(status_urls, status_filenames)


## 2.2. Read in  layers --------------------------------------------------------
corine_change_stack <- read_rasters(change_filenames)
corine_status_stack <- read_rasters(status_filenames)

# 3. CUT AND MASK CHANGE LAYERS TO NORWAY --------------------------------------

## 3.1. Download country shapefile ----
norway <- geodata::gadm(country = "NOR", level = 0, 
                        path = tempdir(),
                        version = "latest")
# Check shapefile
plot(norway)

# Save shapefile
writeVector(norway, here::here("data", "raw_data", "raw_norway_shapefile",
                               "norway.shp"))

## 3.2. Re-project Norway shapefile to match projection of CORINE layers -------

# Check projections
crs(norway, proj = TRUE)
crs(corine_change_stack[[1]], proj = TRUE)

# Reproject Norway shapefile to the CORINE layers
norway_corine_projection <- project(norway, crs(corine_change_stack))

# Check projection
crs(norway_corine_projection, proj = TRUE) #projection correct now

# Save re-projected norway shapefile
writeVector(norway, here::here("data", "derived_data", 
                               "reprojected_norway_shapefile", 
                               "norway_corine_projection.shp"))

## 3.3. Crop and mask CORINE stack to Norway -----------------------------------
norway_corine_change_stack <- crop_mask_to_norway(corine_change_stack, 
                                                  norway_corine_projection)
norway_corine_status_stack <- crop_mask_to_norway(corine_status_stack, 
                                                  norway_corine_projection)

# Save the cropped layers 
terra::writeRaster(norway_corine_change_stack, 
                   here("data", "derived_data", "norway_corine_change_stack.tif"))
terra::writeRaster(norway_corine_status_stack, 
                   here("data", "derived_data", "norway_corine_status_stack.tif"), 
                   overwrite = TRUE)

# 4. MODIFY VALUES OF CHANGE LAYERS TO HELP IDENTIFY LAND COVER CHANGES --------
#The class codes/values are changed to unique numbers which will help identify 
# the land cover transitions between years

## 4.1. Modify classes in change layers ----------------------------------------

norway_corine_change_modified <- modify_class_values(norway_corine_change_stack, 
                                                     class_modifications)

## 4.2. Modify classes in status layers ----------------------------------------
norway_corine_change_modified <- modify_class_values(norway_corine_change_stack, 
                                                     class_modifications)

## 4.3. Save modified layers ---------------------------------------------------
terra::writeRaster(norway_corine_change_modified, 
                   here("data", "derived_data", 
                        "norway_corine_change_modified_stack.tif"), 
                   overwrite = TRUE)
terra::writeRaster(norway_corine_status_modified, 
                   here("data", "derived_data",
                        "norway_corine_status_modified_stack.tif"), 
                   overwrite = TRUE)

# END OF SCRIPT ----------------------------------------------------------------