##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 3.1_turnover_land_cover
# This script contains code which calculates the temporal turnover of GBIF
# records in CORINE land cover pixels
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

## 1.1. 100m CORINE Status Layers ----------------------------------------------

# Forest -> TWS
forest_tws_100m <- rast(here("data", "derived_data", 
                             "clc_status_100m_forest_tws.tif"))

# TWS -> Forest
tws_forest_100m <- rast(here("data", "derived_data", 
                             "clc_status_100m_tws_forest.tif"))

# All -> Urban
all_urban_100m <- rast(here("data", "derived_data", 
                            "clc_status_100m_all_urban.tif"))

## 1.2. Cleaned GBIF Occurrences -----------------------------------------------

# Cleaned occurrence records
occurrences_norway <- fread(here("data", "derived_data", 
                                 "clean_occurrences_100m.txt"))

# 2. SPATIAL REFERENCE GRID ----------------------------------------------------

# Create ID raster for cell extraction with same properties as CORINE
land_cover_id <- forest_tws_100m[[1]]

#Assign each cell a unique ID from 1 to ncell
land_cover_id[] <- 1:ncell(forest_tws_100m[[1]])

# Convert CORINE raster to spatial dataframe with correct cell IDs
land_cover_grid <- as.data.frame(forest_tws_100m, xy = TRUE) |>
  mutate(cell = terra::extract(land_cover_id, cbind(x, y))[,1]) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(forest_tws_100m))

# 3. PREPARE OCCURRENCES FOR ANALYSIS ------------------------------------------

## 3.1. Convert occurrences to spatial data ------------------------------------

# Convert occurrence records df to sf object
occurrences_sf <- st_as_sf(occurrences_norway,
                           coords = c("decimalLongitude", "decimalLatitude"),
                           crs = 4326)

# Re-project to match CORINE crs
occurrences_sf_reprojected <- st_transform(occurrences_sf, 
                                           st_crs(forest_tws_100m))

## 3.2. Add before and after periods -------------------------------------------

# Filtering and classification is done step by step to avoid overwriting of
  # years that are shared between periods

# First period (2000-2006 LC change, 1997-2000 = before, 2006-2012 = after)
period1_occurrences <- occurrences_sf_reprojected |>
  filter(year %in% c(1997:2000, 2006:2009)) |>
  mutate(period = case_when(
    year %in% 1997:2000 ~ "1997-2000",
    year %in% 2006:2009 ~ "2006-2009"
  ))

# Second period (2006-2012 LC change, 2003-2006 = before, 2012-2015 = after)
period2_occurrences <- occurrences_sf_reprojected |>
  filter(year %in% c(2003:2006, 2012:2015)) |>
  mutate(period = case_when(
    year %in% 2003:2006 ~ "2003-2006",
    year %in% 2012:2015 ~ "2012-2015"
  ))

# Third period (2012-2018 LC change, 2008-2012 = before, 2015-2018 = after)
period3_occurrences <- occurrences_sf_reprojected |>
  filter(year %in% c(2008:2012, 2015:2018)) |>
  mutate(period = case_when(
    year %in% 2008:2012 ~ "2008-2012",
    year %in% 2015:2018 ~ "2015-2018"
  ))

# Combine all periods
clean_occurrences_sf <- bind_rows(
  period1_occurrences,
  period2_occurrences,
  period3_occurrences
)
