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
norway <- vect(here("data", "derived_data", "reprojected_norway_shapefile",
                   "norway_corine_projection.shp"))

# 2. GET PIXELS (PARTIALLY) OUTSIDE OF SHAPEFILE -------------------------------

# Create template raster with same extent as resolution from CLC
template_raster <- forest_tws_15km[[1]]

# Set all cell values to 1 (representing complete cells)
# This will help us in the next step to differentiate target pixels from those that are NA or 0
template_raster[] <- 1

# Convert cells to polygons
cell_polygons <- as.polygons(template_raster)

# Convert polygons to sf object
cell_polygons_sf <- st_as_sf(cell_polygons)

# Convert shapefile to sf object
norway_sf <- st_as_sf(norway)

# Calculate area of each cell
cell_polygons_sf$total_area <- st_area(cell_polygons_sf)

# Intersect cells with boundary to get proportion of cells within Norway
cells_within <- st_intersection(cell_polygons_sf, norway_sf)

# Calculate area of the cells within
cells_within$area_within <- st_area(cells_within)

# Calculate % of each cell inside the country
cells_within$percent_within <- as.numeric(cells_within$area_within/cells_within$total_area) * 100

# Calculate % of each cell outside the country
cells_within$percent_outside <- 100 - cells_within$percent_within

# Extract cell coordinates from the centroids
cell_centroids <- st_centroid(cells_within)
coords <- st_coordinates(cell_centroids)
cells_within$x <- coords[,1]
cells_within$y <- coords[,2]

# Exclude cells that are fully inside (percent_outside = 0)
edge_cells <- cells_within |>
  filter(percent_outside > 0) |>
  mutate(cell_id = row_number(),
         is_edge_cell = TRUE,
         is_mostly_outside = percent_outside > 50) |>
  mutate(area_km2 = as.numeric())