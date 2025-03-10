##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 1.4_CORINE_edge_pixels_exploration
# This script contains code which identifies the pixels in the 15km rasters
# which fall partially outside of the Norway boundary
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# Forest -> TWS
forest_tws_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_forest_tws.tif"))

# TWS -> Forest
tws_forest_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_tws_forest.tif"))

# All -> Urban
all_urban_15km <- rast(here("data", "derived_data", 
                            "clc_status_15km_all_urban_combined.tif"))

# Norway shapefile
norway <- vect(here("data", "derived_data", "reprojected_norway_shapefile",
                   "norway_corine_projection.shp"))

# 2. GET PIXELS (PARTIALLY) OUTSIDE OF SHAPEFILE -------------------------------

## 2.1. Prepare data -----------------------------------------------------------

# Create template raster with same extent as resolution from CLC
template_raster <- forest_tws_15km[[1]]

# Set all cell values to 1 (representing complete cells)
# This will help us in the next step to differentiate target pixels from those that are NA or 0
template_raster[] <- 1:ncell(template_raster)

# Convert cells to polygons
cell_polygons <- as.polygons(template_raster)

# Convert polygons to sf object
cell_polygons_sf <- st_as_sf(cell_polygons)

# Convert shapefile to sf object
norway_sf <- st_as_sf(norway)

## 2.2. Calculate area outside and within the boundary -------------------------

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
  select(cell_id, x, y, percent_outside, is_mostly_outside, geometry)

# Save edfe cell data to file
save(edge_cells, file = here("data", "derived_data", "edge_cells_15km.rda"))

# 3. CREATE SUMMARY INFORMATION TABLE ------------------------------------------

# Get the total number of cells in the raster
total_cells <- ncell(template_raster)

# Get the cells that have data
cells_with_data <- global(template_raster, function(x) sum(!is.na(x)))[1,1]

# Calculate % of cells that are edge cells
percent_edge <- (nrow(edge_cells) / cells_with_data) * 100

# Create df with summary statistics
edge_summary_stats_15km <- data.frame(Statistic = c("Total number of edge cells",
                                                    "Mean percent outside boundary",
                                                    "Median percent outside boundary",
                                                    "Minimum percent outside boundary",
                                                    "Maximum percent outside boundary",
                                                    "Number of cells mostly outside (>50%)",
                                                    "% of data cells that are edge cells",
                                                    "Total cells in raster",
                                                    "Cells with data"),
                                      Value = c(nrow(edge_cells),
                                                round(mean(edge_cells$percent_outside), 2),
                                                round(median(edge_cells$percent_outside), 2),
                                                round(min(edge_cells$percent_outside), 2),
                                                round(max(edge_cells$percent_outside), 2),
                                                sum(edge_cells$is_mostly_outside),
                                                round(percent_edge, 2),
                                                total_cells,
                                                cells_with_data))

# Save summary statistics dfs to file
save(edge_summary_stats_15km, file = here("data", "derived_data", "edge_summary_stats_15km.rda"))

# 4. VISUALISE EDGE PIXELS -----------------------------------------------------

## 4.1. Map of edge pixels -----------------------------------------------------

# Create map of the edge pixels
edge_map <- ggplot() +
  geom_sf(data = norway_sf, fill = "lightgrey", color = "black") +
  geom_sf(data = edge_cells, aes(color = percent_outside, size = percent_outside)) +
  scale_color_viridis_c(option = "magma", direction = -1,
                        name = "% Outside\nBoundary") +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_y = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering()) +
  annotation_scale(location = "br", width_hint = 0.35)

# Save as supplementary information (.png)
ggsave(filename = here("figures", "SupplementaryFigure1_Edge_Pixels_15km.png"),
       plot = edge_map, width = 10, height = 6, dpi = 300)

# Save as as supplementary information (.svg)
ggsave(filename = here("figures", "SupplementaryFigure1_Edge_Pixels_15km.svg"),
       plot = edge_map, width = 10, height = 6, dpi = 300)

## 4.2. Histogram of values for % outside of the the boundary ------------------

# Plot histogram
edge_histogram <- ggplot(edge_cells, aes(x=percent_outside)) + 
  geom_histogram(binwidth = 10, boundary = 0, 
                 fill = "#C3B1E1", color = "#e9ecef", alpha = 0.9) +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  labs(x = "% of raster cell outside of boundary",
       y = "Count") + 
  theme_classic()

# Save as supplementary information (.png)
ggsave(filename = here("figures", "SupplementaryFigure2_Edge_Pixels_15km_Histogram.png"),
       plot = edge_histogram, width = 10, height = 6, dpi = 300)

# Save as as supplementary information (.svg)
ggsave(filename = here("figures", "SupplementaryFigure2_Edge_Pixels_15km_Histogram.svg"),
       plot = edge_histogram, width = 10, height = 6, dpi = 300)

# 5. REMOVE CELLS WITH >50% OUTSIDE OF BOUNDARY --------------------------------

## 5.1. Identify cells to mask -------------------------------------------------

# Get coordinates of cells with >50% area outside
cells_to_mask <- edge_cells |>
  filter(is_mostly_outside == TRUE) |>
  st_drop_geometry() |>
  select(x, y)

## 5.2. Mask cells from Forest -> TWS ------------------------------------------

# Create copy of input raster
forest_tws_15km_masked <- forest_tws_15km

# Loop through each raster and set the cells >50% outside to NA
for(i in 1:nlyr(forest_tws_15km_masked)){
  # Get the current layer
  current_layer <- forest_tws_15km[[i]]
  
  # Convert coordinates to cell indices
  cells_idx <- cellFromXY(current_layer, as.matrix(cells_to_mask[, c("x", "y")]))
  
  # Set these cells to NA
  current_layer[cells_idx] <- NA
  
  # Update the layer in the stack
  forest_tws_15km_masked[[i]] <- current_layer
}

## 5.3. Mask cells from TWS -> Forest ------------------------------------------
# Create copy of input raster
tws_forest_15km_masked <- tws_forest_15km

# Loop through each raster and set the cells >50% outside to NA
for(i in 1:nlyr(tws_forest_15km_masked)){
  # Get the current layer
  current_layer <- tws_forest_15km[[i]]
  
  # Convert coordinates to cell indices
  cells_idx <- cellFromXY(current_layer, as.matrix(cells_to_mask[, c("x", "y")]))
  
  # Set these cells to NA
  current_layer[cells_idx] <- NA
  
  # Update the layer in the stack
  tws_forest_15km_masked[[i]] <- current_layer
}

## 5.4. Mask cells from All -> Urban -------------------------------------------

# Create copy of input raster
all_urban_15km_masked <- all_urban_15km

# Loop through each raster and set the cells >50% outside to NA
for(i in 1:nlyr(all_urban_15km_masked)){
  # Get the current layer
  current_layer <- all_urban_15km[[i]]
  
  # Convert coordinates to cell indices
  cells_idx <- cellFromXY(current_layer, as.matrix(cells_to_mask[, c("x", "y")]))
  
  # Set these cells to NA
  current_layer[cells_idx] <- NA
  
  # Update the layer in the stack
  all_urban_15km_masked[[i]] <- current_layer
}

# 6. Validation checks ---------------------------------------------------------

## 6.1. Check that correct cells were removed ----------------------------------

# Compare numer of cells that should be removed with the cells that were actually removed
expecte_masked_cells <- nrow(cells_to_mask) #434

# Check the first layer or each raster
forest_masked_check <- is.na(forest_tws_15km_masked[[1]]) & !is.na(forest_tws_15km[[1]])
tws_masked_check <- is.na(tws_forest_15km_masked[[1]]) & !is.na(tws_forest_15km[[1]])
urban_masked_check <- is.na(all_urban_15km_masked[[1]]) & !is.na(all_urban_15km[[1]])

cat("Actual cells masked in forest_tws:", 
    global(forest_masked_check, "sum", na.rm=TRUE)[1,1], "\n") #426
cat("Actual cells masked in tws_forest:", 
    global(tws_masked_check, "sum", na.rm=TRUE)[1,1], "\n") #426
cat("Actual cells masked in all_urban:", 
    global(urban_masked_check, "sum", na.rm=TRUE)[1,1], "\n") #426 - 8 cells are not masked, why?

## 6.2. Check why the 8 cells are missing --------------------------------------

# Check if any cells_to_mask coordinates fall outside the raster extent
raster_extent <- ext(forest_tws_15km)
outside_extent <- cells_to_mask |>
  mutate(outside = x < raster_extent[1] | x > raster_extent[2] | 
           y < raster_extent[3] | y > raster_extent[4])
cat("Cells outside raster extent:", sum(outside_extent$outside), "\n") 
 # Cells outside raster extent: 0

# Check if any cells are already NA in the unmasked layer
tryCatch({
  already_na_forest <- check_na_values(forest_tws_15km[[1]], cells_to_mask)
  cat("Cells already NA in forest_tws:", already_na_forest, "\n")
}, error = function(e) {
  cat("Error checking forest_tws:", e$message, "\n")
}) # Cells already NA in forest_tws: 8

tryCatch({
  already_na_tws <- check_na_values(tws_forest_15km[[1]], cells_to_mask)
  cat("Cells already NA in tws_forest:", already_na_tws, "\n")
}, error = function(e) {
  cat("Error checking tws_forest:", e$message, "\n")
}) # Cells already NA in tws_forest: 8

tryCatch({
  already_na_urban <- check_na_values(all_urban_15km[[1]], cells_to_mask)
  cat("Cells already NA in all_urban:", already_na_urban, "\n")
}, error = function(e) {
  cat("Error checking all_urban:", e$message, "\n")
}) # Cells already NA in all_urban: 8


## 6.3. Visual check of masked cells -------------------------------------------

# Create a raster with the masked cells
masked_cell_raster <- forest_tws_15km_masked[[1]]

# Set all cells to 0
masked_cell_raster[] <- 0 

# Set the masked cells to 1
masked_cell_raster[forest_masked_check] <- 1

# Plot masked cells
ggplot() +
  geom_sf(data = norway_sf, fill = "lightgrey", color = "black") +
  geom_tile(data = as.data.frame(masked_cell_raster, xy = TRUE),
            aes(x = x, y = y, fill = factor('2000-2006_Non-forest'))) +
  geom_sf(data = edge_cells %>% filter(is_mostly_outside == TRUE),
          color = "red", size = 1, fill = NA) +
  scale_fill_manual(values = c("0" = "transparent", "1" = "blue",
                               name = "Status",
                               labels = c("0" = "Unchanged", "1" = "Masked"))) +
  theme_classic()

## 6.4. Check that unmasked cells are unchanged --------------------------------

# Compare values in the non-masked cells between original and masked rasters
forest_tws_check <- forest_tws_15km[[1]]
forest_tws_masked_check <- forest_tws_15km_masked[[1]]

# Create a mask for non-NA cells in both rasters
valid_cells <- !is.na(forest_tws_check) & !is.na(forest_tws_masked_check)

# Extract values for valid cells for both rasters
original_value <- forest_tws_check[valid_cells]
masked_values <- forest_tws_masked_check[valid_cells]

# Check if values are identical
all_values_match <- all(original_value == masked_values)
cat("All unmasked values remain unchanged:", all_values_match, "\n")

# If they don't all match, show statistics about the differences
if(!all_values_match){
  differences <- original_values - masked_values
  cat("Number of cells with different values:", sum(differences != 0), "\n")
  cat("Mean difference:", mean(differences), "\n")
  cat("Range of differences:", range(differences), "\n")
} 

## 6.5. Save masked layers to file ---------------------------------------------

# Forest -> TWS
writeRaster(forest_tws_15km_masked, 
            here("data", "derived_data", "clc_status_15km_forest_tws_masked.tif"),
            overwrite = TRUE)

# TWS -> Forest
writeRaster(tws_forest_15km_masked, 
            here("data", "derived_data", "clc_status_15km_tws_forest_masked.tif"),
            overwrite = TRUE)

# All -> Urban
writeRaster(all_urban_15km_masked, 
            here("data", "derived_data", "clc_status_15km_all_urban_combined_masked.tif"),
            overwrite = TRUE)

# END OF SCRIPT ----------------------------------------------------------------