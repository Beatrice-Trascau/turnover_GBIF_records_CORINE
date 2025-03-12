##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 3.1_turnover_land_cover
# This script contains code which calculates the temporal turnover of GBIF
# records in CORINE land cover pixels
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

## 1.1. 15km CORINE Status Layers ----------------------------------------------

# Forest -> TWS
forest_tws_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_forest_tws_masked.tif"))

# TWS -> Forest
tws_forest_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_tws_forest_masked.tif"))

# All -> Urban
all_urban_15km <- rast(here("data", "derived_data", 
                            "clc_status_15km_all_urban_combined_masked.tif"))


## 1.2. Cleaned GBIF Occurrences -----------------------------------------------

# Cleaned occurrence records
occurrences_norway <- fread(here("data", "derived_data", 
                                 "clean_occurrences_15km.txt"))

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
  mutate(period = case_when(year %in% 1997:2000 ~ "1997-2000",
                            year %in% 2006:2009 ~ "2006-2009"))

# Second period (2006-2012 LC change, 2003-2006 = before, 2012-2015 = after)
period2_occurrences <- occurrences_sf_reprojected |>
  filter(year %in% c(2003:2006, 2012:2015)) |>
  mutate(period = case_when(year %in% 2003:2006 ~ "2003-2006",
                            year %in% 2012:2015 ~ "2012-2015"))

# Third period (2012-2018 LC change, 2008-2012 = before, 2015-2018 = after)
period3_occurrences <- occurrences_sf_reprojected |>
  filter(year %in% c(2008:2012, 2015:2018)) |>
  mutate(period = case_when(year %in% 2008:2012 ~ "2008-2012",
                            year %in% 2015:2018 ~ "2015-2018"))

# Combine all periods
clean_occurrences_sf <- bind_rows(period1_occurrences,
                                  period2_occurrences,
                                  period3_occurrences)

# 4. PREP DATA FOR TURNOVER CALCULATION  ---------------------------------------

## 4.1. Extract cell IDs for occurrences ---------------------------------------

# This method retains spatial information for later
clean_occurrences_for_turnover <- clean_occurrences_sf |>
  mutate(cell = terra::extract(land_cover_id, 
                               st_coordinates(clean_occurrences_sf))[, 1])

## 4.2. Keep unique species records per cell and period ------------------------

# Make sure to only retain unique species records per cell and period
# Splitting these by time period to hopefully save computing power

# 2000-2006 period
period1_unique <- clean_occurrences_for_turnover |>
  filter(period %in% c("1997-2000", "2006-2009")) |>
  select(period, species, cell) |>
  distinct(cell, species, period)

# 2006-2012 period  
period2_unique <- clean_occurrences_for_turnover |>
  filter(period %in% c("2003-2006", "2012-2015")) |>
  select(period, species, cell) |>
  distinct(cell, species, period)

# 2012-2018 period
period3_unique <- clean_occurrences_for_turnover |>
  filter(period %in% c("2008-2012", "2015-2018")) |>
  select(period, species, cell) |>
  distinct(cell, species, period)

# Combine the results
unique_occurrences <- bind_rows(period1_unique,
                                period2_unique,
                                period3_unique)

# 5. CALCULATE TURNOVER --------------------------------------------------------

# Define period pairs to use in the function
period_pairs <- list(
  list(before = "1997-2000", after = "2006-2009", change_period = "2000-2006"),
  list(before = "2003-2006", after = "2012-2015", change_period = "2006-2012"),
  list(before = "2008-2012", after = "2015-2018", change_period = "2012-2018"))

# Apply turnover function to all periods
temporal_turnover <- map_df(period_pairs, function(pair) {
  calculate_jaccard_for_periods(
    unique_occurrences,
    pair$before,
    pair$after,
    pair$change_period)
})

# 6. COMBINE TURNOVER AND LAND COVER -------------------------------------------

# Add turnover results to the spatial dataframe and land cover
turnover_with_landcover <- temporal_turnover |>
  left_join(st_as_sf(land_cover_grid) |>
      mutate(x = st_coordinates(.)[,1],
             y = st_coordinates(.)[,2]),
      by = "cell") |>
  mutate(landcover_category = case_when(change_period == "2000-2006" ~ period_2000_2006,
                                        change_period == "2006-2012" ~ period_2006_2012,
                                        change_period == "2012-2018" ~ period_2012_2018))

# Save df 
saveRDS(turnover_with_landcover, 
        here("data", "derived_data", "full_turnover_results_100m.rds"))

# 7. QUICK EXPLORATION ---------------------------------------------------------
ggplot() +
  geom_boxplot(data = turnover_with_landcover,
               aes(x = landcover_category, y = jaccard_dissimilarity,
                   fill = "All data")) +
  geom_boxplot(data = turnover_3plus_species,
               aes(x = landcover_category, y = jaccard_dissimilarity,
                   fill = "3+ species")) +
  geom_boxplot(data = turnover_5plus_species,
               aes(x = landcover_category, y = jaccard_dissimilarity,
                   fill = "5+ species")) +
  facet_wrap(~change_period) +
  theme_bw() +
  labs(
    x = "Land Cover Category",
    y = "Jaccard Dissimilarity Index",
    title = "Species Turnover by Land Cover Change Category and Filtering",
    fill = "Dataset"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
