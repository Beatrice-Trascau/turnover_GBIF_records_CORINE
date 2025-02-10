##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 1.3_CORINE_exploration
# This script contains code which explores land cover transitions for each
# grain size considered
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

## 1.1. 100m resolution --------------------------------------------------------

# Forest -> TWS
forest_tws_100m <- rast(here("data", "derived_data", 
                             "clc_status_100m_forest_tws.tif"))

# TWS -> Forest
tws_forest_100m <- rast(here("data", "derived_data", 
                             "clc_status_100m_tws_forest.tif"))

# All -> Urban
all_urban_100m <- rast(here("data", "derived_data", 
                            "clc_status_100m_all_urban.tif"))

## 1.2. 1km resolution ---------------------------------------------------------

# Forest -> TWS
forest_tws_1km <- rast(here("data", "derived_data", 
                            "clc_status_1km_forest_tws.tif"))

# TWS -> Forest
tws_forest_1km <- rast(here("data", "derived_data", 
                            "clc_status_1km_tws_forest.tif"))

# All -> Urban
all_urban_1km <- rast(here("data", "derived_data", 
                           "clc_status_1km_all_urban.tif"))

## 1.3. 5km resolution ---------------------------------------------------------

# Forest -> TWS
forest_tws_5km <- rast(here("data", "derived_data", 
                            "clc_status_5km_forest_tws.tif"))

# TWS -> Forest
tws_forest_5km <- rast(here("data", "derived_data", 
                            "clc_status_5km_tws_forest.tif"))

# All -> Urban
all_urban_5km <- rast(here("data", "derived_data", 
                           "clc_status_5km_all_urban.tif"))

## 1.4. 15km resolution --------------------------------------------------------

# Forest -> TWS
forest_tws_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_forest_tws.tif"))

# TWS -> Forest
tws_forest_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_tws_forest.tif"))

# All -> Urban
all_urban_15km <- rast(here("data", "derived_data", 
                            "clc_status_15km_all_urban.tif"))

# Load Norway shapefile
norway <- vect(here("data", "derived_data", "reprojected_norway_shapefile", 
                    "norway_corine_projection.shp"))

# 2. EXTRACT SUMMARY TABLES ----------------------------------------------------

## 2.1. 100m resolution --------------------------------------------------------

# Extract frequency tables
forest_tws_100m_freq <- freq(forest_tws_100m)
tws_forest_100m_freq <- freq(tws_forest_100m)
all_urban_100m_freq <- freq(all_urban_100m)

# Add column to distinguish transition types
forest_tws_100m_freq$transition <- "Forest to TWS"
tws_forest_100m_freq$transition <- "TWS to Forest"
all_urban_100m_freq$transition <- "All to Urban"

# Bind all 3 into single df
summary_df_100m <- rbind(forest_tws_100m_freq,
                    tws_forest_100m_freq,
                    all_urban_100m_freq)
# Clean it up 
summary_df_100m_clean <- summary_df_100m |>
  mutate(layer = case_when(layer == 1 ~ "2000-2006",
                           layer == 2 ~ "2006-2012",
                           layer == 3 ~ "2012-2018"),
         area = count * 0.01)

# Create table from the data
summary_df_100m_clean |>
  kable(col.names = c("layer", "value", "count", "transition", "area")) |>
  kable_styling(bootstrap_options = c("striped", "hover"))
  
## 2.2. 1km resolution ---------------------------------------------------------

# Forest -> TWS
forest_tws_1km_summary <- create_summary_table(forest_tws_1km, 1000, 
                                               "Forest to TWS")

# TWS -> Forests
tws_forest_1km_summary <- create_summary_table(tws_forest_1km, 1000, 
                                               "TWS to Forest")

# All -> Urban
all_urban_1km_summary <- create_summary_table(all_urban_1km, 1000, 
                                              "All to Urban")

## 2.3. 5km resolution ---------------------------------------------------------

# Forest -> TWS
forest_tws_5km_summary <- create_summary_table(forest_tws_5km, 5000, 
                                               "Forest to TWS")

# TWS -> Forest
tws_forest_5km_summary <- create_summary_table(tws_forest_5km, 5000, 
                                               "TWS to Forest")

# All -> Urban
all_urban_5km_summary <- create_summary_table(all_urban_5km, 5000, 
                                              "All to Urban")

## 2.4. 15km resolution --------------------------------------------------------

# Forest -> TWS
forest_tws_15km_summary <- create_summary_table(forest_tws_15km, 15000, 
                                                "Forest to TWS")

# TWS -> Forest
tws_forest_15km_summary <- create_summary_table(tws_forest_15km, 15000, 
                                                "TWS to Forest")

# All -> Urban
all_urban_15km_summary <- create_summary_table(all_urban_15km, 15000, 
                                               "All to Urban")

## 2.5. Export table -----------------------------------------------------------

# Combine all summaries
all_summaries <- bind_rows(
  forest_tws_100m_summary, tws_forest_100m_summary, all_urban_100m_summary,
  forest_tws_1km_summary, tws_forest_1km_summary, all_urban_1km_summary,
  forest_tws_5km_summary, tws_forest_5km_summary, all_urban_5km_summary,
  forest_tws_15km_summary, tws_forest_15km_summary, all_urban_15km_summary)

# Export table
