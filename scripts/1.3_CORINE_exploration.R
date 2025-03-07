##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 1.3_CORINE_exploration
# This script contains code which explores land cover transitions for the 
# 15km x 15km resolution
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

# 2. EXTRACT SUMMARY TABLES ----------------------------------------------------

# Frequency table for Forest -> TWS
forest_tws_freq <- freq(forest_tws_15km) |>
  mutate(transition = "Forest to TWS",
         layer_name = names(forest_tws_15km)[layer])
    
# Frequency table for TWS -> Forest
tws_forest_freq <- freq(tws_forest_15km) |>
  mutate(transition = "TWS to Forest",
         layer_name = names(tws_forest_15km)[layer])

# Frequency table for all to urban 
all_urban_freq <- freq(all_urban_15km) |>
  mutate(transition = "All to Urban", 
         layer_name = names(all_urban_15km)[layer])

# Combine all frequencies in a single dataframe
summary_15km  <- rbind(forest_tws_freq, tws_forest_freq, all_urban_freq)

# Clean up the dataframe 
summary_15km_clean <- summary_15km |>
  rename(small_pixel_number = value) |>
  mutate(area = small_pixel_number * 0.01,
         resolution = "15km") 

# Get summary table
summary_15km_table <- summary_15km_clean |>
  kable(col.names = c("layer", "small_pixel_number", "count", 
                      "transition", "layer_name", "area", "resolution")) |>
  kable_styling(bootstrap_options = c("striped", "hover"))

# Print summary
print(summary_15km_table)


# 3. MAP CHANGES ---------------------------------------------------------------

# 3.1. Convert rasters to df for plotting --------------------------------------

# Define the list of rasters to be converted
raster_list <- list("forest_tws_3" = forest_tws_15km[[3]],   # Forest -> TWS 2000-2006
                    "forest_tws_7" = forest_tws_15km[[7]],   # Forest -> TWS 2006-2012
                    "forest_tws_11" = forest_tws_15km[[11]], # Forest -> TWS 2012-2018
                    "tws_forest_3" = tws_forest_15km[[3]],   # TWS -> Forest 2000-2006
                    "tws_forest_7" = tws_forest_15km[[7]],   # TWS -> Forest 2006-2012
                    "tws_forest_11" = tws_forest_15km[[11]], # TWS -> Forest 2012-2018
                    "all_urban_2" = all_urban_15km[[2]],     # All -> Urban 2000-2006
                    "all_urban_5" = all_urban_15km[[5]],     # All -> Urban 2006-2012
                    "all_urban_8" = all_urban_15km[[8]])      # All -> Urban 2012-2018

# Define list to store dfs
raster_dfs <- list()

# Loop through rasters and create dfs
for (name in names(raster_list)){
  # Convert raster to df
  df <- as.data.frame(raster_list[[name]], xy = TRUE)
  
  # Rename third column to "value"
  names(df)[3] <- "value"
  
  # Set 0 values to NA for better map visualisation
  df$value[df$value == 0] <- NA
  
  # Calculate % relative to total possible pixels in 15km cells = 150 * 150 = 22 500
  df$percent <- (df$value / 22500) * 100
  
  # Store dfs in the list created earlier
  raster_dfs[[name]] <- df
}
  
## 3.2. F -> TWS panels --------------------------------------------------------

