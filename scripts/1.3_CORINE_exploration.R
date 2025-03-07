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

# Load Norway shapefile
norway <- vect(here("data", "derived_data", "reprojected_norway_shapefile", 
                    "norway_corine_projection.shp"))

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


## 3.4. Plot at 15km resolution ------------------------------------------------

# Forest -> TWS 
p1_15km <- create_aggregated_panel(forest_tws_15km[[3]], norway, 
                                   "Forest to TWS", "2000-2006", TRUE, TRUE)
p2_15km <- create_aggregated_panel(forest_tws_15km[[7]], norway, 
                                   "Forest to TWS", "2006-2012", FALSE, FALSE)
p3_15km <- create_aggregated_panel(forest_tws_15km[[11]], norway,
                                   "Forest to TWS", "2012-2018", FALSE, FALSE)

# TWS -> Forest
p4_15km <- create_aggregated_panel(tws_forest_15km[[3]], norway, 
                                   "TWS to Forest", "2000-2006", TRUE, FALSE)
p5_15km <- create_aggregated_panel(tws_forest_15km[[7]], norway, 
                                   "TWS to Forest", "2006-2012", FALSE, FALSE)
p6_15km <- create_aggregated_panel(tws_forest_15km[[11]], norway, 
                                   "TWS to Forest", "2012-2018", FALSE, FALSE)

# All -> Urban 
p7_15km <- create_aggregated_panel(all_urban_15km[[1]], norway, 
                                   "All to Urban", "2000-2006", TRUE, FALSE)
p8_15km <- create_aggregated_panel(all_urban_15km[[2]], norway, 
                                   "All to Urban", "2006-2012", FALSE, FALSE)
p9_15km <- create_aggregated_panel(all_urban_15km[[3]], norway, 
                                   "All to Urban", "2012-2018", FALSE, FALSE)

# Combine 15km plots into specific rows
top_row_15km <- plot_grid(p1_15km, p4_15km, p7_15km, nrow = 1, align = 'h',
                          rel_widths = c(1.2, 1, 1))
middle_row_15km <- plot_grid(p2_15km, p5_15km, p8_15km, nrow = 1, align = 'h')
bottom_row_15km <- plot_grid(p3_15km, p6_15km, p9_15km, nrow = 1, align = 'h')

# Combine all rows into 1
final_plot_15km <- plot_grid(top_row_15km, middle_row_15km, bottom_row_15km,
                             nrow = 3, align = 'v')

# Save 15km plot
ggsave(filename = here("figures", "Figure1_landcover_transitions_15km.png"),
       plot = final_plot_15km, width = 12, height = 12, dpi = 300)
