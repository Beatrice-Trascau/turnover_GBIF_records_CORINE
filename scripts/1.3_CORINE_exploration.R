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

## 3.1. Convert rasters to df for plotting -------------------------------------

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

# Find global min and max for unified legend scale
all_percent_values <- unlist(lapply(raster_dfs, function(df) df$percent))
min_percent <- min(all_percent_values, na.rm = TRUE)
max_percent <- max(all_percent_values, na.rm = TRUE)

# Round values up and down to 2 decimal places
min_percent_rounded <- ceiling(min_percent * 100) / 100 
max_percent_rounded <- ceiling(max_percent * 100) / 100  # Round up

# Create common scale for all plots
common_fill_scale <- scale_fill_viridis_c(option = "virdis",
                                         name = "% of grid cell\narea changed",
                                         limits = c(min_percent_rounded,
                                                    max_percent_rounded),
                                         breaks = seq(min_percent_rounded, 
                                                      max_percent_rounded, 
                                                      length.out = 5),
                                         labels = function(x) sprintf("%.2f%%", x),
                                         na.value = "white",
                                         guide = guide_colorbar(direction = "vertical",
                                                                barheight = unit(10, "lines"),
                                                                barwidth = unit(1, "lines"),
                                                                title.position = "top",
                                                                title.hjust = 0.5,
                                                                frame.colour = "black",
                                                                ticks.colour = "black"))

  
## 3.2. F -> TWS panels --------------------------------------------------------

# F -> TWS 2000-2006
p1 <- ggplot() +
  geom_tile(data = raster_dfs[["forest_tws_3"]], 
            aes(x = x, y = y, fill = percent)) +
  common_fill_scale +
  geom_sf(data = st_as_sf(norway), fill = NA, color = "black", linewidth = 0.2) +
  coord_sf() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_y = unit(0.8, "cm"), 
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "br", width_hint = 0.35)

# F -> TWS 2006-2012
p2 <- ggplot() +
  geom_tile(data = raster_dfs[["forest_tws_7"]], 
            aes(x = x, y = y, fill = percent)) +
  common_fill_scale +
  guides(fill = "none") +
  geom_sf(data = st_as_sf(norway), fill = NA, color = "black", linewidth = 0.2) +
  coord_sf() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

# F -> TWS 2012-2018
p3 <- ggplot() +
  geom_tile(data = raster_dfs[["forest_tws_11"]], 
            aes(x = x, y = y, fill = percent)) +
  common_fill_scale +
  guides(fill = "none") +
  geom_sf(data = st_as_sf(norway), fill = NA, color = "black", linewidth = 0.2) +
  coord_sf() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

## 3.3. TWS -> F panels --------------------------------------------------------

# TWS -> F 2000-2006
p4 <- ggplot() +
  geom_tile(data = raster_dfs[["tws_forest_3"]], 
            aes(x = x, y = y, fill = percent)) +
  common_fill_scale +
  guides(fill = "none") +
  geom_sf(data = st_as_sf(norway), fill = NA, color = "black", linewidth = 0.2) +
  coord_sf() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

# TWS -> F 2006-2012
p5 <- ggplot() +
  geom_tile(data = raster_dfs[["tws_forest_7"]], 
            aes(x = x, y = y, fill = percent)) +
  common_fill_scale +
  guides(fill = "none") +
  geom_sf(data = st_as_sf(norway), fill = NA, color = "black", linewidth = 0.2) +
  coord_sf() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

# TWS -> F 2012-2018
p6 <- ggplot() +
  geom_tile(data = raster_dfs[["tws_forest_11"]], 
            aes(x = x, y = y, fill = percent)) +
  common_fill_scale +
  guides(fill = "none") +
  geom_sf(data = st_as_sf(norway), fill = NA, color = "black", linewidth = 0.2) +
  coord_sf() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

## 3.4. All -> Artifical Surfaces panels ---------------------------------------

# All -> AF 2000-2006
p7 <- ggplot() +
  geom_tile(data = raster_dfs[["all_urban_2"]], 
            aes(x = x, y = y, fill = percent)) +
  common_fill_scale +
  guides(fill = "none") +
  geom_sf(data = st_as_sf(norway), fill = NA, color = "black", linewidth = 0.2) +
  coord_sf() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

# All -> AF 2006-2012
p8 <- ggplot() +
  geom_tile(data = raster_dfs[["all_urban_5"]], 
            aes(x = x, y = y, fill = percent)) +
  common_fill_scale +
  guides(fill = "none") +
  geom_sf(data = st_as_sf(norway), fill = NA, color = "black", linewidth = 0.2) +
  coord_sf() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

# All -> AF 2012-2018
p9 <- ggplot() +
  geom_tile(data = raster_dfs[["all_urban_8"]], 
            aes(x = x, y = y, fill = percent)) +
  common_fill_scale +
  guides(fill = "none") +
  geom_sf(data = st_as_sf(norway), fill = NA, color = "black", linewidth = 0.2) +
  coord_sf() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

## 3.5. Combine plots into single figure ---------------------------------------

# First row
top_row <- plot_grid(p1, p2, p3, nrow = 1, 
                     labels = c('a)', 'b)', 'c)'),
                     align = "h", rel_widths = c(1.3, 1, 1))

# Second row
middle_row <- plot_grid(p4, p5, p6, 
                        labels = c('d)', 'e)', 'f)'),
                        nrow = 1, align = "h")


# Third row
bottom_row <- plot_grid(p7, p8, p9, 
                        labels = c('g)', 'h)', 'i)'),
                        nrow = 1, align = "h")

# Combine all rows into a single figure
final_map <- plot_grid(top_row, middle_row, bottom_row,
                       nrow = 3, align = "v")

# Save figure to file
ggsave(filename = here("figures", "Figure1_landcover_transitions_15km.png"),
       plot = final_map, width = 12, height = 12, dpi = 300)

ggsave(filename = here("figures", "Figure1_landcover_transitions_15km.svg"),
       plot = final_map, width = 12, height = 12, dpi = 300)

# 4. EXPLORATORY FIGURES -------------------------------------------------------

# Calculate total area affected by each transition type per time period
transition_summary <- summary_15km_clean |>
  group_by(transition, layer_name) |>
  summarize(
    total_area_km2 = sum(area),
    .groups = 'drop'
  ) |>
  arrange(transition, layer_name)

# Print transition summary
print(transition_summary)

# Create a bar plot of transition areas by time period
transition_plot <- ggplot(transition_summary, 
                          aes(x = layer_name, y = total_area_km2, fill = transition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d() +
  labs(x = "Time Period", y = "Total Area (km²)", 
       title = "Land Cover Transitions in Norway (15km resolution)") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Save the transition summary plot
ggsave(filename = here("figures", "Figure2_transition_summary_15km.png"),
       plot = transition_plot, width = 10, height = 6, dpi = 300)

# END OF SCRIPT ----------------------------------------------------------------