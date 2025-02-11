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

# Extract frequency tables
forest_tws_1km_freq <- freq(forest_tws_1km)
tws_forest_1km_freq <- freq(tws_forest_1km)
all_urban_1km_freq <- freq(all_urban_1km)

# Add layer names to the frequency tables
forest_tws_1km_freq <- forest_tws_1km_freq |>
  mutate(transition = "Forest to TWS",
         layer_name = names(forest_tws_1km)[forest_tws_1km_freq$layer])

tws_forest_1km_freq <- tws_forest_1km_freq |>
  mutate(transition = "TWS to Forest",
         layer_name = names(tws_forest_1km)[tws_forest_1km_freq$layer])

all_urban_1km_freq <- all_urban_1km_freq |>
  mutate(transition = "All to Urban",
         layer_name = names(all_urban_1km)[all_urban_1km_freq$layer])

# Bind all 3 into single df
summary_df_1km <- rbind(forest_tws_1km_freq,
                         tws_forest_1km_freq,
                         all_urban_1km_freq)
# Clean it up 
summary_df_1km_clean <- summary_df_1km |>
  mutate(area = count * 0.01) |>
  rename(small_pixel_number = value)

# Create table from the data
summary_df_1km_clean |>
  kable(col.names = c("layer", "small_pixel_number", "count", "transition", 
                      "layer_name", "area")) |>
  kable_styling(bootstrap_options = c("striped", "hover"))


## 2.3. 5km resolution ---------------------------------------------------------

# Extract frequency tables
forest_tws_5km_freq <- freq(forest_tws_5km)
tws_forest_5km_freq <- freq(tws_forest_5km)
all_urban_5km_freq <- freq(all_urban_5km)

# Add layer names to the frequency tables
forest_tws_5km_freq <- forest_tws_5km_freq |>
  mutate(transition = "Forest to TWS",
         layer_name = names(forest_tws_5km)[forest_tws_5km_freq$layer])

tws_forest_5km_freq <- tws_forest_5km_freq |>
  mutate(transition = "TWS to Forest",
         layer_name = names(tws_forest_5km)[tws_forest_5km_freq$layer])

all_urban_5km_freq <- all_urban_5km_freq |>
  mutate(transition = "All to Urban",
         layer_name = names(all_urban_5km)[all_urban_5km_freq$layer])

# Bind all 3 into single df
summary_df_5km <- rbind(forest_tws_5km_freq,
                        tws_forest_5km_freq,
                        all_urban_5km_freq)
# Clean it up 
summary_df_5km_clean <- summary_df_5km |>
  mutate(area = count * 0.01) |>
  rename(small_pixel_number = value)

# Create table from the data
summary_df_5km_clean |>
  kable(col.names = c("layer", "small_pixel_number", "count", "transition", 
                      "layer_name", "area")) |>
  kable_styling(bootstrap_options = c("striped", "hover"))

## 2.4. 15km resolution --------------------------------------------------------

# Extract frequency tables
forest_tws_15km_freq <- freq(forest_tws_15km)
tws_forest_15km_freq <- freq(tws_forest_15km)
all_urban_15km_freq <- freq(all_urban_15km)

# Add layer names to the frequency tables
forest_tws_15km_freq <- forest_tws_15km_freq |>
  mutate(transition = "Forest to TWS",
         layer_name = names(forest_tws_15km)[forest_tws_15km_freq$layer])

tws_forest_15km_freq <- tws_forest_15km_freq |>
  mutate(transition = "TWS to Forest",
         layer_name = names(tws_forest_15km)[tws_forest_15km_freq$layer])

all_urban_15km_freq <- all_urban_15km_freq |>
  mutate(transition = "All to Urban",
         layer_name = names(all_urban_15km)[all_urban_15km_freq$layer])

# Bind all 3 into single df
summary_df_15km <- rbind(forest_tws_15km_freq,
                        tws_forest_15km_freq,
                        all_urban_15km_freq)
# Clean it up 
summary_df_15km_clean <- summary_df_15km |>
  mutate(area = count * 0.01) |>
  rename(small_pixel_number = value)

# Create table from the data
summary_df_15km_clean |>
  kable(col.names = c("layer", "small_pixel_number", "count", "transition", 
                      "layer_name", "area")) |>
  kable_styling(bootstrap_options = c("striped", "hover"))

# 3. MAP CHANGES ---------------------------------------------------------------

## 3.1. Plot at 100m resolution ------------------------------------------------

# Function to create a single panel
create_transition_panel <- function(raster_layer, norway, transition_type, time_period, filter_value) {
  # Convert raster to dataframe
  raster_df <- as.data.frame(raster_layer, xy = TRUE)
  names(raster_df)[3] <- "value"
  
  # Filter based on transition type
  filtered_df <- switch(transition_type,
                        "Forest to TWS" = filter(raster_df, value == "Forest to TWS"),
                        "TWS to Forest" = filter(raster_df, value == "TWS to Forest"),
                        "All to Urban" = filter(raster_df, 
                                                !value %in% c("Urban no change", "No urban conversion"))
  )
  
  # Create base plot
  p <- ggplot() +
    geom_sf(data = st_as_sf(norway), fill = "white", color = "gray50", linewidth = 0.2) +
    geom_point(data = filtered_df, aes(x = x, y = y), 
               size = 0.1, color = "red", alpha = 0.5) +
    coord_sf() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 8, hjust = 0.5)
    ) +
    ggtitle(paste0(time_period, "\n", transition_type))
  
  return(p)
}

# Create all panels
# Forest to TWS panels
p1 <- create_transition_panel(forest_tws_100m[[1]], norway, "Forest to TWS", "2000-2006")
p2 <- create_transition_panel(forest_tws_100m[[2]], norway, "Forest to TWS", "2006-2012")
p3 <- create_transition_panel(forest_tws_100m[[3]], norway, "Forest to TWS", "2012-2018")

# TWS to Forest panels
p4 <- create_transition_panel(tws_forest_100m[[1]], norway, "TWS to Forest", "2000-2006")
p5 <- create_transition_panel(tws_forest_100m[[2]], norway, "TWS to Forest", "2006-2012")
p6 <- create_transition_panel(tws_forest_100m[[3]], norway, "TWS to Forest", "2012-2018")

# All to Urban panels
p7 <- create_transition_panel(all_urban_100m[[1]], norway, "All to Urban", "2000-2006")
p8 <- create_transition_panel(all_urban_100m[[2]], norway, "All to Urban", "2006-2012")
p9 <- create_transition_panel(all_urban_100m[[3]], norway, "All to Urban", "2012-2018")

# Combine all panels
combined_plot <- (p1 + p4 + p7) / 
  (p2 + p5 + p8) / 
  (p3 + p6 + p9) +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 12, hjust = 0.5)))

# Save the plot
ggsave(
  filename = here("figures", "Figure1_landcover_transitions_100m.png"),
  plot = combined_plot,
  width = 12,
  height = 12,
  dpi = 300)

