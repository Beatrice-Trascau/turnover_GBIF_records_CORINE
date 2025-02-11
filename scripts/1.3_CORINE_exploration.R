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
                           "clc_status_1km_all_urban_combined.tif"))

## 1.3. 5km resolution ---------------------------------------------------------

# Forest -> TWS
forest_tws_5km <- rast(here("data", "derived_data", 
                            "clc_status_5km_forest_tws.tif"))

# TWS -> Forest
tws_forest_5km <- rast(here("data", "derived_data", 
                            "clc_status_5km_tws_forest.tif"))

# All -> Urban
all_urban_5km <- rast(here("data", "derived_data", 
                           "clc_status_5km_all_urban_combined.tif"))

## 1.4. 15km resolution --------------------------------------------------------

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

# Create list with all transitions for each resolution
resolution_data <- list(
  "100m" = list(
    forest_tws = forest_tws_100m,
    tws_forest = tws_forest_100m,
    all_urban = all_urban_100m),
  "1km" = list(
    forest_tws = forest_tws_1km,
    tws_forest = tws_forest_1km,
    all_urban = all_urban_1km),
  "5km" = list(
    forest_tws = forest_tws_5km,
    tws_forest = tws_forest_5km,
    all_urban = all_urban_5km),
  "15km" = list(
    forest_tws = forest_tws_15km,
    tws_forest = tws_forest_15km,
    all_urban = all_urban_15km))

# Extract summary statistics tables for all resolutions at once
summary_tables <- process_all_resolutions(resolution_data)

# Check individual tables
summary_100m <- summary_tables[[1]]
summary_1km <- summary_tables[[2]]
summary_5km <- summary_tables[[3]]
summary_15km <- summary_tables[[4]]

# 3. MAP CHANGES ---------------------------------------------------------------

## 3.1. Plot at 100m resolution ------------------------------------------------

# Forest -> TWS panels
p1 <- create_transition_panel(forest_tws_100m[[1]], norway, 
                              "Forest to TWS", "2000-2006")
p2 <- create_transition_panel(forest_tws_100m[[2]], norway, 
                              "Forest to TWS", "2006-2012")
p3 <- create_transition_panel(forest_tws_100m[[3]], norway, 
                              "Forest to TWS", "2012-2018")

# TWS -> Forest
p4 <- create_transition_panel(tws_forest_100m[[1]], norway, "
                              TWS to Forest", "2000-2006")
p5 <- create_transition_panel(tws_forest_100m[[2]], norway,
                              "TWS to Forest", "2006-2012")
p6 <- create_transition_panel(tws_forest_100m[[3]], norway, 
                              "TWS to Forest", "2012-2018")

# All -> Urban 
p7 <- create_transition_panel(all_urban_100m[[1]], norway, 
                              "All to Urban", "2000-2006")
p8 <- create_transition_panel(all_urban_100m[[2]], norway, 
                              "All to Urban", "2006-2012")
p9 <- create_transition_panel(all_urban_100m[[3]], norway, 
                              "All to Urban", "2012-2018")

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

## 3.2. Plot at 1km resolution -------------------------------------------------

# Forest -> TWS
p1_1km <- create_aggregated_panel(forest_tws_1km[[3]], norway, 
                                  "Forest to TWS", "2000-2006", TRUE, TRUE)
p2_1km <- create_aggregated_panel(forest_tws_1km[[7]], norway, 
                                  "Forest to TWS", "2006-2012", FALSE, FALSE)
p3_1km <- create_aggregated_panel(forest_tws_1km[[11]], norway, 
                                  "Forest to TWS", "2012-2018", FALSE, FALSE)

# TWS -> Forest
p4_1km <- create_aggregated_panel(tws_forest_1km[[3]], norway, 
                                  "TWS to Forest", "2000-2006", TRUE, FALSE)
p5_1km <- create_aggregated_panel(tws_forest_1km[[7]], norway, 
                                  "TWS to Forest", "2006-2012", FALSE, FALSE)
p6_1km <- create_aggregated_panel(tws_forest_1km[[11]], norway, 
                                  "TWS to Forest", "2012-2018", FALSE, FALSE)

# All -> Urban panels
p7_1km <- create_aggregated_panel(all_urban_1km[[1]], norway, 
                                  "All to Urban", "2000-2006", TRUE, FALSE)
p8_1km <- create_aggregated_panel(all_urban_1km[[2]], norway, 
                                  "All to Urban", "2006-2012", FALSE, FALSE)
p9_1km <- create_aggregated_panel(all_urban_1km[[3]], norway, "
                                  All to Urban", "2012-2018", FALSE, FALSE)

# Combine 1km plots into specific rows
top_row_1km <- plot_grid(p1_1km, p4_1km, p7_1km, nrow = 1, align = 'h', 
                         rel_widths = c(1.2, 1, 1))
middle_row_1km <- plot_grid(p2_1km, p5_1km, p8_1km, nrow = 1, align = 'h')
bottom_row_1km <- plot_grid(p3_1km, p6_1km, p9_1km, nrow = 1, align = 'h')

# Combine all rows into 1
final_plot_1km <- plot_grid(top_row_1km, middle_row_1km, bottom_row_1km,
                            nrow = 3, align = 'v')

# Save 1km plot
ggsave(filename = here("figures", "Figure1_landcover_transitions_1km.png"),
       plot = final_plot_1km, width = 12, height = 12, dpi = 300)

## 3.3. Plot at 5km resolution -------------------------------------------------

# Forest -> TWS
p1_5km <- create_aggregated_panel(forest_tws_5km[[3]], norway, 
                                  "Forest to TWS", "2000-2006", TRUE, TRUE)
p2_5km <- create_aggregated_panel(forest_tws_5km[[7]], norway, 
                                  "Forest to TWS", "2006-2012", FALSE, FALSE)
p3_5km <- create_aggregated_panel(forest_tws_5km[[11]], norway, 
                                  "Forest to TWS", "2012-2018", FALSE, FALSE)

# TWS -> Forest
p4_5km <- create_aggregated_panel(tws_forest_5km[[3]], norway, 
                                  "TWS to Forest", "2000-2006", TRUE, FALSE)
p5_5km <- create_aggregated_panel(tws_forest_5km[[7]], norway, 
                                  "TWS to Forest", "2006-2012", FALSE, FALSE)
p6_5km <- create_aggregated_panel(tws_forest_5km[[11]], norway, 
                                  "TWS to Forest", "2012-2018", FALSE, FALSE)

# All -> Urban 
p7_5km <- create_aggregated_panel(all_urban_5km[[1]], norway, 
                                  "All to Urban", "2000-2006", TRUE, FALSE)
p8_5km <- create_aggregated_panel(all_urban_5km[[2]], norway, 
                                  "All to Urban", "2006-2012", FALSE, FALSE)
p9_5km <- create_aggregated_panel(all_urban_5km[[3]], norway, 
                                  "All to Urban", "2012-2018", FALSE, FALSE)

# Combine 5km plots into specifc rows
top_row_5km <- plot_grid(p1_5km, p4_5km, p7_5km, nrow = 1, align = 'h',
                         rel_widths = c(1.2, 1, 1))
middle_row_5km <- plot_grid(p2_5km, p5_5km, p8_5km, nrow = 1, align = 'h')
bottom_row_5km <- plot_grid(p3_5km, p6_5km, p9_5km, nrow = 1, align = 'h')

# Combine all rows into one
final_plot_5km <- plot_grid(top_row_5km, middle_row_5km, bottom_row_5km,
                            nrow = 3, align = 'v')

# Save 5km plot
ggsave(filename = here("figures", "Figure1_landcover_transitions_5km.png"),
       plot = final_plot_5km, width = 12, height = 12, dpi = 300)

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
