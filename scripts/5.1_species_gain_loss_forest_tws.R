##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 5.1_species_gain_loss_forest_tws
# This script calculates and visualises species gain and loss for birds and
# vascular plants in areas with forest to TWS land cover changes
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# Turnover data
load(here("data", "derived_data", 
          "vascular_plants_all_periods_turnover_all_land_cover_chanegs_15km.rda"))
load(here("data", "derived_data", 
          "birds_all_periods_turnover_all_land_cover_chanegs_15km.rda"))

# Rename for easier handling
plants_turnover <- vascular_plants_all_periods_turnover_all_land_cover_chanegs_15km
birds_turnover <- birds_all_periods_turnover_all_land_cover_chanegs_15km

# Forest -> TWS raster
forest_tws_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_forest_tws_masked.tif"))

# TWS -> Forest raster
tws_forest_15km <- rast(here("data", "derived_data", 
                             "clc_status_15km_tws_forest_masked.tif"))

# All -> Urban raster
all_urban_15km <- rast(here("data", "derived_data", 
                            "clc_status_15km_all_urban_combined_masked.tif"))

# 2. PREPARE DATA FOR ANALYSIS -------------------------------------------------

## 2.1. Prepare vascular plants data -------------------------------------------

# Process all land cover change information
plants_all_changes <- plants_turnover |>
  # add species change metrics
  mutate(species_change = total_spp_after - total_spp_before,
         change_type = case_when(species_change > 0 ~ "Gain",
                                 species_change < 0 ~ "Loss",
                                 TRUE ~ "No Change"),
         abs_change = abs(species_change),
         taxonomic_group = "Vascular Plants") |>
  # extract land cover change information based on time period
  mutate(
    # forest to TWS transitions
    forest_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest no change`,
                                 lc_time_period == "2006-2012" ~ `2006-2012_Forest no change`,
                                 lc_time_period == "2012-2018" ~ `2012-2018_Forest no change`,
                                 TRUE ~ NA_real_),
    forest_to_tws = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest to TWS`,
                              lc_time_period == "2006-2012" ~ `2006-2012_Forest to TWS`,
                              lc_time_period == "2012-2018" ~ `2012-2018_Forest to TWS`,
                              TRUE ~ NA_real_),
    # TWS to Forest transitions
    tws_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS no change`,
                              lc_time_period == "2006-2012" ~ `2006-2012_TWS no change`,
                              lc_time_period == "2012-2018" ~ `2012-2018_TWS no change`,
                              TRUE ~ NA_real_),
    tws_to_forest = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS to Forest`,
                              lc_time_period == "2006-2012" ~ `2006-2012_TWS to Forest`,
                              lc_time_period == "2012-2018" ~ `2012-2018_TWS to Forest`,
                              TRUE ~ NA_real_),
    # all to Urban transitions
    urban_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Urban_no_change`,
                                lc_time_period == "2006-2012" ~ `2006-2012_Urban_no_change`,
                                lc_time_period == "2012-2018" ~ `2012-2018_Urban_no_change`,
                                TRUE ~ NA_real_),
    all_to_urban = case_when(lc_time_period == "2000-2006" ~ `2000-2006_all_to_urban`,
                             lc_time_period == "2006-2012" ~ `2006-2012_all_to_urban`,
                             lc_time_period == "2012-2018" ~ `2012-2018_all_to_urban`,
                             TRUE ~ NA_real_)) |>
  # create land cover change flags for filtering
  mutate(has_forest_to_tws = !is.na(forest_to_tws) & forest_to_tws > 0,
         has_tws_to_forest = !is.na(tws_to_forest) & tws_to_forest > 0,
         has_all_to_urban = !is.na(all_to_urban) & all_to_urban > 0) |>
  # Select relevant columns
  select(cell_ID, lc_time_period, x, y, taxonomic_group, species_change, 
         change_type, abs_change, total_spp_before, total_spp_after, recorder_effort,
         forest_no_change, forest_to_tws, has_forest_to_tws,
         tws_no_change, tws_to_forest, has_tws_to_forest,
         urban_no_change, all_to_urban, has_all_to_urban)

## 2.2 Prepare birds data ------------------------------------------------------

birds_all_changes <- birds_turnover |>
  # Add species change metrics
  mutate(species_change = total_spp_after - total_spp_before,
         change_type = case_when(species_change > 0 ~ "Gain",
                                 species_change < 0 ~ "Loss",
                                 TRUE ~ "No Change"),
         abs_change = abs(species_change),
         taxonomic_group = "Birds") |>
  # extract land cover change information based on time period
  mutate(
    # forest to TWS transitions
    forest_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest no change`,
                                 lc_time_period == "2006-2012" ~ `2006-2012_Forest no change`,
                                 lc_time_period == "2012-2018" ~ `2012-2018_Forest no change`,
                                 TRUE ~ NA_real_),
    forest_to_tws = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest to TWS`,
                              lc_time_period == "2006-2012" ~ `2006-2012_Forest to TWS`,
                              lc_time_period == "2012-2018" ~ `2012-2018_Forest to TWS`,
                              TRUE ~ NA_real_),
    # TWS to Forest transitions
    tws_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS no change`,
                              lc_time_period == "2006-2012" ~ `2006-2012_TWS no change`,
                              lc_time_period == "2012-2018" ~ `2012-2018_TWS no change`,
                              TRUE ~ NA_real_),
    tws_to_forest = case_when(lc_time_period == "2000-2006" ~ `2000-2006_TWS to Forest`,
                              lc_time_period == "2006-2012" ~ `2006-2012_TWS to Forest`,
                              lc_time_period == "2012-2018" ~ `2012-2018_TWS to Forest`,
                              TRUE ~ NA_real_),
    # all to Urban transitions
    urban_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Urban_no_change`,
                                lc_time_period == "2006-2012" ~ `2006-2012_Urban_no_change`,
                                lc_time_period == "2012-2018" ~ `2012-2018_Urban_no_change`,
                                TRUE ~ NA_real_),
    all_to_urban = case_when(lc_time_period == "2000-2006" ~ `2000-2006_all_to_urban`,
                             lc_time_period == "2006-2012" ~ `2006-2012_all_to_urban`,
                             lc_time_period == "2012-2018" ~ `2012-2018_all_to_urban`,
                             TRUE ~ NA_real_)) |>
  # create land cover change flags for filtering
  mutate(has_forest_to_tws = !is.na(forest_to_tws) & forest_to_tws > 0,
         has_tws_to_forest = !is.na(tws_to_forest) & tws_to_forest > 0,
         has_all_to_urban = !is.na(all_to_urban) & all_to_urban > 0) |>
  # select relevant columns
  select(cell_ID, lc_time_period, x, y, taxonomic_group, species_change, 
         change_type, abs_change, total_spp_before, total_spp_after, recorder_effort,
         forest_no_change, forest_to_tws, has_forest_to_tws,
         tws_no_change, tws_to_forest, has_tws_to_forest,
         urban_no_change, all_to_urban, has_all_to_urban)

# Combine dataframes into 1
all_species_changes <- bind_rows(plants_all_changes, birds_all_changes)

# 3. GET SUMMARY STATISTICS ----------------------------------------------------

## 3.1. Forest -> TWS ----------------------------------------------------------

# Filter for cells with Forest to TWS transitions
forest_tws_data <- all_species_changes |> 
  filter(has_forest_to_tws)

# Create summary statistics table
summary_stats <- forest_tws_data |>
  group_by(lc_time_period, taxonomic_group, change_type) |>
  summarize(count = n(),
            mean_change = mean(abs_change, na.rm = TRUE),
            median_change = median(abs_change, na.rm = TRUE),
            max_change = max(abs_change, na.rm = TRUE),
            .groups = "drop") |>
  # calculate percentage of cells with each change type
  group_by(lc_time_period, taxonomic_group) |>
  mutate(total_cells = sum(count),
         percentage = count / total_cells * 100) |>
  ungroup()

# Print summary table
#print(summary_stats)

## 3.2. TWS -> Forest ----------------------------------------------------------

# Filter for cells with TWS to Forest transitions
tws_forest_data <- all_species_changes |> 
  filter(has_tws_to_forest)

# Create summary statistics table
summary_stats_tws_forest <- tws_forest_data |>
  group_by(lc_time_period, taxonomic_group, change_type) |>
  summarize(count = n(),
            mean_change = mean(abs_change, na.rm = TRUE),
            median_change = median(abs_change, na.rm = TRUE),
            max_change = max(abs_change, na.rm = TRUE),
            .groups = "drop") |>
  # calculate percentage of cells with each change type
  group_by(lc_time_period, taxonomic_group) |>
  mutate(total_cells = sum(count),
         percentage = count / total_cells * 100) |>
  ungroup()

## 3.3. All -> Urban -----------------------------------------------------------

# Filter for cells with TWS to Forest transitions
all_urban_data <- all_species_changes |> 
  filter(has_all_to_urban)

# Create summary statistics table
summary_stats_all_urban <- all_urban_data |>
  group_by(lc_time_period, taxonomic_group, change_type) |>
  summarize(count = n(),
            mean_change = mean(abs_change, na.rm = TRUE),
            median_change = median(abs_change, na.rm = TRUE),
            max_change = max(abs_change, na.rm = TRUE),
            .groups = "drop") |>
  # calculate percentage of cells with each change type
  group_by(lc_time_period, taxonomic_group) |>
  mutate(total_cells = sum(count),
         percentage = count / total_cells * 100) |>
  ungroup()


# 4. VISUALISE GAIN AND LOSS ---------------------------------------------------

# Define time period order
time_period_order <- c("2000-2006", "2006-2012", "2012-2018")

# Define colors for gain and loss
gain_loss_colors <- c("Gain" = "#1a9850", "Loss" = "#d73027", "No Change" = "#4575b4")

## 4.1. Forest -> TWS ----------------------------------------------------------

# Jitter plot of values
jitter_species_change_plot <- ggplot(forest_tws_data, 
                                     aes(x = species_change, y = taxonomic_group)) +
  # add vertical line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  # add jittered points with different shapes for gain and loss
  geom_jitter(aes(color = change_type, shape = change_type, size = forest_to_tws),
              height = 0.3, alpha = 0.7) +
  # set shapes for gain and loss
  scale_shape_manual(values = c("Gain" = 16, "Loss" = 17, "No Change" = 15)) +
  # set colors for gain and loss
  scale_color_manual(values = gain_loss_colors) +
  # set size scale
  scale_size_continuous(name = "Forest → TWS\nTransition Pixels", 
                        range = c(2, 8)) +
  # facet by time period
  facet_wrap(~ factor(lc_time_period, levels = time_period_order), 
             ncol = 3) +
  # labels and theme
  labs(y = "Group",
       color = "Change Type",
       shape = "Change Type") +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "gray80"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        legend.box = "vertical",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(angle=90))

# Density plot for distribution of species change
density_plot <- ggplot(forest_tws_data, 
                       aes(x = species_change, fill = taxonomic_group)) +
  # add vertical line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  # add density curves
  geom_density(alpha = 0.7) +
  # set colors
  scale_fill_manual(values = c("Birds" = "#4393c3", "Vascular Plants" = "#d6604d")) +
  # facet by time period
  facet_wrap(~ factor(lc_time_period, levels = time_period_order), 
             ncol = 3) +
  # labels and theme
  labs(x = "Species Change (After - Before)",
       y = "Density",
       fill = "Group") +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "gray80"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Combine the two plots
gain_loss_forest_tws <- plot_grid(jitter_species_change_plot,
                                  density_plot,
                                  labels = c("a)", "b)"), ncol = 1)
# Save figure
ggsave(filename = here("figures", "Figure5_birds_plants_species_changes_forest_tws.png"),
       plot = gain_loss_forest_tws, width = 12, height = 15, dpi = 300)

ggsave(filename = here("figures", "Figure5_birds_plants_species_changes_forest_tws.svg"),
       plot = gain_loss_forest_tws, width = 12, height = 15, dpi = 300)

## 4.2. TWS -> Forest ----------------------------------------------------------

# Jitter plot of values
jitter_tws_forest <- ggplot(tws_forest_data, 
                                     aes(x = species_change, y = taxonomic_group)) +
  # add vertical line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  # add jittered points with different shapes for gain and loss
  geom_jitter(aes(color = change_type, shape = change_type, size = forest_to_tws),
              height = 0.3, alpha = 0.7) +
  # set shapes for gain and loss
  scale_shape_manual(values = c("Gain" = 16, "Loss" = 17, "No Change" = 15)) +
  # set colors for gain and loss
  scale_color_manual(values = gain_loss_colors) +
  # set size scale
  scale_size_continuous(name = "TWS → Forest\nTransition Pixels", 
                        range = c(2, 8)) +
  # facet by time period
  facet_wrap(~ factor(lc_time_period, levels = time_period_order), 
             ncol = 3) +
  # labels and theme
  labs(y = "Group",
       color = "Change Type",
       shape = "Change Type") +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "gray80"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        legend.box = "vertical",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(angle=90))

# Density plot for distribution of species change
density_plot_tws_forest <- ggplot(tws_forest_data, 
                       aes(x = species_change, fill = taxonomic_group)) +
  # add vertical line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  # add density curves
  geom_density(alpha = 0.7) +
  # set colors
  scale_fill_manual(values = c("Birds" = "#4393c3", "Vascular Plants" = "#d6604d")) +
  # facet by time period
  facet_wrap(~ factor(lc_time_period, levels = time_period_order), 
             ncol = 3) +
  # labels and theme
  labs(x = "Species Change (After - Before)",
       y = "Density",
       fill = "Group") +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "gray80"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Combine the two plots
gain_loss_tws_forest <- plot_grid(jitter_tws_forest,
                                  density_plot_tws_forest,
                                  labels = c("a)", "b)"), ncol = 1)
# Save figure
ggsave(filename = here("figures", "Figure5_birds_plants_species_changes_tws_forest.png"),
       plot = gain_loss_tws_forest, width = 12, height = 15, dpi = 300)

ggsave(filename = here("figures", "Figure5_birds_plants_species_changes_tws_forest.svg"),
       plot = gain_loss_tws_forest, width = 12, height = 15, dpi = 300)

## 4.3. All -> Urban -----------------------------------------------------------

# Jitter plot of values
jitter_all_urban <- ggplot(all_urban_data, 
                            aes(x = species_change, y = taxonomic_group)) +
  # add vertical line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  # add jittered points with different shapes for gain and loss
  geom_jitter(aes(color = change_type, shape = change_type, size = forest_to_tws),
              height = 0.3, alpha = 0.7) +
  # set shapes for gain and loss
  scale_shape_manual(values = c("Gain" = 16, "Loss" = 17, "No Change" = 15)) +
  # set colors for gain and loss
  scale_color_manual(values = gain_loss_colors) +
  # set size scale
  scale_size_continuous(name = "All → Artifical Surfaces\nTransition Pixels", 
                        range = c(2, 8)) +
  # facet by time period
  facet_wrap(~ factor(lc_time_period, levels = time_period_order), 
             ncol = 3) +
  # labels and theme
  labs(y = "Group",
       color = "Change Type",
       shape = "Change Type") +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "gray80"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        legend.box = "vertical",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(angle=90))

# Density plot for distribution of species change
density_plot_all_urban <- ggplot(all_urban_data,
                                 aes(x = species_change, fill = taxonomic_group)) +
  # add vertical line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  # add density curves
  geom_density(alpha = 0.7) +
  # set colors
  scale_fill_manual(values = c("Birds" = "#4393c3", "Vascular Plants" = "#d6604d")) +
  # facet by time period
  facet_wrap(~ factor(lc_time_period, levels = time_period_order), 
             ncol = 3) +
  # labels and theme
  labs(x = "Species Change (After - Before)",
       y = "Density",
       fill = "Group") +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "gray80"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Combine the two plots
gain_loss_all_urban <- plot_grid(jitter_all_urban,
                                 density_plot_all_urban,
                                  labels = c("a)", "b)"), ncol = 1)
# Save figure
ggsave(filename = here("figures", "Figure5_birds_plants_species_changes_all_urban.png"),
       plot = gain_loss_all_urban, width = 12, height = 15, dpi = 300)

ggsave(filename = here("figures", "Figure5_birds_plants_species_changes_all_urban.svg"),
       plot = gain_loss_all_urban, width = 12, height = 15, dpi = 300)

# 5. VISUALISE CHANGE & RECORDER EFFORT ----------------------------------------

# Forest -> TWS
rec_effort_forest_tws <- ggplot(forest_tws_data, 
                                aes(x = recorder_effort, y = species_change)) +
  geom_point(alpha = 0.6) +
  scale_x_continuous(labels = comma) +
  geom_smooth() +
  labs(x = "Recorder Effort", 
       y = "Species Change") +
  theme_classic()

# TWS -> Forest
rec_effort_tws_forest <- ggplot(tws_forest_data, 
                                aes(x = recorder_effort, y = species_change)) +
  geom_point(alpha = 0.6) +
  scale_x_continuous(labels = comma) +
  geom_smooth() +
  labs(x = "Recorder Effort", 
       y = "Species Change") +
  theme_classic()

# All -> Urban
rec_effort_all_urban <- ggplot(all_urban_data, 
                               aes(x = recorder_effort, y = species_change)) +
  geom_point(alpha = 0.6) +
  scale_x_continuous(labels = comma) +
  geom_smooth() +
  labs(x = "Recorder Effort", 
       y = "Species Change") +
  theme_classic()

# Combine into single figure
rec_effort <- plot_grid(rec_effort_forest_tws,
                        rec_effort_tws_forest,
                        rec_effort_all_urban,
                        labels = c("a)", "b)", "c)"), ncol = 1)

# END OF SCRIPT ----------------------------------------------------------------