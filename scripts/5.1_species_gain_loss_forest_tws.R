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

# 2. PREPARE DATA FOR ANALYSIS -------------------------------------------------

# Prepare vascular plants data
plants_forest_tws <- plants_turnover |>
  # determine which rows belong to which time period
  mutate(forest_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest no change`,
                                      lc_time_period == "2006-2012" ~ `2006-2012_Forest no change`,
                                      lc_time_period == "2012-2018" ~ `2012-2018_Forest no change`,
                                      TRUE ~ NA_real_),
    forest_to_tws = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest to TWS`,
                              lc_time_period == "2006-2012" ~ `2006-2012_Forest to TWS`,
                              lc_time_period == "2012-2018" ~ `2012-2018_Forest to TWS`,
                              TRUE ~ NA_real_)) |>
  # calculate species gain and loss
  mutate(species_change = total_spp_after - total_spp_before, # net change in species
    # categorize as gain, loss, or no change
    change_type = case_when(species_change > 0 ~ "Gain",
                            species_change < 0 ~ "Loss",
                            TRUE ~ "No Change"),
    # absolute value of change for plotting
    abs_change = abs(species_change),
    # add taxonomic group identifier
    taxonomic_group = "Vascular Plants") |>
  # filter for cells with forest to TWS transitions
  filter(forest_to_tws > 0) |>
  # keep only necessary columns
  select(cell_ID, lc_time_period, species_change, change_type, abs_change, 
         taxonomic_group, forest_to_tws, total_spp_before, total_spp_after)

# Prepare birds data
birds_forest_tws <- birds_turnover |>
  # determine which rows belong to which time period
  mutate(forest_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest no change`,
                                      lc_time_period == "2006-2012" ~ `2006-2012_Forest no change`,
                                      lc_time_period == "2012-2018" ~ `2012-2018_Forest no change`,
                                      TRUE ~ NA_real_),
    forest_to_tws = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest to TWS`,
                              lc_time_period == "2006-2012" ~ `2006-2012_Forest to TWS`,
                              lc_time_period == "2012-2018" ~ `2012-2018_Forest to TWS`,
                              TRUE ~ NA_real_)) |>
  # calculate species gain and loss
  mutate(species_change = total_spp_after - total_spp_before, # net change in species
    # categorize as gain, loss, or no change
    change_type = case_when(species_change > 0 ~ "Gain",
                            species_change < 0 ~ "Loss",
                            TRUE ~ "No Change"),
    # absolute value of change for plotting
    abs_change = abs(species_change),
    # add taxonomic group identifier
    taxonomic_group = "Birds") |>
  # filter for cells with forest to TWS transitions
  filter(forest_to_tws > 0) |>
  # keep only necessary columns
  select(cell_ID, lc_time_period, species_change, change_type, abs_change, 
         taxonomic_group, forest_to_tws, total_spp_before, total_spp_after)

# Combine dataframes into 1
combined_forest_tws <- bind_rows(plants_forest_tws, birds_forest_tws)

# 3. GET SUMMARY STATISTICS ----------------------------------------------------

# Create summary statistics table
summary_stats <- combined_forest_tws |>
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

# 4. VISUALISE GAIN AND LOSS ---------------------------------------------------

# Define time period order
time_period_order <- c("2000-2006", "2006-2012", "2012-2018")

# Define colors for gain and loss
gain_loss_colors <- c("Gain" = "#1a9850", "Loss" = "#d73027", "No Change" = "#4575b4")

# Jitter plot of values
jitter_species_change_plot <- ggplot(combined_forest_tws, 
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
density_plot <- ggplot(combined_forest_tws, 
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

# END OF SCRIPT ----------------------------------------------------------------