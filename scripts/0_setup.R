##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS 
# 0_setup
# This script contains code which loads/installs necessary packages and defines
# functions used in the analysis
##----------------------------------------------------------------------------##

# 1. LOAD/INSTALL PACKAGES NEEDED FOR ANALYIS ----------------------------------

# Function to check to install/load packages

# Define function
install_load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

# Define list of packages
package_vec <- c("here", "terra", "sf", "geodata", "mapview",
                 "tidyverse", "dplyr", "ggplot2", "gt", "cowplot", 
                 "data.table","tidyterra", "patchwork", "styler", 
                 "scales","plotly", "lme4", "DHARMa", "glmmTMB", 
                 "mgcv", "ggspatial", "htmlwidgets","htmltools",  
                 "webshot2", "rgbif", "CoordinateCleaner", "codyn",
                 "gratia", "lattice", "car", "kableExtra",
                 "betareg", "spdep", "corrplot", "leaflet",
                 "viridis", "DT", "broom", "nlme", "ordbetareg",
                 "climenv", "googledrive", "betapart", "bbmle",
                 "patchwork", "ggpubr", "flextable")

# Execute the function
sapply(package_vec, install_load_package)

# 2. CREATE NECESSARY FILE STRUCTURE -------------------------------------------

# Function to create the file structure needed to run the analysis smoothly
create_project_structure <- function(base_path = "turnover_GBIF_records_CORINE") {
  # define the directory structure
  dirs <- c(
    file.path(base_path),
    file.path(base_path, "data"),
    file.path(base_path, "data", "raw_data"),
    file.path(base_path, "data", "raw_data", "raw_norway_shapefile"),
    file.path(base_path, "data", "raw_data", "worldclim"),
    file.path(base_path, "data", "derived_data"),
    file.path(base_path, "data", "derived_data", "reprojected_norway_shapefile"),
    file.path(base_path, "data", "derived_data", "worldclim"),
    file.path(base_path, "data", "models"),
    file.path(base_path, "data", "models", "exploratory"),
    file.path(base_path, "data", "models", "final"),
    file.path(base_path, "scripts"),
    file.path(base_path, "figures")
  )
  
  # create directories if they don't exist
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    } else {
      cat("Directory already exists:", dir, "\n")
    }
  }
  
  cat("\nProject structure setup complete!\n")
}

# Run function
create_project_structure()

# 3. DOWNLOAD FILES ------------------------------------------------------------

# Function to check if files are in directory and then download them if they
# aren't

download_files <- function(urls, filenames, dir = here("data", "raw_data")) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  for (i in seq_along(urls)) {
    file_path <- file.path(dir, filenames[i])
    if (!file.exists(file_path)) {
      download.file(urls[i], file_path)
    }
  }
}

# 4. READ RASTERS --------------------------------------------------------------

# Function to read rasters (from the raw_data file)

read_rasters <- function(filenames, dir = here("data/raw_data")) {
  rasters <- lapply(filenames, function(x) {
    file_path <- file.path(dir, x)
    if (!file.exists(file_path)) {
      stop(paste("File does not exist:", file_path))
    }
    rast(file_path)
  })
  return(do.call(c, rasters))
}

# 5. CROP AND MASK RASTERS TO NORWAY -------------------------------------------

# Function to crop and mask rasters to Norway

crop_mask_to_norway <- function(raster_stack, norway_shape) {
  return(crop(raster_stack, norway_shape, mask = TRUE))
}

# 6. MODIFY CLASSES IN RASTERS -------------------------------------------------

# Function to modify class values in the rasters

modify_class_values <- function(raster_stack, class_modifications) {
  modified_stack <- raster_stack
  for (mod in class_modifications) {
    modified_stack <- app(modified_stack, fun = function(x) {
      x[x %in% mod$from] <- mod$to
      return(x)
    })
  }
  return(modified_stack)
}

# 7. CALCULATE FOREST -> TWS TRANSITIONS ---------------------------------------

# Function to calculate Forest -> TWS transitions between two time periods
analyse_forest_transition <- function(rast_t1, rast_t2) {
  # Create transition raster
  # 0 = non-forest in t1
  # 1 = forest remained forest
  # 2 = forest converted to TWS
  # 3 = other forest conversion
  
  # Create "dummy" raster with matching spatial characteristics
  transition <- rast_t1
  
  # Give it 0 values to show non-forested areas
  transition[] <- 0
  
  # Identify forest cells in initial layer
  forest_t1 <- rast_t1 == 250
  
  # For forest cells in initial layer, categorize changes:
  transition[forest_t1] <- case_when(
    # Forest remained forest
    rast_t2[forest_t1] == 250 ~ 1,
    # Forest converted to shrubland
    rast_t2[forest_t1] == 590 ~ 2,
    # Forest converted to something else
    TRUE ~ 3
  )
  
  return(transition)
}


# 8. CALCULATE TWS -> FORESTS TRANSITIONS --------------------------------------

# Function to calculate TWS -> Forest transitions between two time periods
analyse_tws_transition <- function(rast_t1, rast_t2) {
  # Create transition raster
  # 0 = non-TWS in t1
  # 1 = TWS remained TWS
  # 2 = TWS converted to forest
  # 3 = other TWS conversion
  
  # Create "dummy" raster with matching spatial characteristics
  transition <- rast_t1
  
  # Give it 0 values to show non-forested areas
  transition[] <- 0
  
  # Identify TWS cells in initial layer
  tws_t1 <- rast_t1 == 590
  
  # For TWS cells in initial layer, categorize changes:
  transition[tws_t1] <- case_when(
    # TWS remained TWS
    rast_t2[tws_t1] == 590 ~ 1,
    # TWS converted to Forest
    rast_t2[tws_t1] == 250 ~ 2,
    # TWS converted to something else
    TRUE ~ 3
  )
  
  return(transition)
}

# 9. CALCULATE ALLL -> URBAN TRANSITIONS ---------------------------------------

# Function to calculate All -> Urban transitions between two time periods
analyse_urban_conversion <- function(rast_t1, rast_t2) {
  # Create transition raster
  # 0 = already urban in t1
  # 1 = forest to urban
  # 2 = shrubland to urban
  # 3 = complex agriculture to urban
  # 4 = agriculture & vegetation to urban
  # 5 = moors, heathland & grassland to urban
  # 6 = sparse vegetation to urban
  # 7 = no conversion to urban
  
  # Create "dummy" raster with matching spatial characteristics
  transition <- rast_t1
  
  # Give it values = 7 = no conversion to urban
  transition[] <- 7
  
  # Mark existing urban areas as 0
  urban_t1 <- rast_t1 == 1
  transition[urban_t1] <- 0
  
  # Identify cells that became urban in t2
  urban_t2 <- rast_t2 == 1
  
  # For cells that became urban, categorize their original land cover
  transition[urban_t2 & !urban_t1] <- case_when(
    # Forest to urban
    rast_t1[urban_t2 & !urban_t1] == 250 ~ 1,
    # Shrubland to urban
    rast_t1[urban_t2 & !urban_t1] == 590 ~ 2,
    # Complex agriculture to urban
    rast_t1[urban_t2 & !urban_t1] == 80 ~ 3,
    # Agriculture & vegetation to urban
    rast_t1[urban_t2 & !urban_t1] == 103 ~ 4,
    # Moors, heathland & grassland to urban
    rast_t1[urban_t2 & !urban_t1] == 380 ~ 5,
    # Sparse vegetation to urban
    rast_t1[urban_t2 & !urban_t1] == 711 ~ 6,
    # This TRUE case should never be reached as we've covered all classes
    TRUE ~ 7
  )
  
  return(transition)
}

# 10. AGGREGATE RASTERS TO LARGER GRIDS -----------------------------------------

# This function aggregates the transition rasters (created with the functions 
#   above) to larger grid sizes

aggregate_transitions <- function(transition_raster, factor) {
  # Create empty list to store aggregated count rasters
  agg_counts <- list()
  
  # Get the levels to know what we're counting
  cats <- levels(transition_raster)[[1]]
  
  # For each layer in the transition raster
  for(layer in 1:nlyr(transition_raster)) {
    # Create empty list for this time period
    layer_counts <- list()
    
    # For each transition category (except "no change" and "no conversion")
    for(val in cats$value) {
      # Create binary raster for this transition
      binary <- transition_raster[[layer]] == val
      
      # Aggregate and sum the occurrences
      count_rast <- terra::aggregate(binary, fact=factor, fun="sum")
      
      # Name the raster based on the category
      names(count_rast) <- paste0(names(transition_raster)[layer], "_",
                                  cats$class[cats$value == val])
      
      # Add to list
      layer_counts[[length(layer_counts) + 1]] <- count_rast
    }
    
    # Combine all counts for this time period
    agg_counts[[layer]] <- rast(layer_counts)
  }
  
  # Combine all time periods
  return(rast(agg_counts))
}

# 11. CHECK IF CELLS ARE NA ----------------------------------------------------

# Function to check if certain cells in a raster were already NA before the masking
#   this function is used to check that masking of cells with >50% of their area
#   outside of the boundary was done correctly
check_na_values <- function(raster_layer, coords) {
  # Extract values directly
  extracted <- terra::extract(raster_layer, as.matrix(coords[, c("x", "y")]))
  
  # Print the structure to understand the output
  print(str(extracted))
  
  # Check if extraction worked and what format it returned
  if (is.null(extracted)) {
    return(NA)
  } else if (is.data.frame(extracted)) {
    # If it's a data frame, find the value column (usually the second column)
    if (ncol(extracted) >= 2) {
      return(sum(is.na(extracted[,2])))
    } else if (ncol(extracted) == 1) {
      return(sum(is.na(extracted[,1])))
    }
  } else if (is.vector(extracted)) {
    # If it's a vector, count NAs directly
    return(sum(is.na(extracted)))
  }
  
  # Fallback
  return(NA)
}

# 12. CALCULATE JACCARD'S DISSIMILARITY INDEX ----------------------------------

# Function to calculate Jaccard dissimilarity for a period pair
calculate_jaccard_for_periods <- function(df, before_period, after_period, change_period) {
  # Get data for before and after periods
  before_data <- df |>
    filter(period == before_period)
  after_data <- df |>
    filter(period == after_period)
  
  # Combine data
  combined_data <- before_data |>
    full_join(after_data, by = "cell", suffix = c("_before", "_after"))
  
  # Calculate Jaccard dissimilarity and species counts
  jaccard_results <- combined_data |>
    group_by(cell) |>
    summarize(
      species_before = list(unique(species_before)),
      species_after = list(unique(species_after)),
      n_species_before = length(unique(species_before)),
      n_species_after = length(unique(species_after)),
      jaccard_dissimilarity = 1 - length(intersect(species_before[[1]], 
                                                   species_after[[1]])) / 
        length(union(species_before[[1]], 
                     species_after[[1]])),
      .groups = 'drop'
    ) |>
    mutate(
      before_period = before_period,
      after_period = after_period,
      change_period = change_period
    )
  
  return(jaccard_results)
}

# 13. DOWNLOAD FOLDERS FROM GOOGLE DRIVE ---------------------------------------

# Function to download a folder and its contents from Google Drive
download_drive_folder <- function(folder_id, folder_name, local_path) {
  
  # Create local folder
  folder_path <- file.path(local_path, folder_name)
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
  
  # List all files in the Google Drive folder
  files_in_folder <- drive_ls(as_id(folder_id))
  
  cat("  Found", nrow(files_in_folder), "files in", folder_name, "\n")
  
  # Download each file
  for (i in seq_len(nrow(files_in_folder))) {
    file_info <- files_in_folder[i, ]
    local_file_path <- file.path(folder_path, file_info$name)
    
    # Skip if file already exists locally
    if (!file.exists(local_file_path)) {
      cat("    Downloading:", file_info$name, "\n")
      drive_download(as_id(file_info$id), 
                     path = local_file_path,
                     overwrite = FALSE)
    } else {
      cat("    Skipping (already exists):", file_info$name, "\n")
    }
  }
}

# 14. MODEL OUTPUT PLOTTING FUNCTIONS ------------------------------------------

## 14.1. Forest to Transitional Woodland Shrub ---------------------------------

# Function to create coefficient plot for a single model
create_coef_plot_forest <- function(model, title, show_y_axis = TRUE) {
  
  # get model summary
  model_summary <- summary(model)
  
  # create a dataframe of the coefficients
  coef_df <- data.frame(term = names(model_summary$tTable[, "Value"]),
                        estimate = model_summary$tTable[, "Value"],
                        std.error = model_summary$tTable[, "Std.Error"],
                        statistic = model_summary$tTable[, "t-value"],
                        p.value = model_summary$tTable[, "p-value"]) |>
    mutate(conf.low = estimate - 1.96 * std.error,
           conf.high = estimate + 1.96 * std.error,
           significance = case_when(p.value < 0.001 ~ "***",
                                    p.value < 0.01 ~ "**", 
                                    p.value < 0.05 ~ "*",
                                    p.value < 0.1 ~ ".",
                                    TRUE ~ ""),
           effect_type = case_when(term == "(Intercept)" ~ "Intercept",
                                   p.value < 0.05 & estimate < 0 ~ "Negative (sig.)",
                                   p.value < 0.05 & estimate > 0 ~ "Positive (sig.)",
                                   TRUE ~ "Non-significant"),
           term_clean = case_when(term == "(Intercept)" ~ "Intercept",
                                  term == "forest_to_tws" ~ "Forest → TWS",
                                  term == "forest_no_change" ~ "Forest (no change)",
                                  term == "delta_recorder_effort" ~ "ΔRecorder effort",
                                  term == "log_recorder_effort" ~ "log(Recorder effort)",
                                  term == "lc_time_period2006-2012" ~ "Period: 2006-2012",
                                  term == "lc_time_period2012-2018" ~ "Period: 2012-2018",
                                  term == "temp_change" ~ "Temperature change",
                                  term == "precip_change" ~ "Precipitation change",
                                  TRUE ~ term))
  
  # split into intercept and effects
  intercept_df <- coef_df |> filter(term == "(Intercept)")
  effects_df <- coef_df |> 
    filter(term != "(Intercept)") |>
    mutate(term_clean = factor(term_clean, 
                               levels = rev(c("Forest → TWS",
                                              "Forest (no change)", 
                                              "Temperature change",
                                              "Precipitation change",
                                              "ΔRecorder effort",
                                              "log(Recorder effort)",
                                              "Period: 2006-2012",
                                              "Period: 2012-2018"))))
  
  # define colors
  effect_colors <- c( "Negative (sig.)" = "#9C27B0",
                      "Positive (sig.)" = "#FF9800",  
                      "Non-significant" = "grey60",
                      "Intercept" = "grey30")
  
  # create intercept plot (top)
  plot_intercept <- ggplot(intercept_df, aes(x = estimate, y = term_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                   height = 0.2, linewidth = 0.7, color = "grey30") +
    geom_point(aes(color = effect_type), size = 3, shape = 15) +
    geom_text(aes(x = conf.high, label = significance), 
              hjust = -0.3, vjust = 0.5, size = 3.5, fontface = "bold") +
    scale_color_manual(values = effect_colors, guide = "none") +
    labs(x = NULL, y = NULL, title = title) +
    theme_classic() +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 10, color = "black", face = "bold"),
          axis.text.x = element_text(size = 9, color = "black"),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
          plot.margin = margin(t = 5, r = 10, b = 2, l = 5)) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)))
  
  # create effects plot (bottom)
  plot_effects <- ggplot(effects_df, aes(x = estimate, y = term_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                   height = 0.2, linewidth = 0.7, color = "grey30") +
    geom_point(aes(color = effect_type), size = 3) +
    geom_text(aes(x = conf.high, label = significance), 
              hjust = -0.3, vjust = 0.5, size = 3.5, fontface = "bold") +
    scale_color_manual(values = effect_colors,
                       name = "Effect type",
                       breaks = c("Positive (sig.)", "Negative (sig.)", "Non-significant"),
                       labels = c("Positive (p < 0.05)", "Negative (p < 0.05)", "Non-significant")) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 9, color = "black"),
          legend.position = "none",
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
          plot.margin = margin(t = 2, r = 10, b = 5, l = 5)) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)))
  
  # show or hide y-axis labels
  if (!show_y_axis) {
    plot_intercept <- plot_intercept + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_effects <- plot_effects + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  # combine intercept and effects for this model
  combined <- plot_intercept / plot_effects + 
    plot_layout(heights = c(1, 4))
  
  return(combined)
}

## 14.2. Transitional Woodland Shrub to Forest ---------------------------------

# Function to create coefficient plot for a single model
create_coef_plot_tws <- function(model, title, show_y_axis = TRUE) {
  
  # get model summary
  model_summary <- summary(model)
  
  # create a dataframe of the coefficients
  coef_df <- data.frame(term = names(model_summary$tTable[, "Value"]),
                        estimate = model_summary$tTable[, "Value"],
                        std.error = model_summary$tTable[, "Std.Error"],
                        statistic = model_summary$tTable[, "t-value"],
                        p.value = model_summary$tTable[, "p-value"]) |>
    mutate(conf.low = estimate - 1.96 * std.error,
           conf.high = estimate + 1.96 * std.error,
           significance = case_when(p.value < 0.001 ~ "***",
                                    p.value < 0.01 ~ "**", 
                                    p.value < 0.05 ~ "*",
                                    p.value < 0.1 ~ ".",
                                    TRUE ~ ""),
           effect_type = case_when(term == "(Intercept)" ~ "Intercept",
                                   p.value < 0.05 & estimate < 0 ~ "Negative (sig.)",
                                   p.value < 0.05 & estimate > 0 ~ "Positive (sig.)",
                                   TRUE ~ "Non-significant"),
           term_clean = case_when(term == "(Intercept)" ~ "Intercept",
                                  term == "tws_to_forest" ~ "TWS → Forest",
                                  term == "tws_no_change" ~ "TWS (no change)",
                                  term == "delta_recorder_effort" ~ "ΔRecorder effort",
                                  term == "log_recorder_effort" ~ "log(Recorder effort)",
                                  term == "lc_time_period2006-2012" ~ "Period: 2006-2012",
                                  term == "lc_time_period2012-2018" ~ "Period: 2012-2018",
                                  term == "temp_change" ~ "Temperature change",
                                  term == "precip_change" ~ "Precipitation change",
                                  TRUE ~ term))
  
  # split into intercept and effects
  intercept_df <- coef_df |> filter(term == "(Intercept)")
  effects_df <- coef_df |> 
    filter(term != "(Intercept)") |>
    mutate(term_clean = factor(term_clean, 
                               levels = rev(c("TWS → Forest",
                                              "TWS (no change)", 
                                              "Temperature change",
                                              "Precipitation change",
                                              "ΔRecorder effort",
                                              "log(Recorder effort)",
                                              "Period: 2006-2012",
                                              "Period: 2012-2018"))))
  
  # define colors
  effect_colors <- c( "Negative (sig.)" = "#9C27B0",
                      "Positive (sig.)" = "#FF9800",  
                      "Non-significant" = "grey60",
                      "Intercept" = "grey30")
  
  # create intercept plot (top)
  plot_intercept <- ggplot(intercept_df, aes(x = estimate, y = term_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                   height = 0.2, linewidth = 0.7, color = "grey30") +
    geom_point(aes(color = effect_type), size = 3, shape = 15) +
    geom_text(aes(x = conf.high, label = significance), 
              hjust = -0.3, vjust = 0.5, size = 3.5, fontface = "bold") +
    scale_color_manual(values = effect_colors, guide = "none") +
    labs(x = NULL, y = NULL, title = title) +
    theme_classic() +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 10, color = "black", face = "bold"),
          axis.text.x = element_text(size = 9, color = "black"),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
          plot.margin = margin(t = 5, r = 10, b = 2, l = 5)) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)))
  
  # create effects plot (bottom)
  plot_effects <- ggplot(effects_df, aes(x = estimate, y = term_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                   height = 0.2, linewidth = 0.7, color = "grey30") +
    geom_point(aes(color = effect_type), size = 3) +
    geom_text(aes(x = conf.high, label = significance), 
              hjust = -0.3, vjust = 0.5, size = 3.5, fontface = "bold") +
    scale_color_manual(values = effect_colors,
                       name = "Effect type",
                       breaks = c("Positive (sig.)", "Negative (sig.)", "Non-significant"),
                       labels = c("Positive (p < 0.05)", "Negative (p < 0.05)", "Non-significant")) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 9, color = "black"),
          legend.position = "none",
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
          plot.margin = margin(t = 2, r = 10, b = 5, l = 5)) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)))
  
  # show or hide y-axis labels
  if (!show_y_axis) {
    plot_intercept <- plot_intercept + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_effects <- plot_effects + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  # combine intercept and effects for this model
  combined <- plot_intercept / plot_effects + 
    plot_layout(heights = c(1, 4))
  
  return(combined)
}

## 14.3. All other classes to Artificial Surfaces ------------------------------

# Function to create coefficient plot for a single model
create_coef_plot_urban <- function(model, title, show_y_axis = TRUE) {
  
  # get model summary
  model_summary <- summary(model)
  
  # create a dataframe of the coefficients
  coef_df <- data.frame(term = names(model_summary$tTable[, "Value"]),
                        estimate = model_summary$tTable[, "Value"],
                        std.error = model_summary$tTable[, "Std.Error"],
                        statistic = model_summary$tTable[, "t-value"],
                        p.value = model_summary$tTable[, "p-value"]) |>
    mutate(conf.low = estimate - 1.96 * std.error,
           conf.high = estimate + 1.96 * std.error,
           significance = case_when(p.value < 0.001 ~ "***",
                                    p.value < 0.01 ~ "**", 
                                    p.value < 0.05 ~ "*",
                                    p.value < 0.1 ~ ".",
                                    TRUE ~ ""),
           effect_type = case_when(term == "(Intercept)" ~ "Intercept",
                                   p.value < 0.05 & estimate < 0 ~ "Negative (sig.)",
                                   p.value < 0.05 & estimate > 0 ~ "Positive (sig.)",
                                   TRUE ~ "Non-significant"),
           term_clean = case_when(term == "(Intercept)" ~ "Intercept",
                                  term == "all_to_urban" ~ "All → Artificial Surfaces",
                                  term == "urban_no_change" ~ "No Change",
                                  term == "delta_recorder_effort" ~ "ΔRecorder effort",
                                  term == "log_recorder_effort" ~ "log(Recorder effort)",
                                  term == "lc_time_period2006-2012" ~ "Period: 2006-2012",
                                  term == "lc_time_period2012-2018" ~ "Period: 2012-2018",
                                  term == "temp_change" ~ "Temperature change",
                                  term == "precip_change" ~ "Precipitation change",
                                  TRUE ~ term))
  
  # split into intercept and effects
  intercept_df <- coef_df |> filter(term == "(Intercept)")
  effects_df <- coef_df |> 
    filter(term != "(Intercept)") |>
    mutate(term_clean = factor(term_clean, 
                               levels = rev(c("All → Artificial Surfaces",
                                              "No Change", 
                                              "Temperature change",
                                              "Precipitation change",
                                              "ΔRecorder effort",
                                              "log(Recorder effort)",
                                              "Period: 2006-2012",
                                              "Period: 2012-2018"))))
  
  # define colors
  effect_colors <- c( "Negative (sig.)" = "#9C27B0",
                      "Positive (sig.)" = "#FF9800",  
                      "Non-significant" = "grey60",
                      "Intercept" = "grey30")
  
  # create intercept plot (top)
  plot_intercept <- ggplot(intercept_df, aes(x = estimate, y = term_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                   height = 0.2, linewidth = 0.7, color = "grey30") +
    geom_point(aes(color = effect_type), size = 3, shape = 15) +
    geom_text(aes(x = conf.high, label = significance), 
              hjust = -0.3, vjust = 0.5, size = 3.5, fontface = "bold") +
    scale_color_manual(values = effect_colors, guide = "none") +
    labs(x = NULL, y = NULL, title = title) +
    theme_classic() +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 10, color = "black", face = "bold"),
          axis.text.x = element_text(size = 9, color = "black"),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
          plot.margin = margin(t = 5, r = 10, b = 2, l = 5)) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)))
  
  # create effects plot (bottom)
  plot_effects <- ggplot(effects_df, aes(x = estimate, y = term_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                   height = 0.2, linewidth = 0.7, color = "grey30") +
    geom_point(aes(color = effect_type), size = 3) +
    geom_text(aes(x = conf.high, label = significance), 
              hjust = -0.3, vjust = 0.5, size = 3.5, fontface = "bold") +
    scale_color_manual(values = effect_colors,
                       name = "Effect type",
                       breaks = c("Positive (sig.)", "Negative (sig.)", "Non-significant"),
                       labels = c("Positive (p < 0.05)", "Negative (p < 0.05)", "Non-significant")) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 9, color = "black"),
          legend.position = "none",
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
          plot.margin = margin(t = 2, r = 10, b = 5, l = 5)) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)))
  
  # show or hide y-axis labels
  if (!show_y_axis) {
    plot_intercept <- plot_intercept + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_effects <- plot_effects + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  # combine intercept and effects for this model
  combined <- plot_intercept / plot_effects + 
    plot_layout(heights = c(1, 4))
  
  return(combined)
}


# END OF SCRIPT ----------------------------------------------------------------