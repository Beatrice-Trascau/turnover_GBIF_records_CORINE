##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 4.1_model_setup_turnover_forest_tws
# This script contains code which sets up the models exploring the imapct of 
# Forest -> TWS land cover transition on temporal turnover
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# Turnover data
load(here("data", "derived_data", 
          "all_periods_turnover_all_land_cover_chanegs_15km.rda"))

# Forest -> TWS raster
clc_status_15km_forest_tws_masked <- rast(here("data", "derived_data", 
                                               "clc_status_15km_forest_tws_masked.tif"))

# 2. PREPARE DATA FOR ANALYSIS -------------------------------------------------

## 2.1. Select only Forest -> TWS columns --------------------------------------

# Check column names
colnames(all_periods_turnover_all_land_cover_chanegs_15km)

# Also rename the columns for easier manipulation of df
turnover_forest_tws_15km <- all_periods_turnover_all_land_cover_chanegs_15km |>
  select(-c('2000-2006_TWS no change', '2000-2006_TWS to Forest',
            '2000-2006_Urban_no_change', '2000-2006_all_to_urban',
            '2006-2012_TWS no change', '2006-2012_TWS to Forest',
            '2006-2012_Urban_no_change', '2006-2012_all_to_urban',
            '2012-2018_TWS no change', '2012-2018_TWS to Forest',
            '2012-2018_Urban_no_change', '2012-2018_all_to_urban')) |>
  # determine which rows belong to which time period
  mutate(forest_no_change = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest no change`,
                                      lc_time_period == "2006-2012" ~ `2006-2012_Forest no change`,
                                      lc_time_period == "2012-2018" ~ `2012-2018_Forest no change`,
                                      TRUE ~ NA_real_),
         forest_to_tws = case_when(lc_time_period == "2000-2006" ~ `2000-2006_Forest to TWS`,
                                   lc_time_period == "2006-2012" ~ `2006-2012_Forest to TWS`,
                                   lc_time_period == "2012-2018" ~ `2012-2018_Forest to TWS`,
                                   TRUE ~ NA_real_)) |>
  # remove columns no longer required
  select(-`2000-2006_Forest no change`, -`2006-2012_Forest no change`, 
         -`2012-2018_Forest no change`,-`2000-2006_Forest to TWS`, 
         -`2006-2012_Forest to TWS`, -`2012-2018_Forest to TWS`)

# Check transformation went ok
head(turnover_forest_tws_15km)

## 2.2. Transform JDI values for beta regression -------------------------------

# JDI is in [0,1] but beta regression requires that values do not touch 0 and 1
# so we will transform the JDI values according to this formula: 
# ( Y * (N - 1) + 0.5 ) / N, Y = response variable, N = sample size

# Get N
N <- nrow(turnover_forest_tws_15km)

# Calculate new JDI values
turnover_forest_tws_15km <- turnover_forest_tws_15km |>
  mutate(JDI_beta = (JDI * (N - 1) + 0.5) / N)

# Check the rest of the values are what you expect them to be
glimpse(turnover_forest_tws_15km) # everything looks ok

# 3. RUN BETA GLM --------------------------------------------------------------

# Run model
model1.1 <- betareg::betareg(JDI_beta ~ forest_to_tws + forest_no_change + 
                               delta_recorder_effort + recorder_effort + lc_time_period,
                             data = turnover_forest_tws_15km)

# Get summary
summary(model1.1)

# Extarct residuals from model
model1.1_residuals <- residuals(model1.1)

# Add residuals to df
turnover_forest_tws_15km$residuals <- model1.1_residuals

# 4. CHECK SPATIAL AUTOCORRELATION OF RESIDUALS --------------------------------

## 4.1. Prepare data for testing of autocorrelation ----------------------------

# Create reference grid
reference_grid <- clc_status_15km_forest_tws_masked[[1]]

# Convert turnover df to sf object
turnover_sf <- st_as_sf(turnover_forest_tws_15km,
                        coords = c("x", "y"),
                        crs = st_crs(reference_grid))

## 4.2. Create spatial neighbours ----------------------------------------------

# Create a neighbour list using k-nearest neighbours (k = 5)
coords_matrix <- st_coordinates(turnover_sf)
neighbours <- knn2nb(knearneigh(coords_matrix, k = 5))

# Convert to spatial weights matrix
weights <- nb2listw(neighbours, style = "W")

## 4.3. Calculate Moran's I ----------------------------------------------------

# Compute Moran's I test for spatial autocorrelation
moran_test <- moran.test(turnover_sf$residuals, listw = weights)

# Print results
print(moran_test) # strong spatial autocorrelation of residuals 

# Create Moran scatterplot to visualise the relationship
moran_plot <- moran.plot(turnover_sf$residuals, listw = weights,
                         xlab = "Model Residuals", 
                         ylab = "Spatially Lagged Residuals")

# 5. PLOT MORAN'S I RESULTS ----------------------------------------------------

## 5.1. Calculate local Moran's I values ---------------------------------------

# Calculate local Moran's I wiht zero.policy
local_moran <- localmoran(turnover_sf$residuals, weights, zero.policy = TRUE)

# Add local Moran's I statistics to the spatial dataframe
turnover_sf$local_moran_i <- local_moran[, 1] # I statistics
turnover_sf$local_moran_p <- local_moran[, 5] # p-value

## 5.2. Handle missing values --------------------------------------------------

# Check for missing values in Moran's I column
missing_values <- sum(is.na(turnover_sf$local_moran_i))
cat("Number of cells with missing Local Moran's I values:", missing_values, "\n") #0

# Filter out NA for mapping
turnover_sf_clean <- turnover_sf |>
  filter(!is.na(local_moran_i))

## 5.3. Map of local Moran I values --------------------------------------------

# Extract coordinates
coords <- st_coordinates(turnover_sf)
turnover_df <- cbind(st_drop_geometry(turnover_sf), coords)

# Plot map
plot1 <- ggplot(turnover_df |> 
                  filter(!is.na(local_moran_i)), 
                aes(x = X, y = Y, fill = local_moran_i)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "Local Moran's I") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

# Save map
# ggsave(here("figures", "SupplementaryFigure4a_Local_Moran_I_Forest_TWS_15km.png"),
#        plot1, width = 10, height = 8)

## 5.3. Map of significance of local Moran's I values --------------------------

# Plot map
plot2 <- ggplot(turnover_df |> 
                  filter(!is.na(local_moran_i)), 
                aes(x = X, y = Y, fill = local_moran_p < 0.05)) +
  geom_tile() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey"),
                    name = "Significant",
                    labels = c("FALSE" = "Not significant", "TRUE" = "p < 0.05")) +
  theme_classic() +
  theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank())

# Save figure 
# ggsave(here("figures", "SupplementaryFigure4b_Local_Moran_I_Significance_Forest_TWS_15km.png"),
#        plot1, width = 10, height = 8)

## 5.4. Combine into single plot -----------------------------------------------

# Combine plots
forest_tws_15km_local_morans_I <- plot_grid(plot1, plot2,
                                            labels = c("a)", "b)"),
                                            nrow = 1, align = "h")

# Save to file
ggsave(here("figures", "SupplementaryFigure4_Local_Moran_I_Forest_TWS_15km.png"),
       forest_tws_15km_local_morans_I, width = 12, height = 8)

# 6. PRINT MORAN'S I SUMMARY STATISTICS ----------------------------------------

# Calculate % of cells with significant spatial autocorrelation
sig_percentage <- mean(turnover_sf$local_moran_p < 0.05, na.rm = TRUE) * 100

# Print summary statistics
cat("\n=== SPATIAL AUTOCORRELATION SUMMARY ===\n")
cat("Global Moran's I:", round(moran_test$estimate[1], 4), "\n")
cat("p-value:", format.pval(moran_test$p.value), "\n")
cat("Percentage of cells with significant local autocorrelation:", round(sig_percentage, 2), "%\n")

# END OF SCRIPT ----------------------------------------------------------------