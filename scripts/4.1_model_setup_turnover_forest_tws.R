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


## 4.4. Map of Moran's I values ------------------------------------------------