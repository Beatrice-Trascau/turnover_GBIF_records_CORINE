##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 2.2_GBIF_data_cleaning
# This script contains code which cleans the downloaded GBIF occurrence records
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# Load data
load(here("data", "raw_data", "occurrence_a.rda"))

# 2. CLEAN RECORDS -------------------------------------------------------------

# 2.1. Remove material citations, fossil specimens, and living specimens -------

clean_occurrences1 <- occurrence_a |>
  filter(!basisOfRecord %in% c("MATERIAL_CITATION", "FOSSIL_SPECIMEN",
                               "LIVING_SPECIMEN"))

# 2.2. Remove records that are not Animalia, Plantae or Fungi ------------------

clean_occurrences2 <- clean_occurrences1 |>
  filter(kingdom %in% c("Animalia", "Animalia", "Fungi"))

# 2.3. Remove records with no registered species-level information -------------
clean_occurrences3 <- clean_occurrences2 |>
  filter(specificEpithet != "")

# 2.4. Remove duplicate records ------------------------------------------------
clean_occurrences4 <- clean_occurrences3 |>
  distinct()

# 2.5. Remove flagged records --------------------------------------------------

# Identify flagged records
coordinate_flags <- clean_coordinates(x = clean_occurrences4,
                                      lon = "decimalLongitude", 
                                      lat = "decimalLatitude",
                                      species = "species",
                                      test = c("equal", "gbif", "zeros"))

# Get a summary of the records
summary(coordinate_flags) # no flagged records

# 2.5. Remove records with specific coordinate uncertainty ---------------------

# Calculate diagonal of 100m x 100m grid cell
sqrt(100*100)

# Remove records with coord uncertainty >100m
clean_occurrences5 <- clean_occurrences4 |>
  filter(coordinateUncertaintyInMeters < 100 |
           is.na(coordinateUncertaintyInMeters))

# 3. PREP DF FOR ANALYSIS ------------------------------------------------------

# Remove unnecessary columns and add 1 more column
clean_occurrences <- clean_occurrences5 |>
  select(gbifID, identifiedBy, basisOfRecord, occurrenceStatus,
         eventDate, year, countryCode, stateProvince, county, municipality,
         locality, decimalLatitude, decimalLongitude, 
         coordinateUncertaintyInMeters, kingdom, phylum, class, order, family,
         genus, specificEpithet, speciesKey, species, organismQuantity,
         occurrenceStatus, scientificName) |>
  mutate(period = case_when(year %in% c(2000:2005) ~ "2000.2005",
                            year %in% c(2006:2011) ~ "2006.2011",
                            year %in% c(2012:2018) ~ "2012.2018"))

# Save cleaned occurrences
write.csv(clean_occurrences,
          here("data", "derived_data", "cleaned_occurrences_july24.txt"))

# END OF SCRIPT ----------------------------------------------------------------