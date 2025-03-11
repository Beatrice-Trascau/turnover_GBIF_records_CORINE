##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 2.3_GBIF_taxonomic_composition
# This script contains code which explores taxonomic composition of occurrence
# records through time
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# Load occurrences
clean_occurrences <- read.csv(here("data", "derived_data", 
                                   "clean_occurrences_15km.txt"))

# 2. PREPARE DATA  -------------------------------------------------------------

## 2.1. Create the groups ------------------------------------------------------

# Create new column for taxonomic group
taxonomic_composition <- clean_occurrences |>
  mutate(taxonomic_group = case_when(kingdom == "Plantae" ~ "Plants",
                                     class == "Aves" ~ "Birds",
                                     phylum == "Arthropoda" ~ "Arthropods",
                                     class == "Mammalia" ~ "Mammals",
                                     TRUE ~ "Other"))
## 2.2. Add before and after periods -------------------------------------------

# Filtering and classification is done step by step to avoid overwriting of
# years that are shared between periods

# First period (2000-2006 LC change, 1997-2000 = before, 2006-2012 = after)
period1_occurrences <- taxonomic_composition |>
  filter(year %in% c(1997:2000, 2006:2009)) |>
  mutate(period = case_when(year %in% 1997:2000 ~ "1997-2000",
                            year %in% 2006:2009 ~ "2006-2009"))

# Second period (2006-2012 LC change, 2003-2006 = before, 2012-2015 = after)
period2_occurrences <- taxonomic_composition |>
  filter(year %in% c(2003:2006, 2012:2015)) |>
  mutate(period = case_when(year %in% 2003:2006 ~ "2003-2006",
                            year %in% 2012:2015 ~ "2012-2015"))

# Third period (2012-2018 LC change, 2008-2012 = before, 2015-2018 = after)
period3_occurrences <- taxonomic_composition |>
  filter(year %in% c(2008:2012, 2015:2018)) |>
  mutate(period = case_when(year %in% 2008:2012 ~ "2008-2012",
                            year %in% 2015:2018 ~ "2015-2018"))

# Combine all periods
occurrences_taxonomy <- bind_rows(period1_occurrences,
                                  period2_occurrences,
                                  period3_occurrences)