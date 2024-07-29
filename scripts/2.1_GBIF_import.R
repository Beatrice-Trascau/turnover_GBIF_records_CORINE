##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 2.1_GBIF_import
# This script contains code which imports the downloaded GBIF occurrence records
##----------------------------------------------------------------------------##

# 1. IMPORT GBIF DOWNLOAD ------------------------------------------------------

# Download key: "0016426-240626123714530"

# Import
occurrence_a <- occ_download_get("0016426-240626123714530") %>%
  occ_download_import()

# 2. SAVE TO FILE --------------------------------------------------------------

# Save file
save(occurrence_a, file = here::here("data","raw_data","occurrence_a.rda"))
