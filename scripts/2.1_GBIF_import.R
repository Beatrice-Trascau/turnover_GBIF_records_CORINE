##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS
# 2.1_GBIF_import
# This script contains code which imports the downloaded GBIF occurrence records
##----------------------------------------------------------------------------##

# 1. IMPORT GBIF DOWNLOAD ------------------------------------------------------

# Download key: "0009294-250802193616735"
# Download link: https://api.gbif.org/v1/occurrence/download/request/0009294-250802193616735.zip
# Download date: 04.08.2025
# DOI: 10.15468/dl.tcgps5

# Import
occurrence <- occ_download_get("0009294-250802193616735") |>
  occ_download_import()

# 2. SAVE TO FILE --------------------------------------------------------------

# Save file
save(occurrence, file = here::here("data","raw_data","occurrence_redownload_4August2025.txt"))
