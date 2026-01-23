##----------------------------------------------------------------------------##
# PAPER 2: CORINE LAND COVER CHANGES AND TURNOVER OF GBIF BIODIVERSITY RECORDS 
# 8_packages
# This script contains code which extracts the version and citation for each
# package used in this project
##----------------------------------------------------------------------------##

# 1. EXPORT REFERENCE FOR ALL PACKAGES -----------------------------------------

# Load setup file
library(renv)

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
                 "climenv", "googledrive")

# Get attached packages
pkgs <- unique(package_vec)

# Exclude base packages
base_pkgs <- c("base", "stats", "utils", "methods",
                   "graphics", "grDevices", "datasets")

# Remove the base packages from the list of used packages
pkgs <- setdiff(pkgs, base_pkgs)

# Freeze/Snapshot the package versions used
renv::snapshot(prompt = FALSE)

# Create BibTex for citation
bib_entries <- lapply(pkgs, function(pkg) {
  cit <- citation(pkg)
  toBibtex(cit)
})

bib_entries  <- unlist(bib_entries)

# Write to BibTex file
writeLines(bib_entries, "references/r_packages.bib")

# 2. EXPORT REFERENCE FOR SPECIFIC PACKAGES ------------------------------------

## 2.1. Terra ------------------------------------------------------------------

# Get citation info
cit <- citation("terra")

# Extract authors
authors <- paste(
  sapply(cit$author, function(a) {
    paste(a$family, paste(a$given, collapse = " "), sep = ", ")
  }),
  collapse = "\nAU  - "
)

# Extract year
year <- cit$year

# Extract title
title <- cit$title

# Extract version
version <- as.character(packageVersion("terra"))

# Extract URL
url <- cit$url

# Build RIS entry
ris <- c("TY  - COMP",
         paste0("TI  - ", title),
         paste0("AU  - ", authors),
         paste0("PY  - ", year),
         paste0("N1  - R package version ", version),
         paste0("UR  - ", url),
         "ER  -")

# Write .RIS reference to file
writeLines(ris, "references/terra.ris")

## 2.2. Sf ------------------------------------------------------------------

# Get citation info
cit <- citation("sf")

# Extract authors
authors <- paste(
  sapply(cit$author, function(a) {
    paste(a$family, paste(a$given, collapse = " "), sep = ", ")
  }),
  collapse = "\nAU  - "
)

# Extract year
year <- cit$year

# Extract title
title <- cit$title

# Extract version
version <- as.character(packageVersion("terra"))

# Extract URL
url <- cit$url

# Build RIS entry
ris <- c("TY  - COMP",
         paste0("TI  - ", title),
         paste0("AU  - ", authors),
         paste0("PY  - ", year),
         paste0("N1  - R package version ", version),
         paste0("UR  - ", url),
         "ER  -")

# Write .RIS reference to file
writeLines(ris, "references/sf.ris")

## 2.3. CoordinateCleaner ------------------------------------------------------

# Get citation info
cit <- citation("CoordinateCleaner")

# Extract authors
authors <- paste(
  sapply(cit$author, function(a) {
    paste(a$family, paste(a$given, collapse = " "), sep = ", ")
  }),
  collapse = "\nAU  - "
)

# Extract year
year <- cit$year

# Extract title
title <- cit$title

# Extract version
version <- as.character(packageVersion("terra"))

# Extract URL
url <- cit$url

# Build RIS entry
ris <- c("TY  - COMP",
         paste0("TI  - ", title),
         paste0("AU  - ", authors),
         paste0("PY  - ", year),
         paste0("N1  - R package version ", version),
         paste0("UR  - ", url),
         "ER  -")

# Write .RIS reference to file
writeLines(ris, "references/CoordinateCleaner.ris")

# END OF SCRIPT ----------------------------------------------------------------