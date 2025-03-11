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
taxonomic_composition_periods <- bind_rows(period1_occurrences,
                                           period2_occurrences,
                                           period3_occurrences)

## 2.3. Calculate counts and proportions by period and taxonomic group ---------
taxonomic_summary <- taxonomic_composition_periods |>
  filter(!is.na(period)) |>
  group_by(period, taxonomic_group) |>
  summarise(count = n(), .groups = "drop") |>
  group_by(period) |>
  mutate(total = sum(count), proportion = count / total) |>
  ungroup()

# 3. CREATE STACKED PLOT -------------------------------------------------------

# Create a colour pallete
taxonomic_colours <- c("Plants" = "#66C2A5",
                     "Birds" = "#FC8D62",
                     "Arthropods" = "#8DA0CB",
                     "Mammals" = "#E78AC3",
                     "Other" = "#A6D854")

# Define order for the periods
period_order <- c("1997-2000", "2006-2009", "2003-2006", 
                  "2012-2015", "2008-2012", "2015-2018")

# Convert period column to a factor and follow the order
taxonomic_summary <- taxonomic_summary |>
  mutate(period = factor(period, levels = period_order))

# Create stacked barplot
taxonomic_plot <- ggplot(taxonomic_summary,
                         aes(x = period, y = proportion, fill = taxonomic_group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = taxonomic_colours) +
  labs(x = "Time Period", y = "Proportion", fill = "Group") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# Save figure as .png
ggsave(filename = here("figures", "Figure2_taxonomic_composition_across_periods.png"),
       width = 10, height = 8, dpi = 300)

# Save figure as .svg
ggsave(filename = here("figures", "Figure2_taxonomic_composition_across_periods.svg"),
       width = 10, height = 8, dpi = 300)

# 4. VALIDATION CHECKS ---------------------------------------------------------

## 4.1. Check that period assignment was correct -------------------------------

# Count records by year and assigned period to verify correct assignment
year_period_check <- taxonomic_composition_periods |>
  group_by(year, period) |>
  summarise(count = n(), .groupd = "drop") |>
  arrange(year)

# View checks
print(year_period_check, n= 25) #looks correct

## 4.2. Check taxonomic group assignment ---------------------------------------

# Count records by kingdom/class/phylum and assigned taxonomic group
taxonomic_group_check <- taxonomic_composition_periods |>
  group_by(kingdom, class, phylum, taxonomic_group) |>
  summarise(count = n(), .groups = "drop") |>
  arrange(taxonomic_group)

# View checks
print(taxonomic_group_check, n = 130) #looks ok

## 4.3. Check total for each period --------------------------------------------

# Calculate total records per period
period_totals <- taxonomic_composition_periods |>
  group_by(period) |>
  summarise(total_records = n(), .groups = "drop")

# View checks
print(period_totals)

## 4.4. Check proportions add up to 1 for each period --------------------------

# Verify that proportions sum to 1 (or close)
proportion_check <- taxonomic_summary |>
  group_by(period) |>
  summarise(sum_proportion = sum(proportion), .group = "drop")

# View check
print(proportion_check) # all add up to 1

## 4.5. Cross-check with other methods -----------------------------------------

# Get only birds for 2012-2015
bird_2012_2015 <- taxonomic_composition_periods |>
  filter(taxonomic_group == "Birds" & period == "2012-2015")

# Check how many occurrences for birds in 2012-2015
nrow(bird_2012_2015) #4 654 633

# Check how many bird records from 2012-2015 there are in the summary df
summary_bird_2012_2015 <- taxonomic_summary |>
  filter(taxonomic_group == "Birds" & period == "2012-2015")

# Print total value from the summary df
print(summary_bird_2012_2015$count) # 4 654 633

# They are matching! Can try this with more combinations of years and groups.

# END OF SCRIPT ----------------------------------------------------------------