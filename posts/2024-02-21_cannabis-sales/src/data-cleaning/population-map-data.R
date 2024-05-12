# Load necessary libraries
library(readxl)
library(dplyr)
library(maps)
library(sf)

# Load population data
pop_data <- read_xlsx("data/raw/colorado-population.xlsx", skip = 2)[2:65,c(2, 4)]

colnames(pop_data) <- c("subregion", "pop")

pop_data$subregion <- tolower(gsub(" COUNTY", "", pop_data$subregion))

# Load map data from package "maps"
usa <- map_data("usa")
state <- map_data("state")

# Choose Colorado
colorado <- subset(state, region == "colorado")
counties <- map_data("county")
colorado_county <- subset(counties, region == "colorado")

# Join the map data and the population data
merged_pop <- inner_join(colorado_county, pop_data, by = "subregion")

# Classify the population into categories for plotting
breaks <- c(0, 200e3, 400e3, 600e3, Inf)
labels <- c("[0, 200]", "(200, 400]", "(400, 600]", "(600, Inf)")

merged_pop$category <- cut(merged_pop$pop,
                           breaks = breaks,
                           labels = labels,
                           include.lowest = TRUE)

# Save
save(colorado, file = "data/derived/colorado.RData")
save(colorado_county, file = "data/derived/colorado-county.RData")
save(merged_pop, file = "data/derived/merged-pop.RData")
