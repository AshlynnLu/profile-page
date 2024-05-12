# Import packages
library(dplyr)

# Load raw data
data_raw <- read.csv("data/raw/marijuana-sales-by-county-2014_to-date_report.csv",
                 skip = 4, header = TRUE)
data_raw <- data_raw[1:5035, ]

# Rename the header
colnames(data_raw) <- c("Month", "Year", "County", "Medical", "Retail")

# Make the Month and Year variable into integer
# Create a new variable which include both Year and Month information
data_raw <- data_raw %>%
  mutate(Month = sprintf("%02d", as.integer(Month)),
         Year = as.integer(Year),
         Year_Month = as.Date(paste(Year, Month, "01", sep = "-")))

# Remove the leading space in front "Total"
data_raw$County[data_raw$County == " Total"] <- "Total"

# Transform the sales into integer
data <- data_raw %>%
  mutate(Medical = ifelse(Medical != "NR" & !is.na(Medical),
                          as.integer(gsub("\\$|,", "", Medical)),
                          Medical))
data <- data_raw %>%
  mutate(Retail = ifelse(Retail != "NR" & !is.na(Retail),
                         as.integer(gsub("\\$|,", "", Retail)),
                         Retail))
# Remove NR and NA
data_non_NR_NA <- data_raw %>%
  filter(Medical != "NR" & !is.na(Medical) & Retail != "NR" & !is.na(Retail)) %>%
  mutate(Medical = as.integer(gsub("\\$|,", "", Medical)),
         Retail = as.integer(gsub("\\$|,", "", Retail)))

save(data_non_NR_NA, file = "data/derived/data-non-NR-NA.RData")

################################################################################
# Total sale for drawing line plot
data_total_sale <- data_non_NR_NA %>%
  filter(County == 'Total') %>%
  mutate(Total_Sale = Medical + Retail) %>%
  select(-County)

save(data_total_sale, file = "data/derived/data-total-sale.RData")

################################################################################
# County data for donut chart
data_county <- data_non_NR_NA %>%
  filter(County != 'Total') %>%
  mutate(Total_Sale = Medical + Retail)


# Compute the proportion of each County sales contribute to the total year sales
data_county_total_year_sale <- data_county %>%
  group_by(Year, County) %>%
  summarise(Total_Year_Sale = sum(Total_Sale),
            Medical = sum(Medical),
            Retail = sum(Retail)) %>%
  ungroup(County) %>%
  mutate(Prop = Total_Year_Sale/sum(Total_Year_Sale))

# Filter the top 4 counties with the largest proportions
top_4_counties <- data_county_total_year_sale %>%
  group_by(Year) %>%
  top_n(4, wt = Prop) %>%
  ungroup()

# Summarize the remaining proportions and combine them into "Other"
other_counties <- data_county_total_year_sale %>%
  group_by(Year) %>%
  anti_join(top_4_counties, by = "County") %>%
  summarise(County = "Other",
            Prop = sum(Prop),
            Total_Year_Sale = sum(Total_Year_Sale),
            Medical = sum(Medical),
            Retail = sum(Retail))

# Combine the top 4 counties and "Other"
combined_county <- bind_rows(top_4_counties, other_counties)

# Normalize the proportions within each year
combined_county <- combined_county %>%
  group_by(Year) %>%
  mutate(Prop = Prop / sum(Prop))

top_counties = c("Denver", "Arapahoe", "Boulder", "Adams", "Jefferson", "Pueblo", "Other")

combined_county <- combined_county %>%
  mutate(County = factor(County, levels = top_counties))

# Save
save(combined_county, file = "data/derived/combined-county.RData")

# Calculate total sales grouped by year
total_sales_by_year <- combined_county %>%
  group_by(Year) %>%
  summarize(Total_Sales = sum(Total_Year_Sale))

# Save
save(total_sales_by_year, file = "data/derived/total-sales-by-year.RData")

################################################################################
# County data for bar chart

# Compute the proportion of each County sales contribute to the total year sales
data_county_total_sale <- data_county %>%
  group_by(County) %>%
  summarise(Total_Sale = sum(Total_Sale),
            Medical = sum(Medical),
            Retail = sum(Retail))

# Filter the top 7 counties with the largest proportions
top_7_counties_total <- data_county_total_sale %>%
  top_n(7, wt = Total_Sale)

# Summarize the remaining proportions and combine them into "Other"
other_counties_total <- data_county_total_sale %>%
  anti_join(top_7_counties_total, by = "County") %>%
  summarise(County = "Other",
            Total_Sale = sum(Total_Sale),
            Medical = sum(Medical),
            Retail = sum(Retail))

# Combine the top 5 counties and "Other"
combined_county_total <- bind_rows(top_7_counties_total, other_counties_total)

# Longer form
combined_county_longer <- combined_county_total %>%
  tidyr::pivot_longer(cols = c('Medical', 'Retail'),
                      names_to = "Sales_Type", values_to = "Sales") %>%
  group_by(County, Sales_Type) %>%
  summarise(Total_Sales = sum(Sales))

bar_plot_order = c("Other", "Jefferson", "Pueblo", "Larimer",
                   "Boulder", "Adams", "Arapahoe", "Denver")

combined_county_longer <- combined_county_longer %>%
  mutate(County = factor(County, levels = bar_plot_order))

# Save
save(combined_county_longer, file = "data/derived/combined-county-longer.RData")

