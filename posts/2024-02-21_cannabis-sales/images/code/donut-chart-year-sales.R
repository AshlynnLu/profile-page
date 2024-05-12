# Load necessary libraries
library(ggplot2)
library(grid)
library(ggnewscale)
library(ggtext)
library(tidyverse)
library(shadowtext)
library(patchwork)
library(dplyr)
library(gridExtra)
library(ggthemes)

# Load colors
source("plots/code/color-setting.R")
# Load function to add watermark
source("src/helper-functions/add-watermark.R")

# Load data
load("data/derived/combined_county.RData")
load("data/derived/total_sales_by_year.RData")

county_sales_donut_chart <- ggplot() +
  geom_bar(data = combined_county,
           aes(x = 0, y = Prop, fill = County),
           stat = "identity",
           width = 3) +
  coord_polar("y", start = 0) +
  xlim(c(-4.5, 1.5)) +
  # Remove background for each subplot
  theme_void() +
  # Facet by year
  facet_wrap(~Year, nrow = 2) +
  # Manually adjust the color of legend
  scale_fill_manual(values = c(BLUE, CYAN, YELLOW, GREEN, DARK_RED, LIGHT_GREEN, OTHER1)) +
  # Text inside donut hole
  geom_text(data = total_sales_by_year,
            aes(label = paste("Total\n", "$", round(Total_Sales/1e6), "m", sep = "")),
            x = -4.5, y = 0, size = 7/.pt, family = "Econ Sans Cnd") +
  # Title, subtitle, and source
  labs(title = "The total cannabis sales by year",
       subtitle = "Colorado, 2014 to 2024",
       caption = "Source: the Colorado Department of Revenue") +
  theme(
    # Plot title
    plot.title = element_text(size = 9.5, face = "bold"),
    # Subtitle
    plot.subtitle = element_text(size = 8, margin = margin(3, 0, 0, 0, unit = "pt")),
    # Source
    plot.caption = element_text(size = 6.5, hjust = 0,
                                family = "Econ Sans Cnd Light",
                                color = source_color,
                                margin = margin(3.5, 0, 0, 0, unit = "pt")),
    # Set background color to white
    plot.background = element_rect(fill = "white", color = "white"),
    # Adjust legend key size
    legend.key.height = unit(7.5, "pt"),
    legend.key.width = unit(5, "pt"),
    # Adjust legend title
    legend.title = element_text(size = 7.5, family = "Econ Sans Cnd Medium"),
    # Adjust legend text size
    legend.text = element_text(size = 7.5, family = "Econ Sans Cnd Light"),
    # Adjust the size of facet titles (year titles)
    strip.text = element_text(size = 8, family = "Econ Sans Cnd",
                              margin = margin(7.5, 0, -15, 0, unit = "pt"), hjust = 0),
    strip.clip = "off",
    # Make the facet closer
    panel.spacing = unit(-1, "pt"),
    # Adjust plot margin
    plot.margin = margin(7.5, 0, 5, 0, unit = "pt")
  )

# Open file to store the plot
png("plots/figure/county_sales_donut_chart.png",
    width = 595*3, height = 290*3, units = "px", res = 300)

county_sales_donut_chart

# Add red box
red_box_on_top()

# Close connection
dev.off()

