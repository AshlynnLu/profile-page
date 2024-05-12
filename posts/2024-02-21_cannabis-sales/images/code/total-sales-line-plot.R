# Load necessary libraries
library(ggplot2)
library(ggtext)
library(tidyverse)
library(dplyr)
library(grid)
library(gridExtra)

# Load colors
source("plots/code/color-setting.R")
# Load function to add watermark
source("src/helper-functions/add-watermark.R")

# Load data
load("data/derived/data_total_sale.RData")

# Line plot
total_sales_line_plot <- data_total_sale %>%
  ggplot() +
  geom_line(aes(x = Year_Month, y = Total_Sale/1e6, color = "Total Sale"), show.legend = TRUE) +
  geom_line(aes(x = Year_Month, y = Medical/1e6, color = "Medical"), show.legend = TRUE) +
  geom_line(aes(x = Year_Month, y = Retail/1e6, color = "Retail"), show.legend = TRUE) +
  scale_color_manual(values = c("Total Sale" = BLUE, "Medical" = CYAN, "Retail" = YELLOW)) +
  # Make the y grid values to the right
  scale_y_continuous(position = "right") +
  # Title, subtitle, source, and legend title
  labs(title = "Cannabis sales",
       subtitle = "Colorado, 2014 to 2024, $m",
       caption = "Source: the Colorado Department of Revenue",
       color = "Sales type") +
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
    # Remove the title for both axes
    axis.title = element_blank(),
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Remove all grid lines
    panel.grid = element_blank(),
    # Add grid lines for the vertical axis, customizing color and size
    panel.grid.major.y = element_line(color = grid_color, linewidth = unit(0.5, "pt")),
    # Remove tick marks on the vertical axis
    axis.ticks.length.y = unit(0, "pt"),
    # Keep tick marks on horizontal axis
    axis.ticks.length.x = unit(3, "pt"),
    # Make vertical grid values on top of the grid line, adjust font size
    axis.text.y = element_text(size = 7, vjust = 0, family = "Econ Sans Cnd Light"),
    # Font size of the axis label
    axis.text.x = element_text(size = 7, family = "Econ Sans Cnd Light"),
    # Only the bottom line of the vertical axis is painted in black
    axis.line.x.bottom = element_line(color = "black", size = 0.4),
    # Adjust legend position
    legend.position = "top",
    legend.margin = margin(7.5, 0,  0, 0, unit = "pt"),
    # Adjust background underneath legend keys
    legend.key = element_rect(fill = "white"),
    # Adjust legend text size
    legend.text = element_text(size = 7.5, family = "Econ Sans Cnd Light"),
    # Adjust legend title
    legend.title = element_text(size = 7.5, family = "Econ Sans Cnd Medium"),
    # Align legend to left
    legend.justification = c(0, 2),
    # Adjust plot margin
    plot.margin = margin(7.5, 0, 5, 0, unit = "pt")
    )

# Open file to store the plot
png("plots/figure/total-sales-line-plot.png",
    width = 595*3, height = 290*3, units = "px", res = 300)

total_sales_line_plot

# Add red box
red_box_on_top()

# Close connection
dev.off()
