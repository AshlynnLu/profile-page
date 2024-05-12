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


total_sales_bar_plot <- ggplot(combined_county_longer) +
  geom_col(aes(x = Total_Sales/1e9, y = County, fill = Sales_Type), width = 0.8, position = "dodge") +
  scale_fill_manual(values=c(BLUE, CYAN)) +
  # Title, subtitle, source, and legend title
  labs(title = "Total cannabis sales",
       subtitle = "Colorado, 2014 to 2024, $bn",
       caption = "Source: the Colorado Department of Revenue",
       fill = "Sales type") +
  # Adjust x axis on top
  scale_x_continuous(expand = c(0, 0), position = "top") +
  theme(
    # Adjust title and caption position align with the plot
    plot.title.position = "plot",
    plot.caption.position = "plot",
    # Plot title
    plot.title = element_text(size = 9.5, hjust = 0, face = "bold"),
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
    panel.grid.major.x = element_line(color = grid_color, linewidth = unit(0.5, "pt")),
    # Remove tick marks on the horizontal axis
    axis.ticks.length.x = unit(0, "pt"),
    # Keep tick marks on vertical axis
    axis.ticks.length.y = unit(3, "pt"),
    # Make vertical grid values on top of the grid line, adjust font size
    axis.text.x = element_text(size = 7, vjust = 0,
                               family = "Econ Sans Cnd Light", color = "black"),
    # Only the left line of the horizontal axis is painted in black
    axis.line.y.left = element_line(color = "black", size = 0.4),
    axis.text.y = element_text(size = 7, hjust = 0,
                               family = "Econ Sans Cnd", color = "black"),
    # Adjust legend position
    legend.position = "top",
    legend.margin = margin(7.5, 0,  0, -32, unit = "pt"),
    # Adjust background underneath legend keys
    legend.key = element_rect(fill = "white"),
    # Adjust legend key size
    legend.key.height = unit(7.5, "pt"),
    legend.key.width = unit(14, "pt"),
    # Adjust legend text size
    legend.text = element_text(size = 7.5, family = "Econ Sans Cnd Light"),
    # Align legend to left
    legend.justification = c(0, 2),
    # Adjust plot margin
    plot.margin = margin(7.5, 0, 5, 0, unit = "pt")
  )  +
  # Adjust legend position to the right and adjust font size
  guides(fill = guide_legend(title.position = "left",
                              title.theme = element_text(hjust = 0, size = 7.5,
                                                         family = "Econ Sans Cnd Medium")))

# Open file to store the plot
png("plots/figure/total-sales-bar-plot.png",
    width = 290*3, height = 290/2*3*3, units = "px", res = 300)
total_sales_bar_plot

# Add red box
red_box_on_top()

# Close connection
dev.off()

