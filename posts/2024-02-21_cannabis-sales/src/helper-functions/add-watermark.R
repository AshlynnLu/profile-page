#' Add a red line and red rectangle on top of the plot
#' like a watermark of visualizations made by The Economist
#'
#' @return Add the watermark to a png
#' @export
#'
#' @examples
#' png("plot.png") # Open file to store the plot
#' plt # Print the plot
#' red_box_on_top() # Apply function
#' dev.off() # Close connection
red_box_on_top <- function(){
  # Add horizontal line on top
  # It goes from x = 0 (left) to x = 1 (right) on the very top of the chart (y = 1)
  # Indicate the line color and width
  grid.lines(
    x = c(0, 1),
    y = 1,
    gp = gpar(col = "#e5001c", lwd = 4)
  )

  # Add rectangle on top-left
  # lwd = 0 means the rectangle does not have an outer line
  # 'just' gives the horizontal and vertical justification
  grid.rect(
    x = 0,
    y = 1,
    width = unit(15, "pt"),
    height = unit(5, "pt"),
    just = c("left", "top"),
    gp = gpar(fill = "#e5001c", lwd = 0)
  )
}

