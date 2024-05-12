# Import the helper function
source("src/helper-functions/CMYK-to-hexadecimal.R")

# Colors used in plotting
BLUE <- "#006BA2"
CYAN <- "#3EBCD2"
YELLOW <- "#EBB434"
GREEN <- "#379A8B"
DARK_RED <- CMYK_to_hexadecimal(0, 75, 35, 45)
LIGHT_GREEN <- "#86E5D4"
OTHER1 <- "#E9EDF0"
OTHER2 <- "#B7C6CF"

background_color <- CMYK_to_hexadecimal(7.5, 0, 0, 5)
grid_color <- CMYK_to_hexadecimal(10, 0, 0, 25)
source_color <- CMYK_to_hexadecimal(0, 0, 0, 75)

# color scales
BLUES <- c("#00588D", "#1270A8", "#3D89C3", "#5DA4DF", "#7BBFFC", "#98DAFF")

