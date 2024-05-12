#' Convert CMYK values to hexadecimal color
#'
#' @param C
#' @param M
#' @param Y
#' @param K
#'
#' @return a string containing hexadecimal color
#' @export
#'
#' @examples
#' CMYK_to_hexadecimal(90, 50, 15, 5)
#' CMYK_to_hexadecimal(67, 0, 18, 0)
#' CMYK_to_hexadecimal(12, 30, 70, 0)
CMYK_to_hexadecimal <- function(C, M, Y, K){
  # -----Input Checks ----------------------------------------------------------
  # Check that each input is a numeric value between 0 and 100
  if (!is.numeric(C) || !is.numeric(M) || !is.numeric(Y) || !is.numeric(K)) {
    stop("All inputs must be numeric values.")
  }
  if (any(C < 0 | C > 100) || any(M < 0 | M > 100) || any(Y < 0 | Y > 100) || any(K < 0 | K > 100)) {
    stop("CMYK values must be between 0 and 100.")
  }

  # ----- Function Body --------------------------------------------------------

  # CMYK values / 100
  C <- C / 100
  M <- M / 100
  Y <- Y / 100
  K <- K / 100

  # Convert CMYK to RGB
  R <- round(255 * (1 - C) * (1 - K))
  G <- round(255 * (1 - M) * (1 - K))
  B <- round(255 * (1 - Y) * (1 - K))

  # Convert RGB to hexadecimal
  hex_color <- rgb(R, G, B, maxColorValue = 255)

  # Output the hexadecimal color
  hex_color
}
