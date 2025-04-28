#' Psuedo Probability Density Function for the Variety Normal Distribution
#'
#' Computes a pseudo probability density function (PDF) for the variety normal distribution.
#' This function supports both the heteroskedastic and homoskedastic cases.
#'
#' @param x Numeric value or vector at which to evaluate the density.
#' @param poly An `mpoly` object representing the polynomial
#' @param sd A positive numeric value representing the standard deviation
#' @param homo Logical, default is `TRUE`. If `TRUE`, computes the homoskedastic
#'  variety normal density. If `FALSE`, computes the heteroskedastic variety
#'  normal density.
#'@param log Logical. If `TRUE`, returns the log of the density.
#' @return A numeric value (or vector) representing the density evaluated at `x`.
#'
#' @examples
#' library(mpoly)
#'
#' # Example 1: Heteroskedastic VN with a simple polynomial
#' poly1 <- mp("x^2 - 1")  # g(x) = x^2 - 1
#' pdvnorm(0.5, poly1, sd = 1)  # Homoskedastic (default)
#' pdvnorm(0.5, poly1, sd = 1, homo = FALSE)  # Heteroskedastic
#'
#' # Example 2: Homoskedastic VN with a cubic polynomial
#' poly2 <- mp("x^3 - x")  # g(x) = x^3 - x
#' pdvnorm(-0.5, poly2, sd = 1)
#' pdvnorm(-0.5, poly2, sd = 1, homo = FALSE)
#'
#' # Example 3: Vectorized evaluation
#' x_vals <- seq(-2, 2, length.out = 5)  # Evaluate at multiple points
#' pdvnorm(x_vals, poly1, sd = 1)
#'
#' # Example 4: Handling a polynomial where gradient norm is small
#' poly3 <- mp("x^4 - 2*x^2")  # g(x) = x^4 - 2x^2
#' pdvnorm(0, poly3, sd = 1)  # Should return an error for homoskedastic case
#'
#'
#'
#' @export
pdvnorm <- function(x, poly, sd, homo = TRUE, log = FALSE) {
  if (!inherits(poly, "mpoly")) {
    stop("Error: 'poly' must be an 'mpoly' object.")
  }

  if (!is.numeric(sd) || length(sd) != 1 || sd <= 0) {
    stop("Error: 'sd' must be a single positive numeric value.")
  }

  x <- as.numeric(x)

  # Evaluate the polynomial g(x | beta) at x
  g_func <- suppressMessages(as.function(poly))
  g_val <- g_func(x)

  if (homo) {
    # Gradient
    grad_g <- gradient(poly)
    grad_g_func <- suppressMessages(as.function(grad_g))
    grad_g_val <- sqrt(sum(grad_g_func(x)^2))

    # Homoskedastic log-density
    g_bar <- g_val / grad_g_val
    log_density <- -0.5 * log(2 * pi) - log(sd) - (g_bar^2) / (2 * sd^2)
  } else {
    # Heteroskedastic log-density
    log_density <- -0.5 * log(2 * pi) - log(sd) - (g_val^2) / (2 * sd^2)
  }

  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}

