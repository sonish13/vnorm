#' Solve a Zero-Dimensional Variety
#'
#' Estimates the solution of a zero-dimensional variety by sampling with
#' [rvnorm()] and averaging the draws. For well-behaved varieties, this method works reliably. When multiple
#' isolated components exist, the sampler may favor some over others, and the
#' posterior mean may fall between true solutions.
#'
#' @param polylist An `mpolyList` or `mpoly` object containing the polynomials.
#' @param sd Numeric scalar, standard deviation used by [rvnorm()] (default: `0.01`).
#' @param n Integer, number of draws used for averaging (default: `1e5`).
#' @param sig_digit Integer, significant digits to round the solution (default: `3`).
#' @param vars Character vector of variable names to average; defaults to `mpoly::vars(polylist)`.
#' @param inc_warmup Logical; include warmup draws when extracting (default: `FALSE`).
#' @param show_message Logical; show Stan sampling messages (default: `FALSE`).
#' @param stanfit Logical; if `TRUE`, return both the Stan fit and the solution (default: `FALSE`).
#' @param ... Additional arguments passed to [rvnorm()].
#'
#' @return If `stanfit = FALSE`, a named numeric vector of rounded posterior means.
#' If `stanfit = TRUE`, a list with `stanfit` (Stan fit) and `results` (the vector).
#'
#'
#' @examples
#' library(mpoly)
#' f1 <- mp("x^2 - y")
#' f2 <- mp("x^2 + y")
#' polylist <- mpolyList(f1, f2)
#'
#' # Posterior-mean solution only
#' variety_solve(polylist, n = 2e4, sd = 0.01, sig_digit = 3)
#'
#' # Return Stan fit and solution
#' out <- variety_solve(polylist, n = 1e4, stanfit = TRUE)
#' out
#'
#'
#'@export

variety_solve <- function(polylist,
                          sd = 0.01,
                          n = 1e5,
                          sig_digit = 3,
                          vars = mpoly::vars(polylist),
                          inc_warmup = FALSE,
                          show_message = FALSE,
                          stanfit = FALSE,
                          ...) {
  if (!inherits(polylist, "mpoly") && !inherits(polylist, "mpolyList"))
    stop("`polylist` must be an mpoly or mpolyList object.", call. = FALSE)

  n_eqns <- length(polylist)
  n_vars <- length(vars)
  # This checks might have to be extended.
  if(n_eqns != n_vars) stop("The variety is not zero-dimensional.")
  samps <- rvnorm(n = n ,
                  poly = polylist,
                  sd = sd,
                  show_messages = show_message,
                  output = "stanfit",
                  ...)
  df <- as.data.frame(samps$draws(format = "df", inc_warmup = inc_warmup))
  df <- df[(nrow(df) - n + 1):nrow(df), mpoly::vars(polylist)]
  row.names(df) <- NULL
  df <- df |> colMeans() |> round(sig_digit)
  if (stanfit) return(list(stanfit = samps, results = df)) else return(df)
}

