#' Rejection Sampler for the Variety Normal Distribution
#'
#' Perform rejection sampling to generate draws from a
#' variety normal distribution.
#'
#' @param n The number of draws desired from each chain after warmup.
#' @param poly An `mpoly` object.
#' @param sd The "standard deviation" component of the normal kernel.
#' @param vars A character vector of the indeterminates in the distribution.
#' @param w A named list of box constraints for vectors to be passed to Stan,
#'   see examples. If a single number, a box window `(-w, w)` is applied to all
#'   variables.
#' @param output Either `"simple"` or `"tibble"` output format.
#' @param dist Either `"norm"` (normal) or `"unif"` (uniform).
#' @param homo If `TRUE`, sampling is done from a homoskedastic variety normal
#'   distribution.
#' @param correct_p_coefficients If `TRUE`, normalize polynomial coefficients.
#' @param correct_dp_coefficients If `TRUE`, normalize derivative coefficients.
#' @param message If `TRUE`, print progress messages showing remaining samples.
#'
#' @return A matrix or tibble containing the accepted samples.



rejection_sampler <- function(n,
                              poly,
                              sd = .01,
                              vars = sort(mpoly::vars(poly)),
                              w = 1.25,
                              output = "simple",
                              dist = c("norm", "unif"),
                              homo = TRUE,
                              correct_p_coefficients = FALSE,
                              correct_dp_coefficients = FALSE,
                              message = FALSE) {
  # Expand scalar/length-2 window specs to per-variable bounds.
  n_vars <- length(vars)
  if (is.numeric(w) && length(w) == 1) {
    w <- replicate(length(vars), c(-w, w), simplify = FALSE)
    names(w) <- vars
  } else if (is.numeric(w) && length(w) == 2) {
    w <- replicate(length(vars), w, simplify = FALSE)
    names(w) <- vars
  }
  dist <- match.arg(dist)
  if (correct_p_coefficients) poly <- normalize_coefficients(poly)

  if (is.mpolyList(poly)) {
    # mpolyList case: g(x) is vector-valued and uses Jacobian-based scaling.
    n_polys <- length(poly)
    pf <- as.function(poly, varorder = vars, silent = TRUE)
    dp <- dpfs <- vector(mode = "list", length = n_polys)
    for (i in seq_len(n_polys)) {
      dp[[i]] <- deriv(poly[[i]], vars)
      if (correct_dp_coefficients) dp[[i]] <- normalize_coefficients(dp[[i]])
      dpfs[[i]] <- as.function(dp[[i]], varorder = vars, silent = TRUE)
    }

    dpf <- function(x) {
      mat <- matrix(NA_real_, nrow = n_polys, ncol = n_vars)
      for (i in seq_len(n_polys)) mat[i, ] <- dpfs[[i]](x)
      mat
    }

    if (is.vector(sd)) {
      # Scalar/diagonal scale case.
      if (homo) {
        pbar <- function(x) {
          g <- pf(x)
          J <- dpf(x)
          as.numeric(g %*% solve(tcrossprod(J), g))
        }
      } else {
        pbar <- function(x) pf(x)^2
      }
      log_ptilde <- function(x) -pbar(x) / (2 * sd^2)
    } else if (is.matrix(sd)) {
      # Full covariance case via spectral decomposition.
      eig_sd <- eigen(sd, symmetric = TRUE)
      la <- eig_sd$values
      P <- eig_sd$vectors
      la_inv <- 1 / la
      sqrt_sd_inv <- diag(sqrt(la_inv)) %*% t(P)

      if (homo) {
        pbar <- function(x) {
          g <- pf(x)
          J <- dpf(x)
          if (n_vars == n_polys) {
            J_inv <- solve(J)
          } else {
            J_inv <- MASS::ginv(J)
          }
          as.numeric(crossprod(sqrt_sd_inv %*% J_inv %*% g))
        }
      } else {
        pbar <- function(x) pf(x)^2
      }

      log_ptilde <- function(x) -pbar(x) / 2
    }
  } else if (is.mpoly(poly)) {
    # Single-polynomial case.
    pf <- as.function(poly, varorder = vars, silent = TRUE)
    dp <- deriv(poly, vars)
    if (correct_dp_coefficients) dp <- normalize_coefficients(dp)

    if (n_vars > 1) {
      ssdp <- mp("0")
      for (i in seq_len(n_vars)) {
        ssdp <- ssdp + dp[[i]]^2
      }
    } else {
      ssdp <- dp^2
    }

    ssdpf <- if (is.constant(ssdp)) function(x) ssdp[[1]][["coef"]] else as.function(ssdp, varorder = vars, silent = TRUE)

    if (homo) {
      # Normalize by gradient magnitude for approximate arc-length scaling.
      pbar <- function(x) pf(x) / sqrt(ssdpf(x))
      log_ptilde <- function(x) {
        -exp((2 * log(abs(pf(x))) - log(ssdpf(x))) - log(2) - 2 * log(sd))
      }
    } else {
      pbar <- pf
      log_ptilde <- function(x) {
        -exp((2 * log(abs(pf(x)))) - log(2) - 2 * log(sd))
      }
    }
  } else {
    stop("`poly` should be a `mpoly` or `mpolyList` object.", call. = FALSE)
  }

  mat <- matrix(nrow = 0, ncol = n_vars, dimnames = list(NULL, vars))
  n_remaining <- n
  while (n_remaining > 0) {
    # Propose points uniformly inside the box, then accept/reject.
    if (message) cat("\r", strrep(" ", 80))
    if (message) cat("\r", scales::number_format(big.mark = ",")(n_remaining), " remaining...", sep = "")

    u <- matrix(nrow = n_remaining, ncol = n_vars, dimnames = list(NULL, vars))
    for (i in seq_len(n_vars)) {
      u[, i] <- runif(n_remaining, w[[i]][1], w[[i]][2])
    }

    if (dist == "norm") {
      gbars <- apply(u, 1, log_ptilde)
      accept_reject_ndcs <- which(log(runif(n_remaining)) <= gbars)
    } else if (dist == "unif") {
      gbars <- apply(u, 1, pbar)
      accept_reject_ndcs <- which(abs(gbars) <= sd)
    }
    n_accept <- length(accept_reject_ndcs)

    if (n_accept > 0) {
      samp <- u[accept_reject_ndcs, , drop = FALSE]
      mat <- rbind(mat, samp)
      n_remaining <- n_remaining - n_accept
    }
  }
  cat("\r", strrep(" ", 80))
  cat("\r")
  row.names(mat) <- NULL
  if (output == "tibble") tibble::as_tibble(mat) else mat

}
