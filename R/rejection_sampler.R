#' Rejection sampler for variety normal distribution
#'
#'
#' This function performs rejection sampling to generate samples from a
#' variety normal distribution.
#'
#' @param n The number of draws desired from each chain after warmup.
#' @param poly An mpoly object.
#' @param sd The "standard deviation" component of the normal kernel.
#' @param vars A character vector of the indeterminates in the distribution.
#' @param w A named list of box constraints for vectors to be passed to Stan,
#'   see examples. A If a single number, a box window (-w,w) is applied to all
#'   variables.
#' @param dist  either `"norm"` (normal) or `"unif"` (uniform).
#' @param homo If \code{TRUE}, the sampling is done from homoskedastic variety
#'  normal distribution.
#' @param correct_p_coefficients If \code{TRUE} normalize the coefficients of the polynomial.
#' @param correct_dp_coefficients If \code{TRUE} to normalize the coefficients of the polynomial's derivative.
#' @param message If \code{TRUE}, message the user about remaining numbers of sampling
#'
#' @return A matrix or tibble containing the accepted samples.



rejection_sampler <- function(n,
                              poly,
                              sd = .01,
                              vars = sort(mpoly::vars(poly)),
                              w = 1.25,
                              output = "simple",
                              dist = c("norm", "unif"),
                              homo= TRUE,
                              correct_p_coefficients = FALSE,
                              correct_dp_coefficients = FALSE,
                              message = FALSE) {



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
    n_polys <- length(poly)
    pf <- as.function(poly, varorder = vars, silent = TRUE)
    dp <- dpfs <- vector(mode = "list", length = n_polys)
    for (i in 1:n_polys) {
      dp[[i]] <- deriv(poly[[i]], vars)
      if (correct_dp_coefficients) dp[[i]] <- normalize_coefficients(dp[[i]])
      dpfs[[i]] <- as.function(dp[[i]], varorder = vars, silent = TRUE)
    }

    dpf <- function(x) {
      mat <- matrix(NA_real_, nrow = n_polys, ncol = n_vars)
      for (i in 1:n_polys) mat[i, ] <- dpfs[[i]](x)
      mat
    }

    if (is.vector(sd)) {
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
      ptilde <- function(x) exp(-pbar(x) / (2 * sd^2))
    } else if (is.matrix(sd)) {
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
      ptilde <- function(x) exp(-pbar(x) / 2)
    }
  } else if (is.mpoly(poly)) {
    pf <- as.function(poly, varorder = vars, silent = TRUE)
    dp <- deriv(poly, vars)
    if (correct_dp_coefficients) dp <- normalize_coefficients(dp)

    if (n_vars > 1) {
      ssdp <- mp("0")
      for (i in 1:n_vars) ssdp <- ssdp + dp[[i]]^2
    } else {
      ssdp <- dp^2
    }

    ssdpf <- if (is.constant(ssdp)) function(x) ssdp[[1]][["coef"]] else as.function(ssdp, varorder = vars, silent = TRUE)

    if (homo) {
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
    ptilde <- function(x) exp(-pbar(x)^2 / (2 * sd^2))
  } else stop("`poly` should be a `mpoly` or `mpolyList` object.", call. = FALSE)

  mat <- matrix(nrow = 0, ncol = n_vars, dimnames = list(NULL, vars))
  n_remaining <- n
  while (n_remaining > 0) {
    if (message) cat("\r", strrep(" ", 80))
    if (message) cat("\r", scales::number_format(big.mark = ",")(n_remaining), " remaining...", sep = "")

    u <- matrix(nrow = n_remaining, ncol = n_vars, dimnames = list(NULL, vars))
    for (i in 1:n_vars) {
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
  out <- if (output == "tibble") tibble::as_tibble(mat) else mat
  out

}

