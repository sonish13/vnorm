#' Psuedo Probability Density Function for the Variety Normal Distribution
#'
#' Computes a pseudo probability density function (PDF) for the variety normal distribution.
#' This function supports both the heteroskedastic and homoskedastic cases.
#'
#' @param x A numeric vector of length equal to the number of variables in
#'   \code{poly}, or a numeric matrix (or data frame) with that many columns.
#' @param poly An \code{mpoly} object (univariate case) or an \code{mpolyList}
#' object (multivariate case).
#' @param sigma A positive numeric value representing the standard deviation for
#' univariate case. For the multivariate case,
#' may be a scalar, vector, or matrix. If \code{homo = TRUE}, \code{sigma}
#' must conform to the number of variables; if \code{homo = FALSE}, it must
#' conform to the number of polynomials.
#' @param homo Logical, default is `TRUE`. If `TRUE`, computes the homoskedastic
#'  variety normal density. If `FALSE`, computes the heteroskedastic variety
#'  normal density.
#' @param log Logical. If `TRUE`, returns the log of the density.
#' @return A numeric value (or vector) representing the density evaluated at `x`.
#'
#' @examples
#'
#' library(mpoly)
#'
#' ## Univariate usage
#' poly_uni <- mp("x")
#' pdvnorm(c(-1, 0, 1), poly_uni, sigma = 1)
#' pdvnorm(0, poly_uni, sigma = 2, log = TRUE)
#'
#' ## Multivariate (square) system: two polynomials in two variables
#' poly_sq <- mpolyList(mp("x"), mp("y"))
#' X <- rbind(c(0, 0), c(1, 2), c(-1, 3))
#'
#' ## Different sigma forms
#' pdvnorm(X, poly_sq, sigma = 1,        homo = TRUE)
#' pdvnorm(X, poly_sq, sigma = c(1, 2),  homo = TRUE)
#' pdvnorm(X, poly_sq, sigma = diag(c(1, 4)), homo = TRUE)
#'
#' ## Multivariate (underdetermined): one polynomial in two variables
#' poly_under <- mpolyList(mp("x + y"))
#' X_under <- rbind(c(1, 1), c(2, -1), c(0, 3))
#' pdvnorm(X_under, poly_under, sigma = diag(2),   homo = TRUE)
#' pdvnorm(X_under, poly_under, sigma = 1,         homo = FALSE)
#'
#' ## Multivariate (overdetermined): three polynomials in two variables
#' poly_over <- mpolyList(mp("x"), mp("y"), mp("x + y"))
#' X_over <- rbind(c(1, 2), c(0, -1), c(2, 2))
#' pdvnorm(X_over, poly_over, sigma = diag(2),    homo = TRUE)
#' pdvnorm(X_over, poly_over, sigma = c(1, 2, 3), homo = FALSE)
#'
#' @export
pdvnorm <- function(x, poly, sigma, homo = TRUE, log = FALSE) {
  # ---- Unified type check ---------------------------------------------------
  is_uni   <- inherits(poly, "mpoly")
  is_multi <- inherits(poly, "mpolyList")
  if (!(is_uni || is_multi)) {
    stop("'poly' must be either an 'mpoly' (univariate) or an 'mpolyList' (multivariate).")
  }

  # ===========================================================================
  # UNIVARIATE (mpoly)
  # ===========================================================================
  if (is_uni) {
    # here 'sigma' acts as sd
    if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0) {
      stop("For univariate case, 'sigma' must be a single positive numeric (sd).")
    }
    sd <- as.numeric(sigma)

    n <- length(mpoly::vars(poly))  # number of variables
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.vector(x)) {
      if (n == 1) {
        x <- matrix(x, ncol = 1)
      } else {
        if (length(x) != n) stop("x must be length n or a matrix with n columns.")
        x <- matrix(x, nrow = 1)
      }
    } else if (!(is.matrix(x) && ncol(x) == n)) {
      stop("x must be length n or a matrix with n columns.")
    }

    g_func <- suppressMessages(as.function(poly))

    if (homo) {
      grad_g_obj <- mpoly::gradient(poly)
      if (mean(is.constant(grad_g_obj)) == 1) {
        const <- unname(unlist(grad_g_obj)); force(const)
        grad_g <- function(...) const
      } else {
        grad_g <- suppressMessages(as.function(grad_g_obj))
      }
    }

    log_density <- apply(x, 1, function(row_vec) {
      g_val <- g_func(row_vec)
      if (homo) {
        grad_g_val <- sqrt(sum(grad_g(row_vec)^2))
        if (grad_g_val == 0) return(Inf)  # preserve your original behavior
        g_val <- g_val / grad_g_val
      }
      z <- g_val / sd
      -(0.5 * z^2 + base::log(sd) + 0.5 * base::log(2 * base::pi))
    })
    return(if (log) log_density else exp(log_density))
  }

  # ===========================================================================
  # MULTIVARIATE (mpolyList)
  # ===========================================================================
  vars <- mpoly::vars(poly)
  n <- length(vars)     # number of variables
  m <- length(poly)     # number of polynomials

  # Coerce x to matrix with n columns
  X <- if (is.null(dim(x))) matrix(as.numeric(x), nrow = 1) else as.matrix(x)
  if (!is.numeric(X) || any(!is.finite(X))) stop("'x' must be finite numeric.")
  if (ncol(X) != n) stop("'x' must have length n (vector) or n columns (matrix).")

  # Sigma handling
  if (length(sigma) == 1) {
    Sigma <- if (homo) diag(sigma, n) else diag(sigma, m)
  } else if (is.vector(sigma)) {
    if (homo && length(sigma) != n) stop("When homo=TRUE, length(sigma) must be n.")
    if (!homo && length(sigma) != m) stop("When homo=FALSE, length(sigma) must be m.")
    Sigma <- diag(sigma)
  } else {
    Sigma <- as.matrix(sigma)
    if (homo) {
      if (!all(dim(Sigma) == c(n, n))) stop("'sigma' matrix must be n x n when homo=TRUE.")
    } else {
      if (!all(dim(Sigma) == c(m, m))) stop("'sigma' matrix must be m x m when homo=FALSE.")
    }
  }

  # g(x | beta) for all rows
  g_fns <- suppressMessages(as.function(poly))
  g_vals_mat <- t(apply(X, 1, g_fns))
  if (m == 1) g_vals_mat <- matrix(g_vals_mat, ncol = 1)

  grad_fun <- vector("list", m)
  for (i in 1:m) {
    grad_fun[[i]] <- deriv(poly[[i]], var = vars)
    if (mean(is.constant(grad_fun[[i]])) == 1) {
      const <- unname(unlist(grad_fun[[i]]))
      force(const)
      grad_fun[[i]] <- function(...) {
        const
      }
    } else {
      grad_fun[[i]] <- suppressMessages(
        as.function(deriv(poly[[i]], var = vars), varorder = vars)
      )
    }
  }
  # --------------------------------------------------------------------------

  L <- chol(Sigma)
  log_det_sigma <- 2 * sum(log(diag(L)))

  out_log <- numeric(nrow(X))
  for (i in seq_len(nrow(X))) {
    xi <- as.numeric(X[i, ])
    g_vals <- as.numeric(g_vals_mat[i, ])

    if (homo) {
      # Build Jacobian J(x): m x n
      J <- matrix(NA_real_, nrow = m, ncol = n)
      for (j in seq_len(m)) J[j, ] <- grad_fun[[j]](xi)

      sv  <- svd(J)
      tol <- max(dim(J)) * .Machine$double.eps * ifelse(length(sv$d) > 0, sv$d[1], 0)
      r   <- sum(sv$d > tol)

      if (n > m && r == m) {
        # underdetermined but full row rank: J^+ = J' (J J')^{-1}
        Jp <- t(J) %*% solve(J %*% t(J))
      } else if (m > n && r == n) {
        # overdetermined but full column rank: J^+ = (J'J)^{-1} J'
        Jp <- solve(t(J) %*% J) %*% t(J)
      } else if (m == n && r == n) {
        # square full rank
        Jp <- solve(J)
      } else {
        # rank-deficient or numerically singular: fall back to SVD pseudoinverse
        dplus <- ifelse(sv$d > tol, 1 / sv$d, 0)
        Jp <- sv$v %*% (dplus * t(sv$u))  # V %*% diag(dplus) %*% t(U)
      }
      # ----------------------------------------------------------------------

      v <- Jp %*% g_vals
      q <- n
    } else {
      v <- g_vals
      q <- m
    }

    # Quadratic form via Cholesky: v' Sigma^{-1} v
    quad <- sum(backsolve(L, v, transpose = TRUE)^2)
    out_log[i] <- -0.5 * (q * log(2 * pi) + log_det_sigma + quad)
  }

  if (log) out_log else exp(out_log)
}
