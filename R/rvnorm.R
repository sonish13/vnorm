#' The Variety Normal Distribution
#'
#' Un-normalized density function and random generation for the variety normal
#' distribution with mean equal to \code{poly} and "standard deviation" equal to
#' \code{sd}. Please see details for caveats.
#'
#' If the variety you are interested in is connected, this strategy should work
#' well out of the box.  If it isn't, you'll likely need to rely on running
#' multiple chains, and it is very likely, if not probable, that the sampling
#' will be biased to one or more of those components and down-sample others.
#' Question: what is the relative likelihood of each component, or an equal unit
#' of length, on different components? How does this generalize to more
#' varieties of varying dimensions?
#'
#' @param n The number of draws desired from each chain after warmup.
#' @param poly An mpoly object.
#' @param sd The "standard deviation" component of the normal kernel.
#' @param output \code{"simple"}, \code{"tibble"}, \code{"stanfit"}.
#' @param chains The number of chains to run for the random number generation,
#'   see [stan()].
#' @param cores The number of CPU cores to distribute the chains across, see
#'   [stan()].
#' @param warmup Number of warmup iterations in [stan()].
#' @param inc_warmup If \code{TRUE}, the MCMC warmup steps are included in the
#'   output.
#' @param thin [stan()] \code{thin} parameter.
#' @param normalized If \code{TRUE}, the polynomial is gradient-normalized. This
#'   is highly recommended.
#' @param inject_direct Directly specify printed polynomial to string inject
#'   into the stan code. Requires you specify \code{vars}, \code{numerator}, and
#'   \code{denominator}.
#' @param verbose \code{TRUE} or \code{FALSE}; determines level of messaging.
#' @param vars A character vector of the indeterminates in the distribution.
#' @param numerator,denominator A character(1) containing the printed numerator
#'   of the variety normal distribution.
#' @param w A named list of box constraints for vectors to be passed to Stan,
#'   see examples. A If a single number, a box window (-w,w) is applied to all
#'   variables.
#' @param refresh The \code{refresh} argument of [stan()], which governs how
#'   much information is provided to the user while sampling.
#' @param pre_compiled Whether to use pre-compiled stan models or not. Available
#'   for polynomials with three indeterminates and three degrees. Defaults to
#'   \code{TRUE}.
#' @param code_only If \code{TRUE}, will only formulate and return Stan code.
#' @param ... Additional parameters to pass to [stan()].
#' @name rvnorm
#' @return Either (1) matrix whose rows are the individual draws from the
#'   distribution, (2) a [tbl_df-class] object with the draws along with
#'   additional information, or (3) an object of class [stanfit-class].
#' @author David Kahle
#' @examples
#'
#' \dontrun{
#'
#' library("tidyverse")
#' options("mc.cores" = parallel::detectCores() - 1)
#'
#' ## basic usage
#' ########################################
#'
#' # single polynomial
#' p <- mpoly::mp("x^2 + y^2 - 1")
#' samps <- rvnorm(1000, p, sd = .05, w = 2)
#' head(samps)
#' str(samps)
#' plot(samps, asp = 1)
#' plot(p, add = TRUE)
#'
#' # returning a data frame
#' (samps <- rvnorm(5000, p, sd = .05, w = 2, output = "tibble"))
#'
#' ggplot(samps, aes(x, y)) +
#'   geom_point(size = .5) +
#'   coord_equal()
#'
#' ggplot(samps, aes(x, y)) +
#'   geom_point(size = .5) +
#'   geom_variety(poly = p) +
#'   coord_equal()
#'
#' ggplot(samps, aes(x, y)) +
#'   geom_point(aes(color = g), size = .5) +
#'   geom_variety(poly = p) +
#'   scale_color_gradient2() +
#'   coord_equal()
#'
#' ggplot(samps, aes(x, y)) +
#'   stat_density2d(
#'     aes(fill = after_stat(density)),
#'     geom = "raster", contour = FALSE
#'    ) +
#'   geom_variety(poly = p) +
#'   coord_equal()
#'
#' library("ggdensity")
#' ggplot(samps, aes(x, y)) +
#'   geom_hdr(xlim = c(-2,2), ylim = c(-2,2)) +
#'   geom_variety(poly = p) +
#'   coord_equal()
#'
#'
#' # in three variables
#' (samps <- rvnorm(20, mp("x^2 + y^2 + z^2 - 1"), sd = .05, w = 2))
#' apply(samps, 1, function(v) sqrt(sum(v^2)))
#'
#'
#' # more than one polynomial, # vars > # eqns, underdetermined system
#' p <- mp(c("x^2 + y^2 + z^2 - 1", "z"))
#' (samps <- rvnorm(500, p, sd = .1, output = "tibble"))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#'
#' ggplot(samps, aes(x, y, color = `g[1]`)) + geom_point() +
#'   scale_color_gradient2(mid = "gray80") + coord_equal()
#'
#' ggplot(samps, aes(x, y, color = `g[2]`)) + geom_point() +
#'   scale_color_gradient2(mid = "gray80") + coord_equal()
#'
#' ggplot(samps, aes(x, z, color = `g[1]`)) + geom_point() +
#'   scale_color_gradient2(mid = "gray80") + coord_equal()
#'
#'
#'
#' # more than one polynomial, # vars < # eqns, overdetermined system
#' p <- mp(c("3 x", "3 y", "2 x + 2 y", "3 (x^2 + y)", "3 (x^2 - y)"))
#' (samps <- rvnorm(500, p, sd = .1, output = "tibble"))
#'
#' samps |>
#'   select(x, y, starts_with("g")) |>
#'   pivot_longer(starts_with("g"), "equation", "value") |>
#'   ggplot(aes(x, y, color = value)) + geom_point() +
#'     scale_color_gradient2(mid = "gray80") + coord_equal() +
#'     facet_wrap(~ equation)
#'
#'
#' ## using refresh to get more info
#' ########################################
#'
#' rvnorm(2000, p, sd = .1, "tibble", verbose = TRUE)
#' rvnorm(2000, p, sd = .1, "tibble", refresh = 500)
#' rvnorm(2000, p, sd = .1, "tibble", refresh = 0) # default
#' rvnorm(2000, p, sd = .1, "tibble", refresh = -1)
#'
#'
#'
#' ## many chains in parallel
#' ########################################
#'
#' options(mc.cores = parallel::detectCores())
#' p <- mp("x^2 + (4 y)^2 - 1")
#' (samps <- rvnorm(1e4, p, sd = .01, "tibble", verbose = TRUE, chains = 8))
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .01*c(1,1)) + coord_equal()
#' # decrease sd to get more uniform sampling
#'
#'
#' ## windowing for unbounded varieties
#' ########################################
#' # windowing is needed for unbounded varieties
#' # in the following, look at the parameters block
#'
#' p <- mp("x y - 1") # unbounded variety, 1 poly
#' p <- mp(c("x y - 1", "y - x")) # 2 polys
#'
#' rvnorm(1e3, p, sd = .01, "tibble", code_only = TRUE)
#' rvnorm(1e3, p, sd = .01, "tibble", code_only = TRUE, w = 1.15)
#'
#' window <- list("x" = c(-1.5, 1.25), "y" = c(-2, 1.5))
#' rvnorm(1e3, p, sd = .01, "tibble", code_only = TRUE, w = window)
#'
#' window <- list("x" = c(-1.5, 1.5))
#' rvnorm(1e3, p, sd = .01, "tibble", code_only = TRUE, w = window)
#'
#'
#'
#' ## the importance of normalizing
#' ########################################
#' # one of the effects of the normalizing is to stabilize variances, making
#' # them roughly equivalent globally over the variety.
#'
#' # lemniscate of bernoulli
#' p <- mp("(x^2 + y^2)^2 - 2 (x^2 - y^2)")
#'
#' # normalized, good
#' (samps <- rvnorm(2000, p, .05, "tibble"))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#'
#' # unnormalized, bad
#' (samps <- rvnorm(2000, p, .05, "tibble", normalized = FALSE))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#'
#'
#'
#' ## semi-algebraic sets
#' ########################################
#' # inside the semialgebraic set x^2 + y^2 <= 1
#' # this is the same as x^2 + y^2 - 1 <= 0, so that
#' # x^2 + y^2 - 1 + s^2 == 0 for some slack variable s
#' # this is the projection of the sphere into the xy-plane.
#'
#' p <- mp("1 - (x^2 + y^2) - s^2")
#' samps <- rvnorm(1e4, p, sd = .1, "tibble", chains = 8, refresh = 1e3)
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#'
#' ggplot(sample_n(samps, 2e3), aes(x, y, color = s)) +
#'   geom_point(size = .5) +
#'   scale_color_gradient2() +
#'   coord_equal()
#'
#' # alternative representation
#' # x^2 + y^2 - 1 <= 0 iff s^2 (x^2 + y^2 - 1) + 1 == 0
#' # note that it's gradient is more complicated.
#' p <- mp("s^2 (x^2 + y^2 - 1) + 1")
#' samps <- rvnorm(1e4, p, sd = .1, "tibble", chains = 8, w = 2, refresh = 1e3)
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#'
#'
#' ## keeping the warmup / the importance of multiple chains
#' ########################################
#'
#' p <- mp("((x + 1.5)^2 + y^2 - 1) ((x - 1.5)^2 + y^2 - 1)")
#' ggvariety(p, xlim = c(-3,3)) + coord_equal()
#'
#' # notice the migration of chains initialized away from the distribution
#' # (it helps to make the graphic large on your screen)
#' samps <- rvnorm(500, p, sd = .05, "tibble", chains = 8, inc_warmup = TRUE)
#' ggplot(samps, aes(x, y, color = iter)) +
#'   geom_point(size = 1, alpha = .5) + geom_path(alpha = .2) +
#'   coord_equal() + facet_wrap(~ factor(chain))
#'
#' samps <- rvnorm(2500, p, sd = .05, "tibble", chains = 8, inc_warmup = TRUE)
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) +
#'   coord_equal() + facet_wrap(~ factor(chain))
#'
#'
#' ## ideal-variety correspondence considerations
#' ########################################
#'
#' p <- mp("x^2 + y^2 - 1")
#'
#' samps_1 <- rvnorm(250, p^1, sd = .1, output = "tibble", chains = 8)
#' samps_2 <- rvnorm(250, p^2, sd = .1, output = "tibble", chains = 8)
#' # samps_3 <- rvnorm(250, p^3, sd = .1, output = "tibble", chains = 8)
#' # samps_4 <- rvnorm(250, p^4, sd = .1, output = "tibble", chains = 8)
#' samps <- bind_rows(mget(apropos("samps_[1-4]")))
#' samps$power <- rep(seq_along(apropos("samps_[1-4]")), each = 2000)
#'
#' ggplot(samps, aes(x, y, color = g < 0)) +
#'   geom_point(size = .5) +
#'   coord_equal(xlim = c(-3,3), ylim = c(-3,3)) +
#'   facet_wrap(~ power)
#'
#'
#' ## neat examples
#' ########################################
#' # an implicit Lissajous region, view in separate window large
#'
#' # x = cos(m t + p)
#' # y = sin(n t + q)
#' (p <- lissajous(3, 2,  -pi/2, 0))
#' (p <- lissajous(4, 3,  -pi/2, 0))
#' (p <- lissajous(5, 4,  -pi/2, 0))
#' (p <- lissajous(3, 3,  0, 0))
#' (p <- lissajous(5, 5,  0, 0))
#' (p <- lissajous(7, 7,  0, 0))
#' ggplot() +
#'   geom_variety(poly = p, xlim = c(-1,1), ylim = c(-1,1)) +
#'   coord_equal()
#'
#' p <- plug(p, "x", mp(".5 x"))
#' p <- plug(p, "y", mp(".5 y"))
#'
#' # algebraic set
#' samps <- rvnorm(5e3, p, sd = .01, "tibble", chains = 8, cores = 8)
#' ggplot(samps, aes(x, y, color = factor(.chain))) +
#'   geom_point(size = .5) +
#'   coord_equal()
#'
#' ggplot(samps, aes(x, y)) + geom_hdr() + coord_equal()
#'
#' ggplot(samps, aes(x, y, color = factor(.chain))) +
#'   geom_point(size = .5) +
#'   coord_equal() +
#'   facet_wrap(~ factor(.chain))
#'
#' # semi-algebraic set
#' samps_normd <- rvnorm(1e4, p + mp("s^2"), sd = .01, "tibble", chains = 8,
#'   cores = 8, normalized = TRUE
#' )
#' samps_unormd <- rvnorm(1e4, p + mp("s^2"), sd = .01, "tibble", chains = 8,
#'   cores = 8, normalized = FALSE
#' )
#'
#' bind_rows(
#'   samps_normd  |> mutate(normd = TRUE),
#'   samps_unormd |> mutate(normd = FALSE)
#' ) |>
#'   ggplot(aes(x, y)) +
#'     geom_point(size = .5) +
#'     # geom_bin2d(binwidth = .05*c(1,1)) +
#'     facet_grid(normd ~ .chain) +
#'     coord_equal()
#'
#'
#' ggplot(samps_normd, aes(x, y)) +
#'   geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#'
#' }
#'













#' @rdname rvnorm
#' @export
rvnorm <- function(
    n,
    poly,
    sd,
    output = "simple",
    chains = 4L,
    warmup = max(500, floor(n/2)),
    inc_warmup = FALSE,
    thin = 1L,
    inject_direct = FALSE,
    verbose = FALSE,
    cores = min(chains, getOption("mc.cores", 1L)),
    normalized = TRUE,
    w,
    vars,
    numerator,
    denominator,
    refresh = 0L,
    code_only = FALSE,
    pre_compiled = TRUE,
    ...
) {

  if ( is.character(poly) ) poly <- mp(poly)
  if ( is.mpoly(poly) ) {
    n_eqs <- 1L
  } else if ( is.mpolyList(poly) ) {
    n_eqs <- length(poly)
  } else {
    stop("`poly` should be either a character vector, mpoly, or mpolyList.", call. = FALSE)
  }
  deg <- mpoly::totaldeg(poly)
  num_of_vars <- length(mpoly::vars(poly))
  if (refresh) if (verbose) refresh <- max(ceiling(n/10), 1L) else refresh <- 0L
  if (!missing(refresh)) stopifnot(is.numeric(refresh), length(refresh) == 1L)

  if (n_eqs > 1 | length(mpoly::vars(poly)) > 2 | mpoly::totaldeg(poly) > 3) pre_compiled = FALSE
  if(code_only) return(create_stan_code(poly, sd, n_eqs, w, normalized))
  #if(!normalized) pre_compiled = FALSE # This line will be removed after normalized stan files are added

  if (pre_compiled) {
    list_for_transformation <- check_and_replace_vars(poly)
    if(!(list_for_transformation$mapping |> unlist() |> is.null())) {
      poly <- list_for_transformation$polynomial
      var_info <- list_for_transformation$mapping
      output_needs_rewriting <- TRUE
    }
    if (normalized) {
      stan_file_name <- paste(num_of_vars, deg, "norm", sep = "_")
    } else {
      stan_file_name <- paste(num_of_vars, deg, sep = "_")
    }
    model <- stan_package_model(name = stan_file_name, package = "vnorm")
    stan_data <- make_coeficients_data(poly, num_of_vars, deg)
    stan_data <- c(stan_data, "w" = w, "si" = sd)

  } else {
    stan_code <- create_stan_code(poly, sd, n_eqs, w, normalized, vars)
    stan_file <- write_stan_file(stan_code)
    model <- cmdstan_model(stan_file)
    stan_data <- list("si" = sd)
  }

  # run stan sampler ---------------------------------------------------------
  samps <- model$sample(
    data = stan_data,
    refresh = refresh,
    iter_warmup = warmup,
    iter_sampling = ceiling(n/chains),
    chains = chains,
    parallel_chains = cores,
    adapt_delta = .999,
    max_treedepth = 20L,
  )

  # parse output and return -------------------------------------------------

  if(output == "simple") {
    df <- samps$draws(format = "df", inc_warmup = inc_warmup) |> as.data.frame()
    df <- df[(nrow(df)-n+1):nrow(df), mpoly::vars(poly)]
    row.names(df) <- NULL
    return( df )
  }

  if (output == "tibble") {
    df <- samps$draws(format = "df", inc_warmup = inc_warmup) |> tibble::as_tibble()
    df <- df[(nrow(df)-n+1):nrow(df),]
    row.names(df) <- NULL
    df <- cbind(df[,-1], df[,1]) # move lp__ to last column
    return( as_tibble(df) )
  }

  if (output == "stanfit") return(samps)

}









mpoly_to_stan <- function (mpoly) {
  p <- get("print.mpoly", asNamespace("mpoly"))
  p(mpoly, stars = TRUE, silent = TRUE, plus_pad = 0L, times_pad = 0L) |>
    stringr::str_replace_all("[*]{2}", "^")
}

mpolyList_to_stan <- function (mpolyList) {
  p <- get("print.mpolyList", asNamespace("mpoly"))
  p(mpolyList, silent = TRUE, stars = TRUE, plus_pad = 0, times_pad = 0) |>
    stringr::str_replace_all("\\*\\*", "^") |>
    stringr::str_c(collapse = ", ")
}

create_stan_code <- function(poly, sd, n_eqs, w, normalized, vars) {
  d <- get("deriv.mpoly", asNamespace("mpoly"))
  if (n_eqs == 1L) {
    # single polynomial provided



    if (!is.mpoly(poly))
      poly <- mp(poly)
    if (missing(vars))
      vars <- mpoly::vars(poly)
    n_vars <- length(vars)
    reorder.mpoly <- get("reorder.mpoly", asNamespace("mpoly"))
    poly <- reorder.mpoly(poly, varorder = sort(vars))

    g_string <- mpoly_to_stan(poly)

    if (normalized) {
      if (n_vars > 1) {
        grad <- deriv(poly, var = mpoly::vars(poly))
        ndg_sq <- Reduce(`+`, grad ^ 2)
      } else {
        ndg_sq <- gradient(poly) ^ 2
      }
      ndg_sq_string <- mpoly_to_stan(ndg_sq)
      ndg_string <- glue::glue("sqrt({ndg_sq_string})")
    } else {
      ndg_string <- "1"
    }

    # set variables
    if (missing(w)) {
      parms <- glue::glue("real {vars};")
    } else {
      if (is.numeric(w) && length(w) == 1L) {
        parms <- glue::glue("real<lower=-{w},upper={w}> {vars};")
      } else if (is.list(w)) {
        stopifnot(all(names(w) %in% vars))
        parms <- vector("character", n_vars)
        for (var_ndx in seq_along(vars)) {
          parms[var_ndx] <- if (vars[var_ndx] %in% names(w)) {
            var_ndx_in_w <- which(names(w) == vars[var_ndx])
            glue::glue(
              "real<lower={w[[var_ndx_in_w]][1]},upper={w[[var_ndx_in_w]][2]}> {vars[var_ndx]};"
            )
          } else {
            glue::glue("real {vars[var_ndx]};")
          }
        }
      } else {
        stop("bound parameter misspecified, see ?rvnorm.", call. = FALSE)
      }
    }
    parms <- parms |> str_c(collapse = "\n  ")

    stan_code <- glue::glue(
      "
data {
  real<lower=0> si;
}

parameters {
  {{parms}}
}

transformed parameters {
  real g = {{g_string}};
  real ndg = {{ndg_string}};
}

model {
  target += normal_lpdf(0.00 | g/ndg, si);
}
      ", .open = "{{", .close = "}}"
    )


  } else {


    # multiple polynomials provided

    vars <- mpoly::vars(poly)
    n_vars <- length(vars)
    # if (n_eqs > n_vars) stop("Overdetermined systems not yet supported.")


    printed_polys <- mpolyList_to_stan(poly)

    printed_jac <- array("", dim = c(n_eqs, n_vars))
    for (i in 1:n_eqs) {
      for (j in 1:n_vars) {
        if (normalized) {
          printed_jac[i,j] <- mpoly_to_stan(d(poly[[i]], vars[j]))
        } else {
          printed_jac[i,j] <- if (i == j) "1" else "0"
        }
      }
    }
    printed_jac <- printed_jac |>
      apply(1L, str_c, collapse = ", ") |>
      str_c("      [", ., "]", collapse = ", \n") |>
      str_c("[\n", ., "\n    ]") |>
      str_replace_all("\\*\\*", "^")


    # set variables
    if (missing(w)) {
      parms <- glue::glue("real {vars};")
    } else {
      if (is.numeric(w) && length(w) == 1L) {
        parms <- glue::glue("real<lower=-{w},upper={w}> {vars};")
      } else if (is.list(w)) {
        stopifnot(all(names(w) %in% vars))
        parms <- vector("character", n_vars)
        for (var_ndx in seq_along(vars)) {
          parms[var_ndx] <- if (vars[var_ndx] %in% names(w)) {
            var_ndx_in_w <- which(names(w) == vars[var_ndx])
            glue::glue("real<lower={w[[var_ndx_in_w]][1]},upper={w[[var_ndx_in_w]][2]}> {vars[var_ndx]};")
          } else {
            glue::glue("real {vars[var_ndx]};")
          }
        }
      } else {
        stop("bound parameter misspecified, see ?rvnorm.", call. = FALSE)
      }
    }
    parms <- parms |> str_c(collapse = "\n  ")

    # normalization matrix
    if (is.numeric(sd) && is.vector(sd) && length(sd) == 1L) {

      gbar_string <- if (n_vars == n_eqs) {
        "J \\ g"
      } else if (n_vars > n_eqs) {
        "J' * ((J*J') \\ g)"
      } else {
        "(J'*J) \\ (J'*g)"
      }

    } else stop("This sd not yet supported.")


    # write stan code
    stan_code <- glue::glue("
data {
  real<lower=0> si;
}

parameters {
  {{parms}}
}

transformed parameters {
  vector[{{n_eqs}}] g = [{{printed_polys}}]';
  matrix[{{n_eqs}},{{n_vars}}] J = {{printed_jac}};
}

model {
  target += normal_lpdf(0.00 | {{gbar_string}}, si);
}
      ",  .open = "{{", .close = "}}"
    )
  }


  stan_code
}



make_coeficients_data <- function(poly, num_of_vars ,deg, basis = c("x", "y","z")) {
  required_coefs<- basis_monomials(basis[seq_along(1:num_of_vars)], deg) |>
    lapply(reorder,varorder = basis) |>
    lapply(coef) |>
    unlist() |>
    get_listed_coeficients() |>
    lapply(function(x) 0)
  available_coef <- get_listed_coeficients( coef(poly) )
  required_coefs[names(available_coef)] <- available_coef
  required_coefs
}



get_listed_coeficients <- function(coefs) {
  convert_names <- function(term) {
    term <- gsub("\\s+", "", term)  # Remove spaces
    if (term == "1") return("b0")   # Constant term should be "b0"
    term <- gsub("\\^", "", term)   # Remove power symbol (^)
    return(paste0("b", term))       # Add "b" at the beginning
  }
  names(coefs) <- sapply(names(coefs), convert_names)
  as.list(coefs)
}



check_and_replace_vars <- function(p) {
  current_vars <- vars(p)
  num_vars <- length(current_vars)
  target_vars <- list(
    c("x"),           # for 1 indeterminate
    c("x", "y"),      # for 2 indeterminates
    c("x", "y", "z")  # for 3 indeterminates
  )
  if (num_vars > 3) {
    stop("The polynomial has more than 3 indeterminates.")
  }
  expected_vars <- target_vars[[num_vars]]
  if (setequal(current_vars, expected_vars)) {
    return(list(polynomial = p, mapping = list()))  # No replacement needed
  }
  var_mapping <- list()

  # Replace current variables with temporary placeholders to avoid conflicts
  temp_vars <- paste0("tmp", seq_along(current_vars))  # Temporary placeholders
  for (i in seq_along(current_vars)) {
    p <- swap(p, current_vars[i], temp_vars[i])
    var_mapping[[temp_vars[i]]] <- current_vars[i]  # Store the original variable
  }

  # Replace the temporary placeholders with the target variables (x, y, z)
  for (i in seq_along(expected_vars)) {
    p <- swap(p, temp_vars[i], expected_vars[i])
    var_mapping[[expected_vars[i]]] <- var_mapping[[temp_vars[i]]]  # Map to target variables
    var_mapping[[temp_vars[i]]] <- NULL  # Remove temp mapping
  }

  list(polynomial = p, mapping = var_mapping)
}




















