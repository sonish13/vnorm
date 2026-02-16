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
#' @param sd The "standard deviation" component of the normal kernel. sd = Sigma
#'   when Sigma is specified.
#' @param Sigma Full "covariance" matrix or a vector specifying the diagonal of
#'   such a matrix.
#' @param output \code{"simple"}, \code{"tibble"}, \code{"stanfit"}.
#' @param rejection If \code{TRUE}, rejection sampling is used.
#' @param chains The number of chains to run for the random number generation,
#'   see \code{stan()}.
#' @param cores The number of CPU cores to distribute the chains across, see
#'   \code{stan()}.
#' @param warmup Number of warmup iterations in \code{stan()}.
#' @param inc_warmup If \code{TRUE}, the MCMC warmup steps are included in the
#'   output.
#' @param thin \code{stan()} \code{thin} parameter.
#' @param homo If \code{TRUE}, the sampling is done from homoskedastic variety
#'  normal distribution.
#' @param verbose \code{TRUE} or \code{FALSE}; determines level of messaging.
#' @param vars A character vector of the indeterminates in the distribution.
#' @param numerator,denominator A character(1) containing the printed numerator
#'   of the variety normal distribution.
#' @param w A named list of box constraints for vectors to be passed to Stan,
#'   see examples. A If a single number, a box window (-w,w) is applied to all
#'   variables.
#' @param refresh The \code{refresh} argument of \code{stan()}, which governs how
#'   much information is provided to the user while sampling.
#' @param pre_compiled Whether to use pre-compiled stan models or not. Available
#'   for polynomials with three indeterminates and three degrees. Defaults to
#'   \code{TRUE}.
#' @param code_only If \code{TRUE}, will only formulate and return Stan code.
#' @param ... Additional parameters to pass to \code{stan()}.
#' @param user_compiled If\code{TRUE}, user compiled stan program made using
#' [compile_stan_code] is used. Defaults to \code{FALSE}
#' @param show_messages If \code{TRUE}, stan sampler messages are shown.
#' Defaults to \code{FALSE}
#' @name rvnorm
#' @return Either (1) matrix whose rows are the individual draws from the
#'   distribution, (2) a [tbl_df-class] object with the draws along with
#'   additional information, or (3) an object of class \code{stanfit}.
#' @examples
#'
#'
#' library("tidyverse")
#' options("mc.cores" = parallel::detectCores() - 1)
#'
#' \dontrun{
#' ## basic usage
#' ########################################
#'
#' # single polynomial
#' p <- mp("x^2 + y^2 - 1")
#' samps <- rvnorm(1000, p, sd = .05)
#' head(samps)
#' str(samps)
#' plot(samps, asp = 1)
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
#'   geom_point(aes(color = factor(.chain))) +
#'   geom_variety(poly = p, linewidth = 1) +
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
#' # in three variables
#' (samps <- rvnorm(20, mp("x^2 + y^2 + z^2 - 1"), sd = .05, w = 2))
#' apply(samps, 1, function(v) sqrt(sum(v^2)))
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
#' # more than one polynomial, # vars < # eqns, overdetermined system
#' p <- mp(c("3 x", "3 y", "2 x + 2 y", "3 (x^2 + y)", "3 (x^2 - y)"))
#' (samps <- rvnorm(500, p, sd = .1, output = "tibble"))
#'
#' samps |>
#'   select(x, y, starts_with("g")) |>
#'   pivot_longer(starts_with("g"), names_to = "equation", values_to = "value") |>
#'   ggplot(aes(x, y, color = value)) + geom_point() +
#'     scale_color_gradient2(mid = "gray80") + coord_equal() +
#'     facet_wrap(~ equation)
#'
#' ## using refresh to get more info
#' ########################################
#'
#' rvnorm(2000, p, sd = .1, "tibble", verbose = TRUE)
#' rvnorm(2000, p, sd = .1, "tibble", refresh = 500)
#' rvnorm(2000, p, sd = .1, "tibble", refresh = 0) # default
#' rvnorm(2000, p, sd = .1, "tibble", refresh = -1)
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
#' ## windowing for unbounded varieties
#' ########################################
#' # windowing is needed for unbounded varieties
#' # in the following, look at the parameters block
#'
#' p <- mp("x y - 1") # unbounded variety, 1 poly
#' p <- mp(c("x y - 1", "y - x")) # 2 polys
#'
#' rvnorm(1e3, p, sd = .01, "tibble", code_only = TRUE)
#' rvnorm(1e3, p, sd = .01, "tibble",  w = 1.15)
#'
#' window <- list("x" = c(-1.5, 1.5))
#' rvnorm(1e3, p, sd = .01, "tibble",  w = window)
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
#' (samps <- rvnorm(2000, p, .05, "tibble", homo = FALSE))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
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
#' ## keeping the warmup / the importance of multiple chains
#' ########################################
#'
#' p <- mp("((x + 1.5)^2 + y^2 - 1) ((x - 1.5)^2 + y^2 - 1)")
#' ggplot() +
#'   geom_variety(poly = p, xlim = c(-3,3)) +
#'   theme(legend.position = "top") +
#'   coord_equal()
#'
#' # notice the migration of chains initialized away from the distribution
#' # (it helps to make the graphic large on your screen)
#' samps <- rvnorm(500, p, sd = .05, "tibble", chains = 8, inc_warmup = TRUE)
#'
#' ## ideal-variety correspondence considerations
#' ########################################
#'
#' p <- mp("x^2 + y^2 - 1")
#'
#' samps_1 <- rvnorm(250, p^1, sd = .1, output = "tibble", chains = 8)
#' samps_2 <- rvnorm(250, p^2, sd = .1, output = "tibble", chains = 8)
#' samps_3 <- rvnorm(250, p^3, sd = .1, output = "tibble", chains = 8)
#' samps_4 <- rvnorm(250, p^4, sd = .1, output = "tibble", chains = 8)
#' samps <- bind_rows(mget(apropos("samps_[1-4]")))
#' samps$power <- rep(seq_along(apropos("samps_[1-4]")), each = 250)
#'
#' ggplot(samps, aes(x, y, color = g < 0)) +
#'   geom_point(size = .5) +
#'   coord_equal(xlim = c(-3,3), ylim = c(-3,3)) +
#'   facet_wrap(~ power)
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
#'   geom_variety(poly = p, xlim = c(-2, 2), ylim = c(-2, 2), n = 201,
#'   show.legend = FALSE) +
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
#'   cores = 8, homo = TRUE
#' )
#' samps_unormd <- rvnorm(1e4, p + mp("s^2"), sd = .01, "tibble", chains = 8,
#'   cores = 8, homo = FALSE
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
#' ggplot(samps_normd, aes(x, y)) +
#'   geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#' }


#' @rdname rvnorm
#' @export
rvnorm <- function(n, poly, sd, output = "simple", Sigma, homo = TRUE,
                   w, vars, numerator, denominator, refresh = 0L,
                   code_only = FALSE, pre_compiled = TRUE, user_compiled = FALSE,
                   show_messages = FALSE, ...,
                   rejection = FALSE,
                   chains = 4L,
                   warmup = max(500, floor(n / 2)),
                   inc_warmup = FALSE,
                   thin = 1L,
                   verbose = FALSE,
                   cores = min(chains, getOption("mc.cores", 1L))
                   ) {

  # Initialization and checks
  if(rejection) {
    if(missing(w)) w = 1
    if (!(is.numeric(sd) && is.vector(sd) && length(sd) == 1L)) stop("This sd is not supported yet for rejection sampler.")
    sample <- rejection_sampler(n = n, poly = poly , sd = sd, vars = sort(mpoly::vars(poly)),
                                  w = w, output = output, homo = homo,
                                  message = verbose)
    return(sample)
  }
  output_needs_rewriting <- FALSE
  if (is.character(poly)) poly <- mp(poly)
  if (is.mpoly(poly)) {
    n_eqs <- 1L
  } else if (is.mpolyList(poly)) {
    n_eqs <- length(poly)
  } else {
    stop("`poly` should be either a character vector, mpoly, or mpolyList.", call. = FALSE)
  }
  if( !missing(Sigma) ) sd <- Sigma
  n_vars <- length(mpoly::vars(poly))
  if(length(sd) == 1){
    sd = sd
  } else if (is.vector(sd) & length(sd) == n_vars){
    sd = diag(sd)
  } else if (is.matrix(sd) & length(diag(sd)) == n_vars ) {
    sd = sd
  } else {
    stop("`Sigma` should be a number, vector of length equal to number of
         variables or matrix with diagonal length equal to number of variables.", call. = FALSE)
  }
  if (refresh) refresh <- if (verbose) max(ceiling(n / 10), 1L) else 0L
  if (n_eqs > 1) pre_compiled <- FALSE
  if (n_eqs > 1 || length(mpoly::vars(poly)) > 3 || base::max(mpoly::totaldeg(poly)) > 3)
    pre_compiled <- FALSE
  if (code_only) {
    stan_code <- create_stan_code(poly, sd, n_eqs, w, homo)
    return(stan_code)
  }

  # Model selection and data preparation
  if (user_compiled) { # incase we want use model that user compiled
    model_name <- generate_model_name(poly = poly, homo = homo, w = !missing(w))
    compiled_stan_info <- get_compiled_stan_info()
    if (nrow(compiled_stan_info) == 0) {
      stop(
        "No compiled model cache found. Run compile_stan_code() first.",
        call. = FALSE
      )
    }
    model_path <- compiled_stan_info$path[compiled_stan_info$name == model_name]
    if (length(model_path) == 0) {
      stop(
        "Requested compiled model not found in cache. ",
        "Run compile_stan_code() for this polynomial first.",
        call. = FALSE
      )
    }
    model_path <- model_path[1]
    model <- cmdstan_model(model_path)
    stan_data <- get_coefficeints_data(poly)
    stan_data <- if (missing(w)) c(stan_data, "si" = sd) else c(stan_data, "w" = w, "si" = sd)
  } else if (pre_compiled) { # data and model if/when pre-compiled model is used
    # replace variable if poly has anything other then x, y or z
    list_for_transformation <- check_and_replace_vars(poly)
    if (!is.null(unlist(list_for_transformation$mapping))) {
      poly <- list_for_transformation$polynomial
      var_info <- list_for_transformation$mapping
      output_needs_rewriting <- TRUE
    }
    # get name of stan model to use
    deg <- mpoly::totaldeg(poly)
    num_of_vars <- length(mpoly::vars(poly))
    stan_file_name <- paste(num_of_vars, deg, if (homo) "vn" else "hvn", sep = "_")
    if (!missing(w)) stan_file_name <- paste(stan_file_name, "w", sep = "_")
    model <- instantiate::stan_package_model(name = stan_file_name, package = "vnorm")
    stan_data <- make_coefficients_data(poly, num_of_vars, deg)
    stan_data <- if (missing(w)) c(stan_data, "si" = sd) else c(stan_data, "w" = w, "si" = sd)
  } else {
    # create code to run stan model from scratch
    stan_code <- create_stan_code(poly, sd, n_eqs, w, homo, vars)
    stan_file <- write_stan_file(stan_code)
    model <- cmdstan_model(stan_file)
    stan_data <- list("si" = sd)
  }

  # Run Stan sampler
  samps <- model$sample(data = stan_data, refresh = refresh, iter_warmup = warmup,
                        iter_sampling = ceiling(n / chains), chains = chains,
                        parallel_chains = cores, adapt_delta = .999, max_treedepth = 20L,
                        save_warmup = inc_warmup, show_messages = show_messages, ... = ...)

  # Parse output and return
  if (output == "simple") {
    df <- as.data.frame(samps$draws(format = "df", inc_warmup = inc_warmup))
    df <- df[(nrow(df) - n + 1):nrow(df), mpoly::vars(poly)]
    row.names(df) <- NULL
    if (output_needs_rewriting) df <- rename_output_df(df, replacement_list = var_info)
    return(df)
  }

  if (output == "tibble") {
    df <- tibble::as_tibble(samps$draws(format = "df", inc_warmup = inc_warmup))
    df <- samps$draws(format = "df", inc_warmup = inc_warmup)
    df <- df[(nrow(df) - n + 1):nrow(df), ]
    row.names(df) <- NULL
    df <- cbind(df[, -1], df[, 1])  # Move lp__ to last column
    if (output_needs_rewriting) df <- rename_output_df(df, replacement_list = var_info)
    return(tibble::as_tibble(df))
  }

  if (output == "stanfit") {
    if (output_needs_rewriting) {
      message("The stanfit object has variable names as x,y,z instead of the ones on the polynomial.
              Don't use pre-compiled model if you need the output to have same names as your polynomial")
    }
    return(samps)
  }
}


create_stan_code <- function(poly, sd, n_eqs, w, homo, vars) {
  d <- get("deriv.mpoly", asNamespace("mpoly"))

  if (n_eqs == 1L) {
    # Single polynomial provided
    if (!is.mpoly(poly)) poly <- mp(poly)
    if (missing(vars)) vars <- mpoly::vars(poly)
    n_vars <- length(vars)
    reorder.mpoly <- get("reorder.mpoly", asNamespace("mpoly"))
    poly <- reorder.mpoly(poly, varorder = sort(vars))

    g_string <- mpoly_to_stan(poly)
    if (homo) {
      grad <- if (n_vars > 1) deriv(poly, var = mpoly::vars(poly)) else gradient(poly)
      ndg_sq <- if (n_vars >1)  Reduce(`+`, grad ^ 2) else grad^2
      ndg_string <- glue::glue("sqrt({mpoly_to_stan(ndg_sq)})")
    } else {
      ndg_string <- "1"
    }

    # Set variables
    if (missing(w)) {
      parms <- glue::glue("real {vars};")
      parms <- paste(parms, collapse = "\n  ")
    } else if (is.numeric(w) && length(w) == 1L) {
      parms <- glue::glue("real<lower=-{w},upper={w}> {vars};")
      parms <- paste(parms, collapse = "\n  ")
    } else {
      parms <- sapply(seq_along(vars), function(i) {
        if (vars[i] %in% names(w)) {
          glue::glue("real<lower={w[[vars[i]]][1]},upper={w[[vars[i]]][2]}> {vars[i]};")
        } else {
          glue::glue("real {vars[i]};")
        }
      })
      parms <- paste(parms, collapse = "\n  ")
    }

    stan_code <- glue::glue(
      "
data {{
  real<lower=0> si;
}}

parameters {{
  {parms}
}}

transformed parameters {{
  real g = {g_string};
  real ndg = {ndg_string};
}}

model {{
  target += normal_lpdf(0.00 | g/ndg, si);
}}
    ")

  } else {
    # Multiple polynomials provided
    vars <- mpoly::vars(poly)
    n_vars <- length(vars)
    printed_polys <- mpolyList_to_stan(poly)

    # Jacobian setup
    printed_jac <- array("", dim = c(n_eqs, n_vars))
    for (i in seq_len(n_eqs)) {
      for (j in seq_len(n_vars)) {
        printed_jac[i, j] <- if (homo) {
          mpoly_to_stan(d(poly[[i]], vars[j]))
        } else if (i == j) {
          "1"
        } else {
          "0"
        }
      }
    }
    printed_jac <- apply(printed_jac, 1, paste, collapse = ", ")
    printed_jac <- paste("      [", printed_jac, "]", collapse = ", \n")
    printed_jac <- paste0("[\n", printed_jac, "\n    ]")

    # Set variables
    if (missing(w)) {
      parms <- paste(sprintf("real %s;", vars), collapse = "\n  ")
    } else if(is.numeric(w) && length(w) == 1L) {
      parms <- glue::glue("real<lower=-{w},upper={w}> {vars};")
      parms <- paste(parms, collapse = "\n  ")
    }else {
      parms <- sapply(seq_along(vars), function(i) {
        if (vars[i] %in% names(w)) {
          sprintf("real<lower=%s,upper=%s> %s;", w[[vars[i]]][1], w[[vars[i]]][2], vars[i])
        } else {
          sprintf("real %s;", vars[i])
        }
      })
      parms <- paste(parms, collapse = "\n  ")
    }
    data_string <- if(length(sd) == 1) "real<lower=0> si" else paste0("cov_matrix[",n_vars,"] si")
    model_string <- if(length(sd) == 1) "normal_lpdf(" else " multi_normal_lpdf("
    mu_string <- if(length(sd) == 1)"0.00" else (paste0("[",paste(rep(0.00, n_vars), collapse = ","),"]'"))
    gbar_string <- if (n_vars == n_eqs) "J \\ g" else if (n_vars > n_eqs) "J' * ((J*J') \\ g)" else "(J'*J) \\ (J'*g)"

    stan_code <- glue::glue(
      "
data {{
  {data_string};
}}

parameters {{
  {parms}
}}

transformed parameters {{
  vector[{n_eqs}] g = [{printed_polys}]';
  matrix[{n_eqs},{n_vars}] J = {printed_jac};
}}

model {{
  target += {model_string}{mu_string} | {gbar_string}, si);
}}
      ")
  }
  stan_code
}











