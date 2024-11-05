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
#' @param homo If \code{TRUE}, the sampling is done from homoskedastic variety
#'  normal distribution.
#' @param normalized If \code{TRUE}, the polynomial is gradient-normalized. This
#'   is highly recommended. Set to \code{TRUE} if homo = \code{TRUE}
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
#' @param user_compiled If\code{TRUE}, user compiled stan program made using
#' [compile_stan_code] is used. Defaults to \code{FALSE}
#' @name rvnorm
#' @return Either (1) matrix whose rows are the individual draws from the
#'   distribution, (2) a [tbl_df-class] object with the draws along with
#'   additional information, or (3) an object of class [stanfit-class].
#' @examples
#' # Example code provided in documentation
#'

#' @rdname rvnorm
#' @export
rvnorm <- function(n, poly, sd, output = "simple", chains = 4L, warmup = max(500, floor(n / 2)),
                   inc_warmup = FALSE, thin = 1L, inject_direct = FALSE, verbose = FALSE,
                   cores = min(chains, getOption("mc.cores", 1L)), homo = TRUE,
                   normalized = homo, w, vars, numerator, denominator, refresh = 0L,
                   code_only = FALSE, pre_compiled = TRUE, user_compiled = FALSE, ...) {

  # Initialization and checks
  output_needs_rewriting <- FALSE
  if (is.character(poly)) poly <- mp(poly)
  if (is.mpoly(poly)) {
    n_eqs <- 1L
  } else if (is.mpolyList(poly)) {
    n_eqs <- length(poly)
  } else {
    stop("`poly` should be either a character vector, mpoly, or mpolyList.", call. = FALSE)
  }
  deg <- mpoly::totaldeg(poly)
  num_of_vars <- length(mpoly::vars(poly))
  if (refresh) refresh <- if (verbose) max(ceiling(n / 10), 1L) else 0L
  if (n_eqs > 1) pre_compiled <- FALSE
  if (n_eqs > 1 || length(mpoly::vars(poly)) > 3 || base::max(mpoly::totaldeg(poly)) > 3)
    pre_compiled <- FALSE
  if (code_only) return(create_stan_code(poly, sd, n_eqs, w, normalized))

  # Model selection and data preparation
  if (user_compiled) { # incase we want use model that user compiled
    model_name <- get_model_name(poly = poly, homo = homo, w = !missing(w))
    model_path <- compiled_stan_info$path[compiled_stan_info$name == model_name]
    model <- cmdstan_model(model_path)
    stan_data <- get_listed_coeficients(coef(poly))
    stan_data <- if (missing(w)) c(stan_data, "si" = sd) else c(stan_data, "w" = w, "si" = sd)
  } else if (pre_compiled) { # data and model if/when pre-compiled model is used
    #replace variable if poly has anything other then x, y or z
    list_for_transformation <- check_and_replace_vars(poly)
    if (!is.null(unlist(list_for_transformation$mapping))) {
      poly <- list_for_transformation$polynomial
      var_info <- list_for_transformation$mapping
      output_needs_rewriting <- TRUE
    }
    # det name of stan model to use
    stan_file_name <- paste(num_of_vars, deg, if (normalized | homo) "hvn" else "vn", sep = "_")
    if (!missing(w)) stan_file_name <- paste(stan_file_name, "w", sep = "_")
    model <- instantiate::stan_package_model(name = stan_file_name, package = "vnorm")
    stan_data <- make_coeficients_data(poly, num_of_vars, deg)
    stan_data <- if (missing(w)) c(stan_data, "si" = sd) else c(stan_data, "w" = w, "si" = sd)
  } else {
    # create code to run stan model from scratch
    stan_code <- create_stan_code(poly, sd, n_eqs, w, normalized, vars)
    stan_file <- write_stan_file(stan_code)
    model <- cmdstan_model(stan_file)
    stan_data <- list("si" = sd)
  }

  # Run Stan sampler
  samps <- model$sample(data = stan_data, refresh = refresh, iter_warmup = warmup,
                        iter_sampling = ceiling(n / chains), chains = chains,
                        parallel_chains = cores, adapt_delta = .999, max_treedepth = 20L)

  # Parse output and return
  if (output == "simple") {
    df <- samps$draws(format = "df", inc_warmup = inc_warmup) |> as.data.frame()
    df <- df[(nrow(df) - n + 1):nrow(df), mpoly::vars(poly)]
    row.names(df) <- NULL
    if (output_needs_rewriting) df <- rename_output_df(df, replacement_list = var_info)
    return(df)
  }

  if (output == "tibble") {
    df <- samps$draws(format = "df", inc_warmup = inc_warmup) |> tibble::as_tibble()
    df <- df[(nrow(df) - n + 1):nrow(df), ]
    row.names(df) <- NULL
    df <- cbind(df[, -1], df[, 1])  # Move lp__ to last column
    if (output_needs_rewriting) df <- rename_output_df(df, replacement_list = var_info)
    return(as_tibble(df))
  }

  if (output == "stanfit") {
    if (output_needs_rewriting) {
      message("The stanfit object has variable names as x,y,z instead of the ones on the polynomial.
              Don't use pre-compiled model if you need the output to have same names as your polynomial")
    }
    return(samps)
  }
}

# Function to create stan code
create_stan_code <- function(poly, sd, n_eqs, w, normalized, vars) {
  d <- get("deriv.mpoly", asNamespace("mpoly"))

  if (n_eqs == 1L) {
    # Single polynomial provided
    if (!is.mpoly(poly)) poly <- mp(poly)
    if (missing(vars)) vars <- mpoly::vars(poly)
    n_vars <- length(vars)
    reorder.mpoly <- get("reorder.mpoly", asNamespace("mpoly"))
    poly <- reorder.mpoly(poly, varorder = sort(vars))

    g_string <- mpoly_to_stan(poly)
    if (normalized) {
      grad <- if (n_vars > 1) deriv(poly, var = mpoly::vars(poly)) else gradient(poly) ^ 2
      ndg_sq <- Reduce(`+`, grad ^ 2)
      ndg_string <- glue::glue("sqrt({mpoly_to_stan(ndg_sq)})")
    } else {
      ndg_string <- "1"
    }

    # Set variables
    parms <- if (missing(w)) {
      glue::glue("real {vars};")
    } else if (is.numeric(w) && length(w) == 1L) {
      glue::glue("real<lower=-{w},upper={w}> {vars};")
    } else {
      parms <- sapply(seq_along(vars), function(i) {
        if (vars[i] %in% names(w)) {
          glue::glue("real<lower={w[[vars[i]]][1]},upper={w[[vars[i]]][2]}> {vars[i]};")
        } else {
          glue::glue("real {vars[i]};")
        }
      }) |> str_c(collapse = "\n  ")
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
        printed_jac[i, j] <- if (normalized) {
          mpoly_to_stan(d(poly[[i]], vars[j]))
        } else if (i == j) {
          "1"
        } else {
          "0"
        }
      }
    }
    printed_jac <- apply(printed_jac, 1, str_c, collapse = ", ") |>
      str_c("      [", ., "]", collapse = ", \n") |>
      str_c("[\n", ., "\n    ]")

    # Set variables
    parms <- if (missing(w)) {
      glue::glue("real {vars};")
    } else {
      sapply(seq_along(vars), function(i) {
        if (vars[i] %in% names(w)) {
          glue::glue("real<lower={w[[vars[i]]][1]},upper={w[[vars[i]]][2]}> {vars[i]};")
        } else {
          glue::glue("real {vars[i]};")
        }
      }) |> str_c(collapse = "\n  ")
    }

    gbar_string <- if (n_vars == n_eqs) "J \\ g" else if (n_vars > n_eqs) "J' * ((J*J') \\ g)" else "(J'*J) \\ (J'*g)"

    stan_code <- glue::glue(
      "
data {{
  real<lower=0> si;
}}

parameters {{
  {parms}
}}

transformed parameters {{
  vector[{n_eqs}] g = [{printed_polys}]';
  matrix[{n_eqs},{n_vars}] J = {printed_jac};
}}

model {{
  target += normal_lpdf(0.00 | {gbar_string}, si);
}}
      ")
  }
  stan_code
}
