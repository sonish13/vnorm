#' Compile Stan model for user-defined polynomials
#'
#' This function helps to avoid recompiling for polynomials
#' with the same indeterminate but different coefficients.
#'
#' @param poly An mpoly object.
#' @param custom_stan_code If `TRUE`, a custom model is compiled even if
#'   the general case of the polynomial is already included during package installation.
#'   Defaults to `FALSE`.
#' @param w A named list of box constraints for vectors to be passed to Stan. See
#'   [rvnorm()] examples. Defaults to `FALSE`.
#' @param homo If `TRUE`, the sampling is done from a homoskedastic variety
#'   normal distribution. Defaults to `TRUE`.
#'
#' This function creates a temporary Stan model for the kind of polynomial the user is
#' interested in. It also creates a variable `compiled_stan_info` in the global
#' environment for compiled Stan models that is needed to run `rvnorm` using user-defined Stan codes.
#'
#' @export
compile_stan_code <- function(poly, custom_stan_code = FALSE, w = FALSE, homo = TRUE) {
  vars <- mpoly::vars(poly)
  num_of_vars <- length(vars)
  deg <- mpoly::totaldeg(poly)

  if (is.mpoly(poly)) {
    n_eqs <- 1L
  } else if (is.mpolyList(poly)) {
    n_eqs <- length(poly)
  }

  if (n_eqs > 1) {
    stop("Cannot compile model for an mpolyList object")
  }

  if (!custom_stan_code) {
    if (length(mpoly::vars(poly)) < 4 & base::max(mpoly::totaldeg(poly) < 4)) {
      stop("Pre-compiled model for the general case already exists from installation.
            Use custom_stan_code = TRUE to use a custom model anyway.")
    }
  }

  stan_code <- make_custom_stan_code(poly = poly, w = w, homo = homo)
  model_name <- get_model_name(poly = poly, w = w, homo = homo)
  model_path <- cmdstanr::write_stan_file(stan_code, getwd())
  compiled_stan_info <- data.frame("name" = model_name, "path" = model_path)

  if (!exists("compiled_stan_info", envir = .GlobalEnv)) {
    assign("compiled_stan_info", compiled_stan_info, envir = .GlobalEnv)
    message("compiled_stan_info variable created in global environment")
  } else {
    compiled_stan_info <- get("compiled_stan_info", envir = .GlobalEnv)
    new_row <- data.frame("name" = model_name, "path" = model_path)
    compiled_stan_info <- rbind(compiled_stan_info, new_row)
    compiled_stan_info <- dplyr::distinct(compiled_stan_info)
    assign("compiled_stan_info", compiled_stan_info, envir = .GlobalEnv)
    message("compiled_stan_info variable updated in global environment")
  }

  cmdstan_model(model_path)
  return("Model Compiled")
}

make_custom_stan_code <- function(poly, w = FALSE, homo = TRUE) {
  vars <- mpoly::vars(poly)
  num_of_vars <- length(vars)

  # Data block
  var_for_data_block <- mpoly::monomials(poly) |>
    lapply(reorder, varorder = vars) |>
    lapply(coef) |>
    unlist() |>
    get_listed_coeficients() |>
    names()

  data_block <- paste(sapply(var_for_data_block, function(x) paste("real", x)), collapse = "; ") |>
    paste0(";")
  if (w) {
    data_block <- paste0(data_block, "real w;")
  }
  data_block <- paste0("data {\nreal si;\n", data_block, "\n}\n")

  # Parameter block
  if (w) {
    params_block <- paste(sapply(vars, function(var) {
      paste0("  real<lower=-", "w", ", upper=", "w", "> ", var, ";")
    }), collapse = "\n")
  } else {
    params_block <- paste(sapply(vars, function(var) {
      paste0("  real ", var, ";")
    }), collapse = "\n")
  }
  params_block <- paste0("parameters {\n", params_block, "\n }\n")

  # Model block
  g <- paste0(
    mpoly::monomials(poly) |> lapply(function(item) {
      item[[1]][names(item[[1]]) == "coef"] <- 1
      item
    }) |>
      lapply(coef) |>
      unlist() |>
      get_listed_coeficients() |>
      names(),
    "*",
    mpoly::monomials(poly) |> lapply(function(item) {
      item[[1]][names(item[[1]]) == "coef"] <- 1
      item
    }) |>
      lapply(mpoly_to_stan) |>
      unlist() |>
      c(),
    collapse = "+"
  )

  g <- gsub("1\\*|\\*1", "", g)
  derivatives <- vars |> lapply(make_derivative_for_custom_function, poly = poly, basis = vars)

  derivative_names <- sapply(seq_along(vars), function(i) {
    paste0("dg", vars[i])
  })

  ndg <- if (homo) {
    paste0("real ndg = sqrt(", paste0(derivative_names, "^2", collapse = " + "), ");")
  } else {
    "real ndg = 1;"
  }

  dg <- sapply(seq_along(vars), function(i)
    paste(derivative_names[i], paste(derivatives[[i]], collapse = ""), sep = " = "))
  dg <- paste0("real ", dg) |> paste0(collapse = ";")

  model_block <- paste0(
    "model {\nreal g = ", g, ";\n", dg, ";\n", ndg,
    "\ntarget += normal_lpdf(0.00 | g/ndg, si); \n}"
  )

  stan_code <- paste0(data_block, params_block, model_block, sep = "")
  model_name <- sprintf(
    "%s_%s%s.stan",
    g,
    if (homo) "vn" else "hvn",
    if (w) "_w" else ""
  )

  stan_code
}

make_derivative_for_custom_function <- function(var, poly, basis = c("x", "y", "z")) {
  d <- get("deriv.mpoly", asNamespace("mpoly"))

  df_for_der <- data.frame(
    num_coef = mpoly::monomials(poly) |>
      lapply(reorder, varorder = basis) |>
      lapply(function(item) {
        item[[1]][names(item[[1]]) == "coef"] <- 1
        item
      }) |>
      lapply(d, var = var) |>
      lapply(mpoly:::coef.mpoly) |>
      unlist() |>
      unname(),

    indeterminates = mpoly::monomials(poly) |>
      lapply(d, var = var) |>
      lapply(function(item) {
        item[[1]][names(item[[1]]) == "coef"] <- 1
        item
      }) |>
      lapply(mpoly_to_stan) |>
      unlist() |>
      c(),

    sym_coef = mpoly::monomials(poly) |>
      lapply(mpoly:::coef.mpoly) |>
      unlist() |>
      get_listed_coeficients() |>
      names()
  )

  df_for_der <- dplyr::filter(df_for_der, num_coef != 0)

  out <- paste0(df_for_der$num_coef, "*", df_for_der$sym_coef, "*", df_for_der$indeterminates, collapse = "+")
  gsub("1\\*|\\*1", "", out)
}
