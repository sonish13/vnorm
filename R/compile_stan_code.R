#' Compile Stan model for user-defined polynomials
#'
#' This function helps to avoid recompiling for polynomials with the same
#' indeterminate but different coefficients.
#'
#' @param poly An mpoly object.
#' @param custom_stan_code If `TRUE`, a custom model is compiled even if the
#'   general case of the polynomial is already included during package
#'   installation. Defaults to `FALSE`.
#' @param w A named list of box constraints for vectors to be passed to Stan.
#'   See [rvnorm()] examples. Defaults to `FALSE`.
#' @param homo If `TRUE`, the sampling is done from a homoskedastic variety
#'   normal distribution. Defaults to `TRUE`.
#'
#'   This function creates a temporary Stan model for the kind of polynomial the
#'   user is interested in. It also creates a variable `compiled_stan_info` in
#'   the global environment for compiled Stan models that is needed to run
#'   `rvnorm` using user-defined Stan codes.
#'
#' @export
#' @examples
#'
#' # compile a model that looks like b0 + bx6 x^6 + by6 y^6 for later input
#' #
#' p <- mp("x^6 + y^6 - 1") # this is the template polynomial not pre-compiled
#' samps <- rvnorm(1000, p, sd = .05)
#' head(samps)
#' compile_stan_code(p) # allows to change coefficients
#' p <- mp("x^6 + 8 y^6 - 1")
#' rvnorm(1e4, p, .05, user_compiled = TRUE)
#'
#'
compile_stan_code <- function(poly, custom_stan_code = FALSE, w = FALSE, homo = TRUE) {
  if (!(is.mpoly(poly) | is.mpolyList(poly))) stop("`poly` should be an mpoly or mpolyList object.")

  if (!custom_stan_code & is.mpoly(poly)) {
    if (length(mpoly::vars(poly)) < 4 & base::max(mpoly::totaldeg(poly)) < 4) {
      stop("Pre-compiled model for the general case already exists from installation.
            Use custom_stan_code = TRUE to use a custom model anyway.")
    }
  }

  stan_code <- get_custom_stan_code(poly = poly, w = w, homo = homo)
  model_name <- generate_model_name(poly = poly, w = w, homo = homo)
  model_path <- cmdstanr::write_stan_file(stan_code, dir = tempdir())
  new_row <- data.frame("name" = model_name, "path" = model_path)

  if (!exists("compiled_stan_info", envir = .GlobalEnv)) {
    assign("compiled_stan_info", new_row, envir = .GlobalEnv)
    message(sprintf("Created registry; registered '%s'", model_name))
  } else {
    compiled_stan_info <- get("compiled_stan_info", envir = .GlobalEnv)
    existing_idx <- which(compiled_stan_info$name == model_name)

    if (length(existing_idx) == 0L) {
      compiled_stan_info <- rbind(compiled_stan_info, new_row)
      message(sprintf("Registered '%s'", model_name))
    } else {
      old_path <- compiled_stan_info$path[existing_idx[1]]
      compiled_stan_info$path[existing_idx[1]] <- model_path
      if (length(existing_idx) > 1L) {
        compiled_stan_info <- compiled_stan_info[-existing_idx[-1], , drop = FALSE]
      }
      if (!identical(old_path, model_path)) {
        message(sprintf("Refreshed path for '%s'", model_name))
      }
    }

    assign("compiled_stan_info", compiled_stan_info, envir = .GlobalEnv)
  }

  cmdstan_model(model_path)

}

get_custom_stan_code <- function(poly, w = FALSE, homo = TRUE) {
  if(is.mpoly(poly)){
    poly <- canonicalize_mpoly(poly)
    vars <- mpoly::vars(poly)
    num_of_vars <- length(vars)

    # Data block
    var_for_data_block <- mpoly::monomials(poly)
    var_for_data_block <- lapply(var_for_data_block, reorder, varorder = vars)
    var_for_data_block <- lapply(var_for_data_block, coef)
    var_for_data_block <- unlist(var_for_data_block)
    var_for_data_block <- get_listed_coeficients(var_for_data_block)
    var_for_data_block <- names(var_for_data_block)

    data_block <- paste(sapply(var_for_data_block, function(x) paste("  real", x)), collapse = "; ")
    data_block <- paste0(data_block, ";")
    if (w) {
      data_block <- paste0(data_block, "  real w;")
    }
    data_block <- paste0("data {\n  real si;\n", data_block, "\n}\n")

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

    # Transformed Parameter Block
    g <- coef_lift(poly)
    g <- mpoly_to_stan(g)

    # Derivatives
    derivatives <- lapply(vars, helper_for_derivative_for_mpoly_stan_code, poly = poly)

    derivative_names <- sapply(seq_along(vars), function(i) {
      paste0("dg", vars[i])
    })

    ndg <- if (homo) {
      paste0("  real ndg = sqrt(", paste0(derivative_names, "^2", collapse = " + "), ");")
    } else {
      "  real ndg = 1;"
    }

    dg <- sapply(seq_along(vars), function(i) {
      paste(derivative_names[i], paste(derivatives[[i]], collapse = ""), sep = " = ")
    })
    dg <- paste0("  real ", dg)
    dg <- paste(dg, collapse = ";")

    model_block <- paste0(
      "model {\n  real g = ", g, ";\n", dg, ";\n", ndg,
      "\n  target += normal_lpdf(0.00 | g/ndg, si); \n}"
    )

    stan_code <- paste0(data_block, params_block, model_block, sep = "")
    model_name <- sprintf(
      "%s_%s%s.stan",
      g,
      if (homo) "vn" else "hvn",
      if (w) "_w" else ""
    )

  }
  else if (is.mpolyList(poly)){
    poly <- canonicalize_mpolylist(poly)
    poly <- sort_mpolylist_lexicographically(poly)
    n_eqs <- length(poly)
    n_vars <- length(vars(poly))
    vars <- list()
    for (i in seq_along(1:length(poly))) {
      vars[[i]] <- vars(poly[[i]])
    }
    var_for_data_block <- list()
    for(i in seq_along(1:length(poly))){
      var_for_data_block[[i]] <- mpoly::monomials(poly[[i]])
      var_for_data_block[[i]] <- lapply(var_for_data_block[[i]], reorder, varorder = vars[[i]])
      var_for_data_block[[i]] <- lapply(var_for_data_block[[i]], coef)
      var_for_data_block[[i]] <- unlist(var_for_data_block[[i]])
      var_for_data_block[[i]] <- get_listed_coeficients(var_for_data_block[[i]])
      var_for_data_block[[i]] <- names(var_for_data_block[[i]])
      var_for_data_block[[i]] <- paste0(var_for_data_block[[i]],"_" ,i)
    }

    var_for_data_block <- unlist(var_for_data_block)
    data_block <- paste(sapply(var_for_data_block, function(x) paste("  real", x)), collapse = "; ")
    data_block <- paste0(data_block, ";")

    if (w) {
      data_block <- paste0(data_block, "\n  real w;")
    }

    data_block <- paste0("data {\n  real si;\n", data_block, "\n}\n")
    vars_for_params <- unique(unlist(vars))
    if (w) {
      params_block <- paste(sapply(vars_for_params, function(var) {
        paste0("  real<lower=-", "w", ", upper=", "w", "> ", var, ";")
      }), collapse = "\n")
    } else {
      params_block <- paste(sapply(vars_for_params, function(var) {
        paste0("  real ", var, ";")
      }), collapse = "\n")
    }
    params_block <- paste0("\nparameters {\n", params_block, "\n}\n")

    # Jacobain Calculation Stuff
    g <- list()
    g_coef <- list()
    g_terms <- list()
    derivatives_pre <- list()
    derivatives <- list()

    for(i in seq_along(1:length(poly))){
      g[[i]] <- helper_for_coef_lift_for_mpolylist(poly[[i]],i)
      g[[i]] <- mpoly_to_stan(g[[i]])
      # Transformed parameter block
      derivatives_pre[[i]] <- lapply(vars(poly), helper_for_derivative_for_mpolylist_stan_code, poly = poly[[i]], i = i)
    }

    g <- unlist(g)
    g <- paste0("  vector[", length(g), "] g = [", paste(g, collapse = ","), "]';")

    if(homo){
      for(i in seq_along(1:length(poly))){
        derivatives[[i]] <- unlist(derivatives_pre[[i]])
      }
      jac <- paste(
        sapply(derivatives, function(v) paste0("      [", paste(v, collapse = ","), "]")),
        collapse = ",\n"
      )
    }else{
      jac <- array("", dim = c(n_eqs, n_vars))
      for (i in 1:n_eqs) {
        for (j in 1:n_vars) {
          jac[i,j] <- if (i == j) "1" else "0"
        }
      }
      jac <- apply(jac, 1L, paste, collapse = ", ")
      jac <- paste("      [", jac, "]", collapse = ", \n")
    }

    dg <- paste0("  matrix[",n_eqs,"," ,n_vars,"] J = [ \n" , jac,"\n    ];")
    trans_block <- paste0("\ntransformed parameters {\n", g, "\n",dg, "\n}\n")


    # Model block
    gbar_string <- if (n_vars == n_eqs) "J \\ g" else if (n_vars > n_eqs) "J' * ((J*J') \\ g)" else "(J'*J) \\ (J'*g)"
    model_block <- paste0("\nmodel {\n  target += normal_lpdf(0.00 |", gbar_string, ", si);\n}")
    stan_code <- paste0(data_block, params_block, trans_block, model_block)
  }else{
    stop("`poly` should either be an mpoly, or an mpolyList object.")
  }
  stan_code
}
