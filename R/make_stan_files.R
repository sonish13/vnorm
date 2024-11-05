# This function creates stan files we need for pre-compiling stan models.
# This is for univariety models

make_stan_files <- function(num_of_vars,
                           totaldeg,
                           homo = TRUE,
                           w = TRUE,
                           basis = c("x", "y", "z")) {

  vars <- basis[seq_along(1:num_of_vars)]

  # Data block
  # Get the necesasary basis monomials for the required number of variables and
  # their total degrees and create the names of the coefficients of the general
  # polynomial
  var_for_data_block <- mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], totaldeg) |>
    lapply(reorder, varorder = basis) |>
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
  # Parameter block to include the indeterminates and use the window w if required
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
  # make the polynomail for stan codes

  # This creates the coefficients
  g <- paste0(
    mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], totaldeg) |>
      lapply(reorder, varorder = basis) |>
      lapply(coef) |>
      unlist() |>
      get_listed_coeficients() |>
      names(),
    "*" ,
    # This creates the terms
    mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], totaldeg) |>
      lapply(reorder, varorder = basis) |>
      lapply(mpoly_to_stan) |>
      unlist() |>
      c(),
    collapse = "+"
  )

  # Remove unnecessaru 1's
  g <- gsub("1\\*|\\*1", "", g)
  # Get the derivatives w.r.t all the variables
  derivatives <- vars |> lapply(make_derivative, num_of_vars = num_of_vars, deg = totaldeg, basis = vars)

  # names of the derivatives in stan files, dgx for variable x, dgy for variable y and so on
  derivative_names <- sapply(seq_along(vars), function(i) {
    paste0("dg", vars[i])
  })

  # make the derivate values and add it with derivative names
  dg <- sapply(seq_along(vars), function(i)
    paste(derivative_names[i], paste(derivatives[[i]], collapse = ""), sep = " = "))

  dg <- paste0("real ", dg) |> paste0(collapse = ";")

  if (homo) {
    ndg <- paste0("real ndg = sqrt(",
                  paste0(derivative_names, "^2", collapse = " + "),
                  ");")
  } else {
    ndg <- "real ndg = 1;"
  }


  # create the complete model block
  model_block <- paste0(
    "model {\nreal g = ", g, ";\n", dg, ";\n", ndg,
    "\ntarget += normal_lpdf(0.00 | g/ndg, si); \n}"
  )

  # complete stan code
  stan_code <- paste0(data_block, params_block, model_block, sep = "")

  filename <- sprintf(
    "%i_%i_%s%s.stan",
    length(vars),
    totaldeg,
    if (homo) "vn" else "hvn",
    if (w) "_w" else ""
  )

  # set path and write files
  path <- here::here("src", "stan", filename)
  readr::write_lines(stan_code, file = path)
}
