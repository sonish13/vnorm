make_stan_file <- function(num_of_vars,
                           totaldeg,
                           homo = TRUE,
                           w = TRUE,
                           basis = c("x", "y", "z")) {
  vars <- basis[seq_along(1:num_of_vars)]
  # Data block
  var_for_data_block <- mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], totaldeg) |>
    lapply(reorder, varorder = basis) |>
    lapply(coef) |>
    unlist() |>
    get_listed_coeficients() |>
    names()
  data_block <- paste(sapply(var_for_data_block, function(x)
    paste("real", x)), collapse = "; ") |>
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
  g = paste0(
    mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], totaldeg) |>
      lapply(reorder, varorder = basis) |>
      lapply(coef) |>
      unlist() |>
      get_listed_coeficients() |>
      names(),
    "*" ,
    mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], totaldeg) |>
      lapply(reorder, varorder = basis) |>
      lapply(mpoly_to_stan) |>
      unlist() |>
      c(),
    collapse = "+"
  )
  g <- gsub("1\\*|\\*1", "", g)
  derivatives <- vars |> lapply(make_derivative, num_of_vars = num_of_vars, deg = totaldeg, basis = vars)

  derivative_names <- sapply(seq_along(vars), function(i) {
    paste0("dg", vars[i])
  })
  if (homo) {
    ndg = paste0("real ndg = sqrt(",
                 paste0(derivative_names, "^2", collapse = " + "),
                 ");")
  } else{
    ndg = "real ndg =1;"
  }

  dg <- sapply(seq_along(vars), function(i)
    paste(derivative_names[i], paste(derivatives[[i]], collapse = ""), sep = " = "))
  dg <- paste0("real ", dg) |> paste0(collapse = ";")
  model_block <- paste0("model {\nreal g = ",
                        g,
                        ";\n",
                        dg,
                        ";\n",
                        ndg,
                        "\ntarget += normal_lpdf(0.00 | g/ndg, si); \n}")
  stan_code <- paste0(data_block, params_block, model_block, sep = "")
  filename <- sprintf("%i_%i_%s%s.stan",
                      length(vars),
                      totaldeg,
                      if (homo)
                        "vn"
                      else
                        "hvn",
                      if (w)
                        "_w"
                      else
                        "")
  path <- here::here("src", "stan", filename)
  readr::write_lines(stan_code, file = path)
  #stan_code

}


