make_coeficients_data <- function(poly, num_of_vars ,deg, basis = c("x", "y","z")) {
  required_coefs<- mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], deg) |>
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

rename_output_df <- function(df, replacement_list) {
  names(df) <- sapply(names(df), function(col) {
    if (col %in% names(replacement_list)) {
      replacement_list[[col]]
    } else {
      col
    }
  })
  return(df)
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

make_derivative_for_custom_function <- function(var, poly, basis = c("x","y","z")) {
  d <- get("deriv.mpoly", asNamespace("mpoly"))
  df_for_der <- data.frame(num_coef = mpoly::monomials(poly) |>
                             lapply(reorder,varorder = basis) |>
                             lapply(function(item) {
                               item[[1]][names(item[[1]]) == "coef"] <- 1
                               item
                             }) |>
                             lapply(d, var = var) |>
                             lapply(mpoly:::coef.mpoly) |>
                             unlist() |>
                             unname(),
                           # For indeterminates after derivatives
                           indeterminates = mpoly::monomials(poly) |>
                             lapply(d, var = var) |>
                             lapply(function(item) {
                               item[[1]][names(item[[1]]) == "coef"] <- 1
                               item
                             }) |>
                             lapply(mpoly_to_stan) |>
                             unlist() |>
                             c(),
                           # For symbolic coefficients
                           sym_coef = mpoly::monomials(poly) |>
                             lapply(mpoly:::coef.mpoly) |>
                             unlist() |>
                             get_listed_coeficients() |>
                             names()
  )
  df_for_der <- dplyr::filter(df_for_der,num_coef != 0)
  out <- paste0(df_for_der$num_coef,"*",
                df_for_der$sym_coef, "*",
                df_for_der$indeterminates ,collapse = "+")
  gsub("1\\*|\\*1", "", out)
}

make_derivative <- function(var, num_of_vars, deg , basis = c("x","y","z")) {
  d <- get("deriv.mpoly", asNamespace("mpoly"))
  df_for_der <- data.frame(num_coef = mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], deg) |>
                             lapply(reorder,varorder = basis) |>
                             lapply(d, var = var) |>
                             lapply(mpoly:::coef.mpoly) |>
                             unlist() |>
                             unname(),
                           # For indeterminates after derivatives
                           indeterminates = mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], deg) |>
                             lapply(reorder,varorder = basis) |>
                             lapply(d, var = var) |>
                             lapply(function(item) {
                               item[[1]][names(item[[1]]) == "coef"] <- 1
                               item
                             }) |>
                             lapply(mpoly_to_stan) |>
                             unlist() |>
                             c(),
                           # For symbolic coefficients
                           sym_coef = mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], deg) |>
                             lapply(reorder,varorder = basis) |>
                             lapply(mpoly:::coef.mpoly) |>
                             unlist() |>
                             get_listed_coeficients() |>
                             names()
  )
  df_for_der <- dplyr::filter(df_for_der,num_coef != 0)
  out <- paste0(df_for_der$num_coef,"*",
                df_for_der$sym_coef, "*",
                df_for_der$indeterminates ,collapse = "+")
  gsub("1\\*|\\*1", "", out)
}

get_model_name <- function(poly, w =T, homo=T) {
  g = paste0(
    mpoly::monomials(poly) |> lapply(function(item) {
      item[[1]][names(item[[1]]) == "coef"] <- 1
      item
    }) |>
      lapply(coef) |>
      unlist() |>
      get_listed_coeficients() |>
      names(),
    "*" ,
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
  sprintf("%s_%s%s.stan",
          g,
          if (homo)
            "vn"
          else
            "hvn",
          if (w)
            "_w"
          else
            "")
}

