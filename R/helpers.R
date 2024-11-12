make_coeficients_data <- function(poly, num_of_vars, deg, basis = c("x", "y", "z")) {
  required_coefs <- mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], deg)
  required_coefs <- lapply(required_coefs, reorder, varorder = basis)
  required_coefs <- lapply(required_coefs, coef)
  required_coefs <- unlist(required_coefs)
  required_coefs <- get_listed_coeficients(required_coefs)
  required_coefs <- lapply(required_coefs, function(x) 0)

  available_coef <- get_listed_coeficients(coef(poly))
  required_coefs[names(available_coef)] <- available_coef
  required_coefs
}

get_listed_coeficients <- function(coefs) {
  convert_names <- function(term) {
    term <- gsub("\\s+", "", term)  # Remove spaces
    if (term == "1") return("b0")   # Constant term should be "b0"
    term <- gsub("\\^", "", term)   # Remove power symbol (^)
    paste0("b", term)               # Add "b" at the beginning
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
  temp_vars <- paste0("tmp", seq_along(current_vars))  # Temporary placeholders

  for (i in seq_along(current_vars)) {
    p <- swap(p, current_vars[i], temp_vars[i])
    var_mapping[[temp_vars[i]]] <- current_vars[i]
  }

  for (i in seq_along(expected_vars)) {
    p <- swap(p, temp_vars[i], expected_vars[i])
    var_mapping[[expected_vars[i]]] <- var_mapping[[temp_vars[i]]]
    var_mapping[[temp_vars[i]]] <- NULL
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
  df
}

mpoly_to_stan <- function(mpoly) {
  p <- get("print.mpoly", asNamespace("mpoly"))
  result <- p(mpoly, stars = TRUE, silent = TRUE, plus_pad = 0L, times_pad = 0L)
  result <- gsub("[*]{2}", "^", result)
  result
}

mpolyList_to_stan <- function(mpolyList) {
  p <- get("print.mpolyList", asNamespace("mpoly"))
  result <- p(mpolyList, silent = TRUE, stars = TRUE, plus_pad = 0, times_pad = 0)
  result <- gsub("\\*\\*", "^", result)
  result <- paste(result, collapse = ", ")
  result
}



get_derivative <- function(var, num_of_vars, deg, basis = c("x", "y", "z")) {
  d <- get("deriv.mpoly", asNamespace("mpoly"))

  # Calculate num_coef
  num_coef <- mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], deg)
  num_coef <- lapply(num_coef, reorder, varorder = basis)
  num_coef <- lapply(num_coef, d, var = var)
  num_coef <- lapply(num_coef, mpoly:::coef.mpoly)
  num_coef <- unlist(num_coef)
  num_coef <- unname(num_coef)

  # Calculate indeterminates
  indeterminates <- mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], deg)
  indeterminates <- lapply(indeterminates, reorder, varorder = basis)
  indeterminates <- lapply(indeterminates, d, var = var)
  indeterminates <- lapply(indeterminates, function(item) {
    item[[1]][names(item[[1]]) == "coef"] <- 1
    item
  })
  indeterminates <- lapply(indeterminates, mpoly_to_stan)
  indeterminates <- unlist(indeterminates)
  indeterminates <- c(indeterminates)

  # Calculate sym_coef
  sym_coef <- mpoly::basis_monomials(basis[seq_along(1:num_of_vars)], deg)
  sym_coef <- lapply(sym_coef, reorder, varorder = basis)
  sym_coef <- lapply(sym_coef, mpoly:::coef.mpoly)
  sym_coef <- unlist(sym_coef)
  sym_coef <- get_listed_coeficients(sym_coef)
  sym_coef <- names(sym_coef)

  # Create data frame
  df_for_der <- data.frame(
    num_coef = num_coef,
    indeterminates = indeterminates,
    sym_coef = sym_coef
  )

  # Filter and format output
  df_for_der <- dplyr::filter(df_for_der, num_coef != 0)
  out <- paste0(df_for_der$num_coef, "*", df_for_der$sym_coef, "*", df_for_der$indeterminates, collapse = "+")
  gsub("1\\*|\\*1", "", out)
}


get_model_name <- function(poly, w = TRUE, homo = TRUE) {
  monomials <- mpoly::monomials(poly)

  coef_names <- lapply(monomials, function(item) {
    item[[1]][names(item[[1]]) == "coef"] <- 1
    item
  })
  coef_names <- lapply(coef_names, coef)
  coef_names <- unlist(coef_names)
  coef_names <- get_listed_coeficients(coef_names)
  coef_names <- names(coef_names)

  indeterminates <- lapply(monomials, function(item) {
    item[[1]][names(item[[1]]) == "coef"] <- 1
    item
  })
  indeterminates <- lapply(indeterminates, mpoly_to_stan)
  indeterminates <- unlist(indeterminates)
  indeterminates <- c(indeterminates)

  g <- paste0(coef_names, "*", indeterminates, collapse = "+")
  g <- gsub("1\\*|\\*1", "", g)

  sprintf(
    "%s_%s%s.stan",
    g,
    if (homo) "vn" else "hvn",
    if (w) "_w" else ""
  )
}

