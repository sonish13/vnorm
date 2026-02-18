make_coefficients_data <- function(
    poly, num_of_vars, deg, basis = c("x", "y", "z")
  ) {
  # Build a complete coefficient list (missing terms filled with zero).
  required_coefs <- mpoly::basis_monomials(basis[seq_len(num_of_vars)], deg)
  required_coefs <- lapply(required_coefs, reorder, varorder = basis)
  required_coefs <- lapply(required_coefs, coef)
  required_coefs <- unlist(required_coefs)
  required_coefs <- get_listed_coeficients(required_coefs)
  required_coefs <- lapply(required_coefs, function(x) 0L)

  available_coef <- get_listed_coeficients(coef(poly))
  required_coefs[names(available_coef)] <- available_coef
  required_coefs
}

get_listed_coeficients <- function(coefs) {
  # Convert monomial names to Stan data names (e.g., x^2 y -> bx2y).
  convert_names <- function(term) {
    term <- gsub("\\s+", "", term)
    if (term == "1") return("b1")
    term <- gsub("\\^", "", term)
    paste0("b", term)
  }
  names(coefs) <- sapply(names(coefs), convert_names)
  as.list(coefs)
}

get_coefficeints_data <- function(poly) {
  # Collect named coefficients for either a single polynomial
  # or polynomial list.
  if (is.mpoly(poly)) {
    data <- get_listed_coeficients(coef(poly))
  } else if (is.mpolyList(poly)) {
    poly <- canonicalize_mpolylist(poly)
    poly <- sort_mpolylist_lexicographically(poly)
    convert_names <- function(term, i) {
      term <- gsub("\\s+", "", term)
      if (term == "1") return(paste0("b1_", i))
      term <- gsub("\\^", "", term)
      paste0("b", term, "_", i)
    }
    coefs <- list()
    for (i in seq_along(poly)) {
      coefs[[i]] <- coef(poly[[i]])
      names(coefs[[i]]) <- sapply(names(coefs[[i]]), convert_names, i = i)
    }
    coefs <- unlist(coefs)
    data <- as.list(coefs)
  }
  data
}

check_and_replace_vars <- function(p) {
  # Map arbitrary variable names to x/y/z for precompiled Stan templates.
  current_vars <- mpoly::vars(p)
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
  # Restore original variable names in output after x/y/z remapping.
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
  # Print polynomial in Stan-compatible syntax.
  p <- get("print.mpoly", asNamespace("mpoly"))
  result <- p(mpoly, stars = TRUE, silent = TRUE, plus_pad = 0L, times_pad = 0L)
  result <- gsub("[*]{2}", "^", result)
  result
}

mpolyList_to_stan <- function(mpolyList) {
  # Print a polynomial list as a Stan vector expression payload.
  p <- get("print.mpolyList", asNamespace("mpoly"))
  result <- p(
    mpolyList,
    silent = TRUE,
    stars = TRUE,
    plus_pad = 0,
    times_pad = 0
  )
  result <- gsub("\\*\\*", "^", result)
  result <- paste(result, collapse = ", ")
  result
}



get_derivative <- function(var, num_of_vars, deg, basis = c("x", "y", "z")) {
  # Construct symbolic derivative in Stan syntax for generated template models.
  d <- get("deriv.mpoly", asNamespace("mpoly"))

  num_coef <- mpoly::basis_monomials(basis[seq_len(num_of_vars)], deg)
  num_coef <- lapply(num_coef, reorder, varorder = basis)
  num_coef <- lapply(num_coef, d, var = var)
  num_coef <- lapply(num_coef, mpoly:::coef.mpoly)
  num_coef <- unlist(num_coef)
  num_coef <- unname(num_coef)

  indeterminates <- mpoly::basis_monomials(basis[seq_len(num_of_vars)], deg)
  indeterminates <- lapply(indeterminates, reorder, varorder = basis)
  indeterminates <- lapply(indeterminates, d, var = var)
  indeterminates <- lapply(indeterminates, function(item) {
    item[[1]][names(item[[1]]) == "coef"] <- 1
    item
  })
  indeterminates <- lapply(indeterminates, mpoly_to_stan)
  indeterminates <- unlist(indeterminates)
  indeterminates <- c(indeterminates)

  sym_coef <- mpoly::basis_monomials(basis[seq_len(num_of_vars)], deg)
  sym_coef <- lapply(sym_coef, reorder, varorder = basis)
  sym_coef <- lapply(sym_coef, mpoly:::coef.mpoly)
  sym_coef <- unlist(sym_coef)
  sym_coef <- get_listed_coeficients(sym_coef)
  sym_coef <- names(sym_coef)

  df_for_der <- data.frame(
    num_coef = num_coef,
    indeterminates = indeterminates,
    sym_coef = sym_coef
  )

  df_for_der <- dplyr::filter(df_for_der, num_coef != 0)
  gsub(
    "1\\*|\\*1",
    "",
    paste0(
      df_for_der$num_coef,
      "*",
      df_for_der$sym_coef,
      "*",
      df_for_der$indeterminates,
      collapse = "+"
    )
  )
}
