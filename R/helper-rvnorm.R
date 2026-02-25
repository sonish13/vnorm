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
