test_that("Compiles model for a simple polynomial", {
  poly <- mpoly::mp("x^2 + y^2 - 1")
  result <- compile_stan_code(poly, custom_stan_code = TRUE, w = FALSE, homo = TRUE)
  expect_equal(result, "Model Compiled")
})

test_that("Stops if pre-compiled model already exists", {
  poly <- mp("x^2 + y^2")
  expect_error(
    compile_stan_code(poly, custom_stan_code = FALSE, w = FALSE, homo = TRUE),
    "Pre-compiled model for the general case already exists from installation.
            Use custom_stan_code = TRUE to use a custom model anyway."
  )
})

test_that("Errors when compiling for mpolyList object", {
  poly <- mpoly::mp(c("x^2 + y^2 - 1","y - x"))
  expect_error(
    compile_stan_code(poly, custom_stan_code = TRUE, w = FALSE, homo = TRUE),
    "Cannot compile model for an mpolyList object"
  )
})

test_that("Creates or updates global environment variable `compiled_stan_info`", {
  poly <- mp("x^2 + y^2 - 1")
  compile_stan_code(poly, custom_stan_code = TRUE, w = FALSE, homo = TRUE)
  expect_true(exists("compiled_stan_info", envir = .GlobalEnv))
  compiled_info <- get("compiled_stan_info", envir = .GlobalEnv)
  expect_true("name" %in% names(compiled_info))
  expect_true("path" %in% names(compiled_info))
})

test_that("Compiles model with box constraints (w)", {
  poly <- mp("x^2 + y^2 - 1")
  result <- compile_stan_code(poly, custom_stan_code = TRUE, w = 1, homo = TRUE)
  expect_equal(result, "Model Compiled")
})

test_that("Compiles model with homoskedastic option", {
  poly <- mp("x^2 + y^2 - 1")
  result <- compile_stan_code(poly, custom_stan_code = TRUE, w = FALSE, homo = TRUE)
  expect_equal(result, "Model Compiled")
})

