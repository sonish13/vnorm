test_that("Compiles model correctly for an mpoly object", {
  p <- mp("x^2 + y^2 - 1")
  result <- compile_stan_code(p, custom_stan_code = TRUE, w = FALSE, homo = TRUE) |>
    suppressMessages()
  expect_equal(
    result$code(),
    c(
      "data {",
      "  real si;",
      "  real bx2;   real by2;   real b1;",
      "}",
      "parameters {",
      "  real x;",
      "  real y;",
      " }",
      "model {",
      "  real g = bx2*x^2+by2*y^2+b1;",
      "  real dgx = 2*bx2*x;  real dgy = 2*by2*y;",
      "  real ndg = sqrt(dgx^2 + dgy^2);",
      "  target += normal_lpdf(0.00 | g/ndg, si); ",
      "}")
  )
})

test_that("Stops if pre-compiled model already exists", {
  p <- mp("x^2 + y^2")
  expect_error(
    compile_stan_code(p, custom_stan_code = FALSE, w = FALSE, homo = TRUE),
    "Pre-compiled model for the general case already exists from installation.
            Use custom_stan_code = TRUE to use a custom model anyway."
  )
})

test_that("Compiles model correctly for an mpolyList object", {
  p <- mp(c("x^2 + y^2 - 1","y - x"))
  result <- compile_stan_code(p, custom_stan_code = TRUE, w = FALSE, homo = TRUE) |>
    suppressMessages()
  expect_equal(
    result$code(),
    c(
      "data {",
      "  real si;",
      "  real bx_1;   real by_1;   real bx2_2;   real by2_2;   real b1_2;",
      "}",
      "",
      "parameters {",
      "  real x;",
      "  real y;",
      "}",
      "",
      "transformed parameters {",
      "  vector[2] g = [bx_1*x+by_1*y,bx2_2*x^2+by2_2*y^2+b1_2]';",
      "  matrix[2,2] J = [ ",
      "      [bx_1,by_1],",
      "      [2*bx2_2*x,2*by2_2*y]",
      "    ];",
      "}",
      "",
      "model {",
      "  target += normal_lpdf(0.00 |J \\ g, si);",
      "}"
    )
  )
})

test_that("Creates global environment variable `compiled_stan_info`", {
  p <- mp("x^2 + y^2 - 1")
  compile_stan_code(p, custom_stan_code = TRUE, w = FALSE, homo = TRUE) |>
    suppressMessages()

  # check compiled_stan_info is created and the right kind of object
  expect_true( exists("compiled_stan_info", envir = .GlobalEnv) )
  compiled_info <- get("compiled_stan_info", envir = .GlobalEnv)
  expect_true( is.data.frame(compiled_info) )
  expect_equal( names(compiled_info), c("name", "path") )
})

test_that("Compiles model with box constraints w", {
  p <- mp("x^2 + y^2 - 1")
  result <- compile_stan_code(p, custom_stan_code = TRUE, w = 1, homo = TRUE)
  expect_equal(
    result$code(),
    c(
      "data {",
      "  real si;",
      "  real bx2;   real by2;   real b1;  real w;",
      "}",
      "parameters {",
      "  real<lower=-w, upper=w> x;",
      "  real<lower=-w, upper=w> y;",
      " }",
      "model {",
      "  real g = bx2*x^2+by2*y^2+b1;",
      "  real dgx = 2*bx2*x;  real dgy = 2*by2*y;",
      "  real ndg = sqrt(dgx^2 + dgy^2);",
      "  target += normal_lpdf(0.00 | g/ndg, si); ",
      "}"
    )
  )
})

test_that("Compiles model with homoskedastic option", {
  p <- mp("x^2 + y^2 - 1")
  result <- compile_stan_code(p, custom_stan_code = TRUE, w = FALSE, homo = TRUE)
  expect_equal(
    result$code(),
    c(
      "data {",
      "  real si;",
      "  real bx2;   real by2;   real b1;",
      "}",
      "parameters {",
      "  real x;", "  real y;", " }",
      "model {",
      "  real g = bx2*x^2+by2*y^2+b1;",
      "  real dgx = 2*bx2*x;  real dgy = 2*by2*y;",
      "  real ndg = sqrt(dgx^2 + dgy^2);",
      "  target += normal_lpdf(0.00 | g/ndg, si); ",
      "}"
    )
  )
})

