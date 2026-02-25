test_that("create_stan_code builds single-polynomial Stan code", {
  p <- mp("x^2 + y^2 - 1")

  code_homo <- create_stan_code(p, sd = 0.1, n_eqs = 1L, homo = TRUE)
  code_hetero_w <- create_stan_code(p, sd = 0.1, n_eqs = 1L, w = 5, homo = FALSE)

  expect_type(code_homo, "character")
  expect_length(code_homo, 1)
  expect_match(code_homo, "normal_lpdf\\(")
  expect_match(code_homo, "real g =")
  expect_match(code_homo, "real ndg = sqrt\\(")

  expect_match(code_hetero_w, "real<lower=-5,upper=5> x;")
  expect_match(code_hetero_w, "real<lower=-5,upper=5> y;")
  expect_match(code_hetero_w, "real ndg = 1;")
})

test_that("create_stan_code supports named bounds for single-polynomial case", {
  p <- mp("x^2 + y^2 - 1")
  w <- list(x = c(-2, 1), y = c(-3, 3))

  code <- create_stan_code(p, sd = 0.1, n_eqs = 1L, w = w, homo = TRUE)

  expect_match(code, "real<lower=-2,upper=1> x;")
  expect_match(code, "real<lower=-3,upper=3> y;")
})

test_that("create_stan_code builds multi-polynomial Stan code", {
  p <- mp(c("x^2 + y^2 - 1", "y"))

  code_scalar <- create_stan_code(p, sd = 0.1, n_eqs = 2L, homo = TRUE)
  code_matrix <- create_stan_code(p, sd = diag(c(1, 2)), n_eqs = 2L, homo = FALSE)

  expect_match(code_scalar, "vector\\[2\\] g =")
  expect_match(code_scalar, "matrix\\[2,2\\] J =")
  expect_match(code_scalar, "normal_lpdf\\(")

  expect_match(code_matrix, "cov_matrix\\[2\\] si")
  expect_match(code_matrix, "multi_normal_lpdf\\(")
  expect_match(code_matrix, "matrix\\[2,2\\] J =")
})
