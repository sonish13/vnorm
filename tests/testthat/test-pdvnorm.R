test_that("pdvnorm matches normal density for p(x)=x", {
  p <- mp("x")
  x <- c(-1, 0, 1)

  dens <- pdvnorm(x, p, sigma = 1)
  log_dens <- pdvnorm(x, p, sigma = 1, log = TRUE)

  expect_equal(dens, dnorm(x))
  expect_equal(log_dens, dnorm(x, log = TRUE))
})

test_that("pdvnorm supports matrix and data.frame input for mpolyList", {
  p <- mp(c("x", "y"))
  x <- rbind(c(0, 0), c(1, 2), c(-1, 1))

  out_homo <- pdvnorm(x, p, sigma = c(1, 2), homo = TRUE)
  out_hetero <- pdvnorm(as.data.frame(x), p, sigma = diag(c(1, 4)), homo = FALSE)

  expect_type(out_homo, "double")
  expect_length(out_homo, nrow(x))
  expect_true(all(is.finite(out_homo)))

  expect_type(out_hetero, "double")
  expect_length(out_hetero, nrow(x))
  expect_true(all(is.finite(out_hetero)))
})

test_that("pdvnorm validates sigma shape and positivity", {
  p <- mp(c("x", "y"))
  x <- c(0, 0)

  expect_error(
    pdvnorm(x, p, sigma = c(1, 2, 3), homo = FALSE),
    "length\\(sigma\\) must be m"
  )
  expect_error(
    pdvnorm(x, p, sigma = c(1, -1), homo = TRUE),
    "vector entries must be positive"
  )
  expect_error(
    pdvnorm(x, p, sigma = matrix(c(1, 2, 2, 4), 2, 2), homo = TRUE),
    "positive definite"
  )
})
