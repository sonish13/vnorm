fake_cmdstan_fit <- function(draws_df) {
  list(
    draws = function(format = "df", inc_warmup = FALSE) {
      draws_df
    }
  )
}

test_that("variety_solve validates input class and dimensionality", {
  expect_error(
    variety_solve(123),
    "`polylist` must be an mpoly or mpolyList object."
  )

  expect_error(
    variety_solve(mp(c("x", "y", "x + y")), n = 10),
    "not zero-dimensional"
  )
})

test_that("variety_solve returns rounded posterior means using rvnorm samples", {
  p <- mp("x")
  fake_df <- data.frame(x = c(0.11, 0.09, 0.10), lp__ = c(0, 0, 0))

  testthat::local_mocked_bindings(
    rvnorm = function(...) fake_cmdstan_fit(fake_df),
    .package = "vnorm"
  )

  out <- variety_solve(p, n = 3, sig_digit = 2)

  expect_type(out, "double")
  expect_named(out, "x")
  expect_equal(unname(out), 0.1)
})

test_that("variety_solve can return both stanfit and results", {
  p <- mp("x")
  fake_df <- data.frame(x = c(1.01, 0.99), lp__ = c(0, 0))
  fake_fit <- fake_cmdstan_fit(fake_df)

  testthat::local_mocked_bindings(
    rvnorm = function(...) fake_fit,
    .package = "vnorm"
  )

  out <- variety_solve(p, n = 2, sig_digit = 2, stanfit = TRUE)

  expect_true(is.list(out))
  expect_true(all(c("stanfit", "results") %in% names(out)))
  expect_identical(out$stanfit, fake_fit)
  expect_equal(unname(out$results), 1.0)
})
