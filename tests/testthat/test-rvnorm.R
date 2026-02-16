test_that("rvnorm() respects RNG seeds", {
  skip_on_cran()
  set.seed(123)
  samps1 <- rvnorm(1000, mp("x^2 + y^2 -1"), sd = 0.1)
  set.seed(123)
  samps2 <- rvnorm(1000, mp("x^2 + y^2 -1"), sd = 0.1)
  expect_equal(samps1, samps2)
})

test_that("rvnorm() has stated output", {
  skip_on_cran()
  samps1 <- rvnorm(1000, mp("x^2 + y^2 -1"), sd = 0.1)
  expect_true(is.data.frame(samps1))
  samps2 <- rvnorm(1000, mp("x^2 + y^2 -1"), sd = 0.1, output = "tibble")
  expect_true(tibble::is_tibble(samps2))
  samps3 <- rvnorm(1000, mp("x^2 + y^2 -1"), sd = 0.1, output = "stanfit")
  expect_true(inherits(samps3, "CmdStanMCMC"))
})

test_that("rvnorm(): code_only = TRUE works",{
  skip_on_cran()
  stan_code <- rvnorm(1000, mp("x^2 + y^2 -1"), sd = 0.1, w = 5, code_only = TRUE)
  expected_result <- "data {
  real<lower=0> si;
}

parameters {
  real<lower=-5,upper=5> x;
  real<lower=-5,upper=5> y;
}

transformed parameters {
  real g = x^2+y^2-1;
  real ndg = sqrt(4*x^2+4*y^2);
}

model {
  target += normal_lpdf(0.00 | g/ndg, si);
}"
  expect_equal(stan_code, expected_result)
})

test_that("rvnorm(): pre_compiled == TRUE works", {
  poly <- mp("x^2 + y^2 - 1")

  # With pre_compiled = TRUE
  result_precompiled <- rvnorm(n = 10, poly = poly, sd = 0.05, pre_compiled = TRUE)
  expect_true(is.data.frame(result_precompiled))

})

test_that("rvnorm() correctly replaces variables in Stan output if poly has non-standard vars", {
  poly <- mp("a^2 + b^2 - 1")
  result <- rvnorm(n = 10, poly = poly, sd = 0.05)

  # Check if output data frame column names match variables in the polynomial
  expect_true(all(mpoly::vars(poly) %in% colnames(result)))
})

test_that("rvnorm() errors for invalid poly input", {
  expect_error(
    rvnorm(n = 10, poly = 123, sd = 0.05),
    "`poly` should be either a character vector, mpoly, or mpolyList."
  )
})


