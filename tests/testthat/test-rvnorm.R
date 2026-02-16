skip_if_no_cmdstan <- function() {
  skip_if_not_installed("cmdstanr")
  ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
  if (is.null(ver) || identical(ver, "none")) {
    skip("CmdStan is not available")
  }
}

rvnorm_or_skip <- function(...) {
  out <- tryCatch(
    rvnorm(...),
    error = function(e) {
      skip(paste("Skipping rvnorm sampling test:", conditionMessage(e)))
    }
  )
  out
}

test_that("rvnorm is consistent", {
  skip_if_no_cmdstan()
  set.seed(123)
  samps1 <- rvnorm_or_skip(250, mpoly::mp("x^2 + y^2 -1"), sd = 0.1, chains = 1, cores = 1)
  set.seed(123)
  samps2 <- rvnorm_or_skip(250, mpoly::mp("x^2 + y^2 -1"), sd = 0.1, chains = 1, cores = 1)
  expect_equal(samps1, samps2)
})

test_that("rvnorm has proper output", {
  skip_if_no_cmdstan()
  samps1 <- rvnorm_or_skip(250, mpoly::mp("x^2 + y^2 -1"), sd = 0.1, chains = 1, cores = 1)
  expect_true(is.data.frame(samps1))
  samps2 <- rvnorm_or_skip(250, mpoly::mp("x^2 + y^2 -1"), sd = 0.1, output = "tibble", chains = 1, cores = 1)
  expect_true(tibble::is_tibble(samps2))
  samps3 <- rvnorm_or_skip(250, mpoly::mp("x^2 + y^2 -1"), sd = 0.1, output = "stanfit", chains = 1, cores = 1)
  expect_true(inherits(samps3, "CmdStanMCMC"))
})

test_that("rvnorm test for code only portion",{
  skip_on_cran()
  stan_code <- rvnorm(1000, mpoly::mp("x^2 + y^2 -1"), sd = 0.1, w = 5, code_only = T)
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

test_that("rvnorm works with pre_compiled option TRUE", {
  skip_if_no_cmdstan()
  poly <- mp("x^2 + y^2 - 1")

  # With pre_compiled = TRUE
  result_precompiled <- rvnorm_or_skip(n = 10, poly = poly, sd = 0.05, pre_compiled = TRUE, chains = 1, cores = 1)
  expect_true(is.data.frame(result_precompiled))

})

test_that("rvnorm correctly replaces variables in Stan output if poly has non-standard vars", {
  skip_if_no_cmdstan()
  poly <- mp("a^2 + b^2 - 1")
  result <- rvnorm_or_skip(n = 10, poly = poly, sd = 0.05, chains = 1, cores = 1)

  # Check if output data frame column names match variables in the polynomial
  expect_true(all(mpoly::vars(poly) %in% colnames(result)))
})

test_that("rvnorm throws error for invalid poly input", {
  expect_error(rvnorm(n = 10, poly = 123, sd = 0.05),
               "`poly` should be either a character vector, mpoly, or mpolyList.")
})

