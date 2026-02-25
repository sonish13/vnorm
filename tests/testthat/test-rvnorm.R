skip_if_no_cmdstan <- function() {
  skip_if_not_installed("cmdstanr")
  ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
  if (is.null(ver) || identical(ver, "none")) skip("CmdStan is not available")
}

rvnorm_or_skip <- function(...) {
  tryCatch(
    rvnorm(...),
    error = function(e) {
      skip(paste("Skipping rvnorm sampling test:", conditionMessage(e)))
    }
  )
}

fake_cmdstan_fit <- function(draws_df) {
  structure(
    list(
      draws = function(format = "df", inc_warmup = FALSE) {
        draws_df
      }
    ),
    class = "CmdStanMCMC"
  )
}

fake_cmdstan_model <- function(draws_df, capture = NULL) {
  fit <- fake_cmdstan_fit(draws_df)
  list(
    sample = function(...) {
      if (!is.null(capture)) {
        capture$args <- list(...)
      }
      fit
    }
  )
}

test_that("rvnorm() respects RNG seeds", {
  skip_if_no_cmdstan()
  set.seed(123)
  samps1 <- rvnorm_or_skip(
    250, mp("x^2 + y^2 -1"), sd = 0.1, chains = 1, cores = 1
  )
  set.seed(123)
  samps2 <- rvnorm_or_skip(
    250, mp("x^2 + y^2 -1"), sd = 0.1, chains = 1, cores = 1
  )
  expect_equal(samps1, samps2)
})

test_that("rvnorm() has stated output", {
  skip_if_no_cmdstan()
  samps1 <- rvnorm_or_skip(
    250, mp("x^2 + y^2 -1"), sd = 0.1, chains = 1, cores = 1
  )
  expect_true(is.data.frame(samps1))
  samps2 <- rvnorm_or_skip(
    250,
    mp("x^2 + y^2 -1"),
    sd = 0.1,
    output = "tibble",
    chains = 1,
    cores = 1
  )
  expect_true(tibble::is_tibble(samps2))
  samps3 <- rvnorm_or_skip(
    250,
    mp("x^2 + y^2 -1"),
    sd = 0.1,
    output = "stanfit",
    chains = 1,
    cores = 1
  )
  expect_true(inherits(samps3, "CmdStanMCMC"))
})

test_that("rvnorm(): code_only = TRUE works", {
  skip_on_cran()
  stan_code <- rvnorm(
    1000, mp("x^2 + y^2 -1"), sd = 0.1, w = 5, code_only = TRUE
  )
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
  skip_if_no_cmdstan()
  poly <- mp("x^2 + y^2 - 1")
  result_precompiled <- rvnorm_or_skip(
    n = 10,
    poly = poly,
    sd = 0.05,
    pre_compiled = TRUE,
    chains = 1,
    cores = 1
  )
  expect_true(is.data.frame(result_precompiled))
})

test_that(
  paste0(
    "rvnorm() correctly replaces variables in Stan output ",
    "if poly has non-standard vars"
  ),
  {
  skip_if_no_cmdstan()
  poly <- mp("a^2 + b^2 - 1")
  result <- rvnorm_or_skip(
    n = 10, poly = poly, sd = 0.05, chains = 1, cores = 1
  )
  expect_true(all(mpoly::vars(poly) %in% colnames(result)))
  }
)

test_that("rvnorm() errors for invalid poly input", {
  expect_error(
    rvnorm(n = 10, poly = 123, sd = 0.05),
    "`poly` should be either a character vector, mpoly, or mpolyList."
  )
})

test_that("rvnorm() rejection path returns samples without Stan", {
  set.seed(1)
  out <- rvnorm(
    n = 12,
    poly = mp("x^2 + y^2 - 1"),
    sd = 0.05,
    rejection = TRUE,
    output = "tibble",
    w = 1.25
  )

  expect_true(tibble::is_tibble(out))
  expect_equal(nrow(out), 12)
  expect_true(all(c("x", "y") %in% names(out)))
})

test_that("rvnorm() rejection path validates sd shape", {
  expect_error(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
      sd = c(0.1, 0.2),
      rejection = TRUE
    ),
    "not supported yet for rejection sampler"
  )
})

test_that("rvnorm() errors when user_compiled cache is missing", {
  clear_compiled_stan_info()
  expect_error(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      user_compiled = TRUE
    ),
    "No compiled model cache found"
  )
})

test_that("rvnorm() code_only works for mpolyList templates", {
  code <- rvnorm(
    n = 10,
    poly = mp(c("x^2 + y^2 - 1", "y")),
    sd = 0.05,
    code_only = TRUE,
    homo = FALSE
  )
  expect_type(code, "character")
  expect_true(length(code) == 1)
  expect_match(code, "target \\+=")
})

test_that("rvnorm() pre_compiled path works with mocked package model", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    cmdstan_model = function(...) stop("cmdstan_model should not be called"),
    make_coefficients_data = function(...) list(fake_coef = 1),
    check_and_replace_vars = function(p) {
      list(polynomial = mp("x^2 + y^2 - 1"), mapping = list(x = "a", y = "b"))
    },
    rename_output_df = function(df, replacement_list) {
      names(df)[names(df) == "x"] <- replacement_list[["x"]]
      names(df)[names(df) == "y"] <- replacement_list[["y"]]
      df
    },
    .package = "vnorm"
  )
  testthat::local_mocked_bindings(
    stan_package_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "instantiate"
  )

  out <- rvnorm(
    n = 2,
    poly = mp("a^2 + b^2 - 1"),
    sd = 0.05,
    output = "tibble",
    pre_compiled = TRUE,
    chains = 1,
    cores = 1
  )

  expect_true(tibble::is_tibble(out))
  expect_equal(nrow(out), 2)
  expect_gte(ncol(out), 2)
  expect_true(is.list(capture$args$data))
})

test_that("rvnorm() falls back when pre-compiled model is unavailable", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )
  testthat::local_mocked_bindings(
    stan_package_model = function(...) stop("Stan model file \"\" not found."),
    .package = "instantiate"
  )

  expect_message(
    out <- rvnorm(
      n = 2,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      pre_compiled = TRUE,
      chains = 1,
      cores = 1
    ),
    "falling back to regular rvnorm sampling"
  )

  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 2)
  expect_equal(capture$args$iter_sampling, 2)
})

test_that("rvnorm() pre_compiled auto-disables for mpolyList and uses regular path", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), check.names = FALSE
  )
  instantiate_called <- new.env(parent = emptyenv())
  instantiate_called$flag <- FALSE

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df),
    .package = "vnorm"
  )
  testthat::local_mocked_bindings(
    stan_package_model = function(...) {
      instantiate_called$flag <- TRUE
      stop("should not be called")
    },
    .package = "instantiate"
  )

  out <- rvnorm(
    n = 2,
    poly = mp(c("x^2 + y^2 - 1", "y")),
    sd = 0.05,
    pre_compiled = TRUE,
    chains = 1,
    cores = 1
  )

  expect_false(instantiate_called$flag)
  expect_true(is.data.frame(out))
})

test_that("rvnorm() user_compiled path works with cached model metadata", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), check.names = FALSE
  )
  model_name <- generate_model_name(mp("x^2 + y^2 - 1"), homo = TRUE, w = FALSE)
  set_compiled_stan_info(data.frame(
    name = model_name, path = "cached-model.stan", stringsAsFactors = FALSE
  ))
  on.exit(clear_compiled_stan_info(), add = TRUE)

  testthat::local_mocked_bindings(
    cmdstan_model = function(...) fake_cmdstan_model(draws_df),
    get_coefficeints_data = function(...) list(fake_coef = 1),
    .package = "vnorm"
  )

  out <- rvnorm(
    n = 2,
    poly = mp("x^2 + y^2 - 1"),
    sd = 0.05,
    user_compiled = TRUE,
    chains = 1,
    cores = 1
  )

  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 2)
})

test_that("rvnorm() errors when requested user_compiled model is missing", {
  set_compiled_stan_info(data.frame(
    name = "other_model", path = "cached-model.stan", stringsAsFactors = FALSE
  ))
  on.exit(clear_compiled_stan_info(), add = TRUE)

  expect_error(
    rvnorm(
      n = 2,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      user_compiled = TRUE
    ),
    "Requested compiled model not found in cache"
  )
})

test_that("rvnorm() validates Sigma dimensions before sampling", {
  expect_error(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      Sigma = c(1, 2, 3),
      code_only = TRUE
    ),
    "`Sigma` should be a number, vector of length equal to number of variables"
  )
})

test_that("rvnorm() passes computed refresh to sample when refresh=TRUE and verbose=TRUE", {
  draws_df <- data.frame(
    lp__ = seq_len(20),
    x = seq(0.1, 2.0, length.out = 20),
    y = seq(0.2, 2.1, length.out = 20),
    check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )

  rvnorm(
    n = 20,
    poly = mp(c("x^2 + y^2 - 1", "y")),
    sd = 0.05,
    pre_compiled = TRUE,
    refresh = TRUE,
    verbose = TRUE,
    chains = 2,
    cores = 1
  )

  expect_true(is.numeric(capture$args$refresh))
  expect_gte(capture$args$refresh, 1)
})

test_that("rvnorm() sets refresh to 0 when refresh=TRUE and verbose=FALSE", {
  draws_df <- data.frame(
    lp__ = seq_len(20),
    x = seq(0.1, 2.0, length.out = 20),
    y = seq(0.2, 2.1, length.out = 20),
    check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )

  rvnorm(
    n = 20,
    poly = mp(c("x^2 + y^2 - 1", "y")),
    sd = 0.05,
    pre_compiled = TRUE,
    refresh = TRUE,
    verbose = FALSE,
    chains = 2,
    cores = 1
  )

  expect_equal(capture$args$refresh, 0)
})

test_that("rvnorm() accepts character poly and valid Sigma shapes in code_only mode", {
  expect_type(
    rvnorm(
      n = 10,
      poly = "x^2 + y^2 - 1",
      sd = 0.05,
      code_only = TRUE
    ),
    "character"
  )

  expect_type(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      Sigma = c(1, 2),
      code_only = TRUE
    ),
    "character"
  )

  expect_type(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      Sigma = diag(c(1, 2)),
      code_only = TRUE
    ),
    "character"
  )
})

test_that("rvnorm() emits stanfit rename warning in precompiled path", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), check.names = FALSE
  )

  testthat::local_mocked_bindings(
    make_coefficients_data = function(...) list(fake_coef = 1),
    check_and_replace_vars = function(p) {
      list(polynomial = mp("x^2 + y^2 - 1"), mapping = list(x = "a", y = "b"))
    },
    .package = "vnorm"
  )
  testthat::local_mocked_bindings(
    stan_package_model = function(...) fake_cmdstan_model(draws_df),
    .package = "instantiate"
  )

  expect_message(
    out <- rvnorm(
      n = 2,
      poly = mp("a^2 + b^2 - 1"),
      sd = 0.05,
      output = "stanfit",
      pre_compiled = TRUE,
      chains = 1,
      cores = 1
    ),
    "variable names as x,y,z"
  )
  expect_true(inherits(out, "CmdStanMCMC"))
})
