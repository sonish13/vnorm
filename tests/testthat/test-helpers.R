test_that("get_listed_coeficients converts coefficient names to Stan-style names", {
  coefs <- c(1, 2, 3)
  names(coefs) <- c("1", "x", "x^2 y")

  out <- get_listed_coeficients(coefs)

  expect_true(is.list(out))
  expect_equal(names(out), c("b1", "bx", "bx2y"))
  expect_equal(unname(unlist(out)), unname(coefs))
})

test_that("make_coefficients_data fills missing basis coefficients with zero", {
  p <- mp("x^2 + 2 y")
  out <- make_coefficients_data(p, num_of_vars = 2, deg = 2)

  expect_true(is.list(out))
  expect_true(all(c("bx2", "by", "bxy", "by2", "b1") %in% names(out)))
  expect_equal(out[["bx2"]], 1)
  expect_equal(out[["by"]], 2)
  expect_equal(out[["bxy"]], 0L)
  expect_equal(out[["by2"]], 0L)
})

test_that("get_coefficeints_data supports mpoly and mpolyList", {
  p1 <- mp("x^2 + y")
  p2 <- mp(c("x^2 + y", "x - 1"))

  out1 <- get_coefficeints_data(p1)
  out2 <- get_coefficeints_data(p2)

  expect_true(is.list(out1))
  expect_true(all(c("bx2", "by") %in% names(out1)))

  expect_true(is.list(out2))
  expect_true(any(grepl("_1$", names(out2))))
  expect_true(any(grepl("_2$", names(out2))))
})

test_that("variable remapping helpers round-trip names", {
  p <- mp("a^2 + b - 1")
  remapped <- check_and_replace_vars(p)

  expect_true(inherits(remapped$polynomial, "mpoly"))
  expect_true(setequal(mpoly::vars(remapped$polynomial), c("x", "y")))
  expect_equal(unname(unlist(remapped$mapping[c("x", "y")])), c("a", "b"))

  df <- data.frame(x = 1, y = 2, lp__ = 3, check.names = FALSE)
  renamed <- rename_output_df(df, remapped$mapping)
  expect_true(all(c("a", "b", "lp__") %in% names(renamed)))
})

test_that("Stan string helpers return expected formats", {
  expect_match(mpoly_to_stan(mp("x^2 + y")), "x\\^2")

  ml <- mpolyList(mp("x^2 + y"), mp("x - 1"))
  ml_str <- mpolyList_to_stan(ml)
  expect_type(ml_str, "character")
  expect_match(ml_str, ",")

  d_str <- get_derivative("x", num_of_vars = 2, deg = 2)
  expect_type(d_str, "character")
  expect_true(nchar(d_str) > 0)
})
