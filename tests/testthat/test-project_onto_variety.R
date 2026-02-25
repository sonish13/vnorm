test_that("project_onto_variety reduces residual for a simple curve", {
  p <- mp("x^2 + y^2 - 1")
  pf <- as.function(p, varorder = c("x", "y"), silent = TRUE)
  x0 <- c(1.4, 0.3)

  x_proj <- project_onto_variety(
    x0, p, al = c(1, 1), adaptive = TRUE, dt = 0.02, tol = 1e-4
  )

  expect_type(x_proj, "double")
  expect_length(x_proj, 2)
  expect_true(all(is.finite(x_proj)))
  expect_lt(abs(pf(x_proj)), abs(pf(x0)))
})

test_that("project_onto_variety supports matrix and data.frame batch input", {
  p <- mp("x^2 + y^2 - 1")
  pts <- rbind(c(1.2, 0.1), c(-1.3, 0.4), c(0.2, -1.4))
  colnames(pts) <- c("x", "y")

  mat_out <- project_onto_variety(
    pts, p, al = c(1, 1), adaptive = FALSE, dt = 0.02, tol = 1e-4
  )
  df_out <- project_onto_variety(
    as.data.frame(pts), p, al = c(1, 1), adaptive = FALSE, dt = 0.02, tol = 1e-4
  )

  expect_true(is.matrix(mat_out))
  expect_equal(dim(mat_out), dim(pts))
  expect_true(all(is.finite(mat_out)))

  expect_true(is.data.frame(df_out))
  expect_equal(dim(df_out), dim(pts))
  expect_true(all(is.finite(as.matrix(df_out))))
})

test_that("project_onto_variety warns and rebuilds helpers when partial helpers are supplied", {
  p <- mp("x^2 + y^2 - 1")
  gfunc <- as.function(p, varorder = c("x", "y"), silent = TRUE)

  expect_warning(
    project_onto_variety(c(1.1, 0.2), p, al = c(1, 1), gfunc = gfunc, tol = 1e-3),
    "computed internally"
  )
})

test_that("alternative projection methods return finite points", {
  p <- mp("x^2 + y^2 - 1")
  x0 <- c(1.3, 0.2)

  x_lag <- project_onto_variety_lagrange(x0, p, tol = 1e-4)
  x_gd <- suppressWarnings(project_onto_variety_gradient_descent(
    x0, p, method = "fixed", ga = 0.01, maxit = 200, tol = 1e-4
  ))
  x_newt <- suppressWarnings(project_onto_variety_newton(
    x0, p, method = "fixed", ga = 0.1, maxit = 50, tol = 1e-4
  ))

  expect_true(all(is.finite(x_lag)))
  expect_true(all(is.finite(x_gd)))
  expect_true(all(is.finite(x_newt)))
  expect_equal(length(x_lag), 2)
  expect_equal(length(x_gd), 2)
  expect_equal(length(x_newt), 2)
})

test_that("project_onto_variety validates required inputs and patch vector", {
  p <- mp("x^2 + y^2 - 1")

  expect_error(project_onto_variety(poly = p), "`x0` must be supplied")
  expect_error(project_onto_variety(c(1, 0)), "`poly` must be supplied")
  expect_error(
    project_onto_variety(c(1, 0), p, al = c(NA_real_, 1)),
    "`al` must be a numeric vector"
  )
})

test_that("project_onto_variety supports fixed-step path with messages", {
  p <- mp("x^2 + y^2 - 1")
  expect_message(
    project_onto_variety(
      c(1.2, 0.1), p, al = c(1, 1), adaptive = FALSE, dt = 0.05, message = TRUE
    )
  )
})

test_that("lagrange method supports optim-based branch and tibble batch input", {
  p <- mp("x^2 + y^2 - 1")
  x0 <- c(1.15, 0.2)

  x_lag_bfgs <- suppressWarnings(
    project_onto_variety_lagrange(x0, p, method = "BFGS", maxit = 50)
  )
  expect_true(all(is.finite(x_lag_bfgs)))

  pts_tbl <- tibble::tibble(x = c(1.1, -1.2), y = c(0.2, 0.3))
  out_tbl <- project_onto_variety_lagrange(pts_tbl, p, method = "newton")
  expect_true(tibble::is_tibble(out_tbl))
  expect_equal(dim(out_tbl), dim(pts_tbl))
})

test_that("gradient descent and newton methods cover line-search branches", {
  p <- mp("x^2 + y^2 - 1")
  x0 <- c(1.2, 0.25)

  x_gd_line <- suppressWarnings(project_onto_variety_gradient_descent(
    x0, p, method = "line", maxit = 150, tol = 1e-4
  ))
  x_gd_opt <- suppressWarnings(project_onto_variety_gradient_descent(
    x0, p, method = "optimal", maxit = 150, tol = 1e-4
  ))
  x_newt_line <- suppressWarnings(project_onto_variety_newton(
    x0, p, method = "line", maxit = 50, tol = 1e-4
  ))

  expect_true(all(is.finite(x_gd_line)))
  expect_true(all(is.finite(x_gd_opt)))
  expect_true(all(is.finite(x_newt_line)))
})

test_that("gradient descent and newton support data frame batch inputs", {
  p <- mp("x^2 + y^2 - 1")
  pts <- data.frame(x = c(1.1, -1.2), y = c(0.15, 0.4))

  gd_out <- suppressWarnings(project_onto_variety_gradient_descent(
    pts, p, method = "fixed", maxit = 100, tol = 1e-4
  ))
  newt_out <- suppressWarnings(project_onto_variety_newton(
    pts, p, method = "fixed", maxit = 50, tol = 1e-4
  ))

  expect_true(is.data.frame(gd_out))
  expect_true(is.data.frame(newt_out))
  expect_equal(dim(gd_out), dim(pts))
  expect_equal(dim(newt_out), dim(pts))
})
