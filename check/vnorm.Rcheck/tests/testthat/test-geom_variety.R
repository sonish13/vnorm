poly_to_df <- function(poly, xlim, ylim, nx, ny) {
  if (!is.mpoly(poly)) poly <- mp(poly)
  f <- as.function(x = poly, varorder = c("x", "y"), silent = TRUE)

  df <- expand.grid(
    "x" = seq(xlim[1], xlim[2], length.out = nx),
    "y" = seq(ylim[1], ylim[2], length.out = ny)
  )

  df$z <- with(df, f(cbind(x, y)))
  df
}

test_that("geom_variety works with mpoly objects", {
  poly <- mp("y - x")
  p <- ggplot() + geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2))
  expect_true(inherits(p, "ggplot"))
})

test_that("geom_variety works with mpolyList objects", {
  poly1 <- mp("x^2 + y^2 - 1")
  poly2 <- mp("y - x")
  poly_list <- mpolyList(poly1, poly2)
  p <- ggplot() + geom_variety(poly = poly_list, xlim = c(-2, 2), ylim = c(-2, 2))
  expect_true(inherits(p, "ggplot"))
})

test_that("geom_variety handles xlim and ylim parameters", {
  poly <- mp("x^2 + y^2 - 1")
  p <- ggplot() + geom_variety(poly = poly, xlim = c(-1, 1), ylim = c(-1, 1))
  expect_equal(p$layers[[1]]$stat_params$xlim, c(-1, 1))
  expect_equal(p$layers[[1]]$stat_params$ylim, c(-1, 1))
})

test_that("geom_variety computes correct paths for simple varieties", {
  poly <- mp("y - x")
  p <- ggplot() + geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2))
  ggplot2::ggplot_build(p)
  df <- poly_to_df(poly, xlim = c(-2, 2), ylim = c(-2, 2), nx = 100, ny = 100)
  expect_true("x" %in% names(df) && "y" %in% names(df) && "z" %in% names(df))
})

test_that("geom_variety generates valid non-empty contour data", {
  poly <- mp("x^2 + y^2 - 1")
  plot_valid <- ggplot() +
    geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2)) +
    ggplot2::coord_equal() +
    labs(title = "Circle Variety Plot")

  built <- ggplot2::ggplot_build(plot_valid)
  dat <- built$data[[1]]
  expect_gt(nrow(dat), 0)
  expect_true(all(c("x", "y", "group") %in% names(dat)))
})

test_that("geom_variety shift generates continuous contours for squared polynomials", {
  poly <- mp("x^2 + y^2 - 1")^2
  p <- ggplot() +
    geom_variety(
      poly = poly,
      xlim = c(-2, 2),
      ylim = c(-2, 2),
      shift = -0.000576
    )

  built <- ggplot2::ggplot_build(p)
  dat <- built$data[[1]]

  expect_gt(nrow(dat), 0)

  points_per_group <- table(dat$group)
  expect_lte(length(points_per_group), 2)
  expect_true(all(points_per_group >= 20))
})

test_that("fragmented shifted contours trigger refinement", {
  poly <- mp("x^2 + y^2 - 1")^2
  df <- variety_paths_with_refinement(
    poly = poly,
    rangex = c(-2, 2),
    rangey = c(-2, 2),
    nx = 101,
    ny = 101,
    shift = -0.000576,
    group = 1
  )

  expect_gt(nrow(df), 0)
  expect_false(is_fragmented_paths(df))
})

test_that("duplicate shifted contours are collapsed", {
  poly <- mp("x^2 + y^2 - 1")^2
  p <- ggplot() +
    geom_variety(
      poly = poly,
      xlim = c(-2, 2),
      ylim = c(-2, 2),
      shift = -0.000576
    )
  dat <- ggplot2::ggplot_build(p)$data[[1]]
  expect_lte(length(unique(dat$group)), 2)
})

test_that("shift suggestion is emitted as message when shift is zero", {
  poly <- mp("(x^2 + y^2 - 1)^2")
  p <- ggplot() + geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2))
  expect_message(ggplot2::ggplot_build(p), "try shift =")
})
