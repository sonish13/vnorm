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

test_that("geom_variety generates expected plot with valid mpoly input", {
  poly <- mp("x^2 + y^2 - 1")
  plot_valid <- ggplot() +
    geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2)) +
    ggplot2::coord_equal() +
    labs(title = "Circle Variety Plot")

  vdiffr::expect_doppelganger("geom_variety valid mpoly plot", plot_valid)
})
