


poly_to_df <- function(
    poly,
    xlim = c(-1,1),
    ylim = c(-1,1),
    n = 101,
    nx = n,
    ny = n) {
  if (!is.mpoly(poly)) poly <- mp(poly)
  f <- as.function( x = poly, varorder = c("x", "y"), silent = TRUE )
  df <- expand.grid(
    "x" = seq(xlim[1], xlim[2], length.out = nx),
    "y" = seq(ylim[1], ylim[2], length.out = ny)
  )

  df$z <- with(df, f(cbind(x, y)))
  df
}



StatVariety <- ggproto(
  "StatVariety",
  StatContour,
  setup_data = function(data, params) {
    if (!is.mpoly(poly)) poly <- mp(poly)
    f <- as.function( x = poly, varorder = c("x", "y"), silent = TRUE )
    df <- expand.grid(
      "x" = seq(xlim[1], xlim[2], length.out = nx),
      "y" = seq(ylim[1], ylim[2], length.out = ny)
    )
    df$z <- with(df, f(cbind(x, y)))
    df
  }
)

stat_variety <- function(
  mapping = NULL,
  data = NULL,
  geom = GeomVariety, position = "identity",
  poly = NULL,
  n = 101,
  nx = n,
  ny = n,
  xlim = c(-1,1),
  ylim = c(-1,1),
  na.rm = FALSE,
  ...,
  show.legend = NA,
  inherit.aes = TRUE
) {

  if (is.null(data)) data <- get("ensure_nonempty_data", envir = asNamespace("ggplot2"))

  layer(
    stat = StatVariety,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      poly = poly,
      n = 101,
      nx = n,
      ny = n,
      xlim = c(-1,1),
      ylim = c(-1,1),
      na.rm = na.rm,
      ...
    )
  )
}


GeomVariety <- ggproto(
  "GeomVariety",
  GeomContour,
  default_aes = aes(
    weight = 1,
    colour = "#123456",
    linewidth = 0.5,
    linetype = 1,
    alpha = NA,
    breaks = 0
  )
)

geom_variety <- function(
    mapping = NULL,
    data = NULL,
    stat = StatVariety, position = "identity",
    poly = NULL,
    n = 101,
    nx = n,
    ny = n,
    xlim = c(-1,1),
    ylim = c(-1,1),
    ...,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE
) {

  if (is.null(data)) data <- get("ensure_nonempty_data", envir = asNamespace("ggplot2"))

  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomVariety,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      poly = poly,
      n = 101,
      nx = n,
      ny = n,
      xlim = c(-1,1),
      ylim = c(-1,1),
      na.rm = na.rm,
      ...
    )
  )
}

# p <- mp("x^2 + y^2 - 1")
# ggplot() +
#   geom_variety(poly = p)
#
# dft <- poly_to_df(p)
# dft |>
#   ggplot(aes(x , y , z = z )) +
#   geom_contour(breaks = 0)
