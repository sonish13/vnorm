#' 1d Varieties in 2d
#'
#' Draw a variety using the biroot and ggplot2.
#'
#' @section Aesthetics: `geom_variety()` understands the following aesthetics
#'   (required aesthetics are in bold):
#'
#'   - **x**
#'   - **y**
#'   - alpha
#'   - color
#'   - group
#'   - linetype
#'   - linewidth
#'   - subgroup
#'
#' @section Computed variables:
#'
#'   \describe{ \item{probs}{The probability associated with the highest density region, specified
#'   by `probs` argument.} }
#'
#' @inheritParams ggplot2::geom_path
#' @seealso [geom_path()]
#' @name geom_variety
#' @examples
#' # Basic plot
#' p <- mp("x^2 + y^2 - 1")
#' ggplot() +
#'   geom_variety(poly = p, xlim = c(-2,2), ylim = c(-2,2)) +
#'   coord_equal()
#'
#' p <- mpoly::lissajous(5, 5, 0, 0)
#' ggplot() +
#'   geom_variety(poly = p, xlim = c(-2,2), ylim = c(-2,2), n = 251) +
#'   coord_equal()
#'






#' @rdname geom_variety
#' @export
stat_variety <- function(
  mapping = NULL,
  data = NULL,
  geom = GeomVariety,
  position = "identity",
  ...,
  poly = NULL,
  n = 101,
  nx = n,
  ny = n,
  xlim = NULL,
  ylim = NULL,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
) {

  if (is.null(data)) data <- ensure_nonempty_data

  layer(
    data = data,
    mapping = mapping,
    stat = StatVariety,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      poly = poly,
      n = n,
      nx = nx,
      ny = ny,
      xlim = xlim,
      ylim = ylim,
      na.rm = na.rm,
      ...
    )
  )
}




#' @rdname geom_variety
#' @format NULL
#' @usage NULL
#' @export
StatVariety <- ggproto(
  "StatVariety",
  Stat,

  compute_group = function(
    self, data, scales, na.rm = FALSE,
    poly, n = 100, nx = n, ny = n, xlim = NULL, ylim = NULL
  ) {

    rangex <- xlim %||% scales$x$dimension()
    rangey <- ylim %||% scales$y$dimension()
    dfxyz <- poly_to_df(poly, rangex, rangey, nx, ny)
    isolines <- xyz_to_isolines(dfxyz, 0)
    iso_to_path(isolines, data$group[1])

  }
)




#' @rdname geom_variety
#' @format NULL
#' @usage NULL
#' @export
GeomVariety <- ggproto(
  "GeomVariety",
  GeomPath,
  default_aes = aes(
    colour = "#3366FF",
    linewidth = 0.5,
    linetype = 1,
    alpha = NA
  )
)





#' @rdname geom_variety
#' @export
geom_variety <- function(
    mapping = NULL,
    data = NULL,
    stat = "variety",
    position = "identity",
    ...,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE
) {

  if (is.null(data)) data <- ensure_nonempty_data

  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomVariety,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}




poly_to_df <- function( poly, xlim, ylim, nx, ny ) {

  if (!is.mpoly(poly)) poly <- mp(poly)
  f <- as.function( x = poly, varorder = c("x", "y"), silent = TRUE )
  df <- expand.grid(
    "x" = seq(xlim[1], xlim[2], length.out = nx),
    "y" = seq(ylim[1], ylim[2], length.out = ny)
  )

  df$z <- with(df, f(cbind(x, y)))
  df
}
