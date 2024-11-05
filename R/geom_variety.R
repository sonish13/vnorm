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
#'   \describe{
#'     \item{probs}{The probability associated with the highest density region, specified
#'     by `probs` argument.}
#'   }
#'
#' @inheritParams ggplot2::geom_path
#' @seealso [geom_path()]
#' @name geom_variety
#' @examples
#' # Basic plots
#' p <- mpoly::mp("y - x")
#' ggplot() +
#'   geom_variety(poly = p, xlim = c(-2, 2), ylim = c(-2, 2))
#'
#' p <- mpoly::mp("y - x^2")
#' ggplot() +
#'   geom_variety(poly = p, xlim = c(-2, 2), ylim = c(-2, 2))
#'
#' p <- mpoly::mp("x^2 + y^2 - 1")
#' ggplot() +
#'   geom_variety(poly = p, xlim = c(-2, 2), ylim = c(-2, 2))
#'
#' p1 <- mp("x^2 + y^2 - 1")
#' p2 <- mp("y - x")
#' poly <- mpolyList(p1, p2)
#' ggplot() +
#'   geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   coord_equal()
#'
#' (p <- lissajous(7, 7, 0, 0))
#' ggplot() +
#'   geom_variety(poly = p, xlim = c(-2, 2), ylim = c(-2, 2), n = 201)
#'
#' ## Setting limits
#' p <- mpoly::mp("y - x^2")
#' ggplot() +
#'   geom_variety(poly = p) +
#'   coord_equal()
#'
#' ggplot() +
#'   geom_variety(poly = p, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   coord_equal()

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
  # workaround for data = NULL
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
    poly, n = 101, nx = n, ny = n, xlim = NULL, ylim = NULL
  ) {
    rangex <- if (is.null(xlim)) {
      if (!is.null(scales$x)) scales$x$dimension() else c(-1, 1)
    } else xlim

    rangey <- if (is.null(ylim)) {
      if (!is.null(scales$y)) scales$y$dimension() else c(-1, 1)
    } else ylim
      # calculation for mpoly
    if (is.mpoly(poly)) {
      dfxyz <- poly_to_df(poly, rangex, rangey, nx, ny)
      isolines <- xyz_to_isolines(dfxyz, 0)
      df <- iso_to_path(isolines, data$group[1])
      df$Variety <- poly |> mpoly_to_stan() |> as.character()
      return(df)
    } else if (is.mpolyList(poly)) {
      #calculation for mpolyList
      data_list <- lapply(seq_along(poly), function(i) {
        dfxyz <- poly_to_df(poly[[i]], rangex, rangey, nx, ny)
        isolines <- xyz_to_isolines(dfxyz, 0)
        df <- iso_to_path(isolines, paste0(data$group[1], "_", i))
        df$Variety <- poly[[i]] |> mpoly_to_stan() |> as.character()
        return(df)
      })
      combined_data <- dplyr::bind_rows(data_list)
      return(combined_data)
    } else {
      stop("Input must be either an mpoly or mpolyList object.")
    }
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
    colour = "red",
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
    poly,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE
) {
  # workaround for data = NULL
  if (is.null(data)) data <- ensure_nonempty_data

  mapping <- if (is.null(mapping)) {
    aes(group = after_stat(group), colour = after_stat(Variety))
  } else {
    modifyList(mapping, aes(group = after_stat(group), colour = after_stat(Variety)))
  }

  layer(
    stat = stat,
    data = data,
    mapping = mapping,
    geom = GeomVariety,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      poly = poly,
      na.rm = na.rm,
      ...
    )
  )
}

# Function to get the value of polynomial in a grid
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

# Used for printing Variety names
mpoly_to_stan <- function(mpoly) {
  p <- get("print.mpoly", asNamespace("mpoly"))
  p(mpoly, stars = TRUE, silent = TRUE, plus_pad = 0L, times_pad = 0L) |>
    stringr::str_replace_all("[*]{2}", "^")
}
