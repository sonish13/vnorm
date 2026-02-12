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
#' @rdname geom_variety
#'
#' @examples
#'
#' library(ggplot2)
#' library(mpoly)
#'
#' # 1) Lemniscate of Bernoulli (∞-shaped)
#' p1 <- mp("(x^2 + y^2)^2 - 2 x^2")
#' ggplot() +
#'   geom_variety(poly = p1, xlim = c(-2, 2), ylim = c(-2, 2), n = 401) +
#'   coord_equal()
#'
#' # 2) Folium of Descartes (loop with an asymptote)
#' p2 <- mp("x^3 + y^3 - 3 x y")
#' ggplot() +
#'   geom_variety(poly = p2, xlim = c(-2, 3), ylim = c(-2, 3), n = 401) +
#'   coord_equal()
#'
#' # 3) "Heart" curve (classic implicit heart)
#' p3 <- mp("(x^2 + y^2 - 1)^3 - x^2 y^3")
#' ggplot() +
#'   geom_variety(poly = p3, xlim = c(-2, 2), ylim = c(-2, 2), n = 401) +
#'   coord_equal()
#'
#' # 4) A 2-polynomial system (mpolyList): circle ∩ hyperbola-like branch overlay
#' #    Useful example showing multiple contours from one call.
#' ps <- mp( c( "x^2 + y^2 - 1", "x y - 0.25" ) )
#' ggplot() +
#'   geom_variety(poly = ps, xlim = c(-2, 2), ylim = c(-2, 2), n = 401) +
#'   coord_equal()
#'
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
    shift = 0,
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
      shift = shift,
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
    poly, n = 101, nx = n, ny = n, xlim = NULL, ylim = NULL, shift = 0
  ) {
    rangex <- if (is.null(xlim)) {
      if (!is.null(scales$x)) scales$x$dimension() else c(-1, 1)
    } else xlim

    rangey <- if (is.null(ylim)) {
      if (!is.null(scales$y)) scales$y$dimension() else c(-1, 1)
    } else ylim

    if (is.mpoly(poly)) {
      dfxyz <- poly_to_df(poly, rangex, rangey, nx, ny, shift)
      check_sign_warning(dfxyz$z, shift)
      isolines <- xyz_to_isolines(dfxyz, 0)
      df <- iso_to_path(isolines, data$group[1])
      df$Polynomial <- as.character(mpoly_to_stan(poly))
      return(df)
    } else if (is.mpolyList(poly)) {
      data_list <- lapply(seq_along(poly), function(i) {
        dfxyz <- poly_to_df(poly[[i]], rangex, rangey, nx, ny, shift)
        check_sign_warning(dfxyz$z, shift)
        isolines <- xyz_to_isolines(dfxyz, 0)
        df <- iso_to_path(isolines, paste0(data$group[1], "_", i))
        df$Polynomial <- as.character(mpoly_to_stan(poly[[i]]))
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
    shift = 0,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE
) {
  if (is.null(data)) {
    data <- ensure_nonempty_data
  }

  mapping <- if (is.null(mapping)) {
    aes(group = after_stat(group), colour = after_stat(Polynomial))
  } else {
    modifyList(mapping, aes(group = after_stat(group), colour = after_stat(Polynomial)))
  }

  layer_obj <- layer(
    stat = stat,
    data = data,
    mapping = mapping,
    geom = GeomVariety,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      poly = poly,
      shift = shift,
      na.rm = na.rm,
      ...
    )
  )

  list(
    layer_obj,
    scale_color_discrete(name = NULL, labels = function(l) parse(text = l))# ,
    # theme(legend.position = "top")
  )
}

poly_to_df <- function(poly, xlim, ylim, nx, ny, shift = 0) {
  if (!is.mpoly(poly)) poly <- mp(poly)
  f <- as.function(x = poly, varorder = c("x", "y"), silent = TRUE)

  df <- expand.grid(
    "x" = seq(xlim[1], xlim[2], length.out = nx),
    "y" = seq(ylim[1], ylim[2], length.out = ny)
  )

  df$z <- with(df, f(cbind(x, y))) + shift
  df
}

check_sign_warning <- function(zvals, shift) {
  signs <- sign(zvals[zvals != 0])
  if (length(signs) == 0) return()
  if (all(signs == 1)) {
    q1 <- quantile(zvals, 0.01)
    warning("All values positive; try setting shift = ", format(-q1, digits = 3))
  } else if (all(signs == -1)) {
    q99 <- quantile(zvals, 0.99)
    warning("All values negative; try setting shift = ", format(-q99, digits = 3))
  }
}

mpoly_to_stan <- function(mpoly) {
  p <- get("print.mpoly", asNamespace("mpoly"))
  result <- p(mpoly, stars = TRUE, silent = TRUE, plus_pad = 0L, times_pad = 0L)
  result <- gsub("[*]{2}", "^", result)
  result
}
