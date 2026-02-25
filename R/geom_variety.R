#' One-Dimensional Varieties in Two Dimensions
#'
#' Plot implicit polynomial varieties with `ggplot2`.
#'
#' @section Aesthetics: [geom_variety()] understands the following aesthetics.
#'   `x` and `y` are computed by the stat, so users typically do not map them
#'   manually:
#'
#'   - x
#'   - y
#'   - alpha
#'   - color
#'   - group
#'   - linetype
#'   - linewidth
#'   - subgroup
#'
#' @section Computed variables:
#'   \describe{
#'     \item{Polynomial}{A parseable label for the polynomial, useful for
#'     `after_stat(Polynomial)` mappings (for example, linetype or colour).}
#'     \item{group}{Contour path group identifier used internally by the layer.}
#'   }
#'
#' @param mapping Aesthetic mappings created with [ggplot2::aes()].
#' @param data Layer data.
#' @param geom The geometric object used to display data; defaults to
#'   [GeomVariety].
#' @param position Position adjustment for the layer.
#' @param ... Additional parameters forwarded to [ggplot2::layer()].
#' @param n Number of grid points used in both x and y directions when `nx` and
#'   `ny` are not supplied.
#' @param nx,ny Number of grid points in the x and y directions.
#' @param xlim,ylim Length-2 numeric vectors giving plotting limits. If not
#'   supplied, limits are taken from the plot scales.
#' @param poly An `mpoly` or `mpolyList` object describing the implicit
#'   polynomial(s) to plot.
#' @param shift A numeric constant added to the evaluated surface before
#'   contouring, i.e. the plotted level set is `poly + shift = 0`.
#'   This is mainly useful when the polynomial does not cross zero on the
#'   plotting grid (for example, `p^2`). If `shift = 0` and all sampled values
#'   have one sign, `geom_variety()` prints a message suggesting a shift value.
#' @param vary_colour Logical. If `TRUE`, map colour to the polynomial label so
#'   users can control per-polynomial colours with `scale_colour_*()`.
#'   Defaults to `FALSE`, which keeps a constant line colour and varies only
#'   linetype across an `mpolyList`.
#'
#' @inheritParams ggplot2::geom_path
#' @seealso [geom_path()]
#' @name geom_variety
#' @rdname geom_variety
#'
#' @examples
#'
#' library("ggplot2")
#'
#' # 1) Ellipse
#' p1 <- mp("x^2 + 4 y^2 - 1")
#' ggplot() +
#'   geom_variety(poly = p1, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   coord_equal()
#'
#' # Works with standard ggplot2 styling
#' ggplot() +
#'   geom_variety(
#'     poly = p1, xlim = c(-2, 2), ylim = c(-2, 2),
#'     colour = "steelblue", linewidth = 0.5
#'   ) +
#'   coord_equal() +
#'   theme_minimal()
#'
#' # 2) Folium of Descartes (singular variety)
#' p2 <- mp("x^3 + y^3 - 3 x y")
#' ggplot() +
#'   geom_variety(poly = p2, xlim = c(-2, 3), ylim = c(-2, 3)) +
#'   coord_equal()
#'
#' # 3) "Heart" curve (classic implicit heart)
#' p3 <- mp("(x^2 + y^2 - 1)^3 - x^2 y^3")
#' ggplot() +
#'   geom_variety(poly = p3, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   coord_equal() +
#'   theme(legend.position = "top")
#'
#' # 4) A 2-polynomial system (mpolyList): circle and xy = 0.25
#' p4 <- mp(c("x^2 + y^2 - 1", "x y - 0.25"))
#' # By default, polynomials differ by linetype (not color).
#' ggplot() +
#'   geom_variety(poly = p4, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   coord_equal()
#'
#' # With different colors (optional)
#' ggplot() +
#'   geom_variety(poly = p4, xlim = c(-2, 2), ylim = c(-2, 2), vary_colour = TRUE) +
#'   coord_equal() +
#'   scale_colour_manual(values = c("steelblue", "firebrick"))
#'
#' # You can also customize linetypes and legend placement with ggplot2 scales/themes
#' ggplot() +
#'   geom_variety(poly = p4, xlim = c(-2, 2), ylim = c(-2, 2), vary_colour = TRUE) +
#'   coord_equal() +
#'   scale_colour_manual(values = c("steelblue", "firebrick")) +
#'   scale_linetype_manual(values = c("solid", "22"), guide = "none") +
#'   theme(legend.position = "top")
#'
#' ## common contouring situations
#' ########################################
#'
#' # 5) Squared polynomial (same zero set, but no sign change on the grid)
#' # geom_variety() will suggest a negative shift when no contour is found.
#' p5 <- mp("x^2 + y^2 - 1")
#' ggplot() +
#'   geom_variety(poly = p5^2, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   coord_equal()
#'
#' # Use the suggested shift (your printed value may differ slightly).
#' ggplot() +
#'   geom_variety(poly = p5^2, xlim = c(-2, 2), ylim = c(-2, 2), shift = -0.001) +
#'   coord_equal()
#'
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
  # Thin wrapper that wires StatVariety into ggplot2::layer().
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
    poly, n = 101, nx = n, ny = n, xlim = NULL, ylim = NULL,
    shift = 0, mul = .05
  ) {
    # Use a denser grid for nonzero shift to reduce fragmented contours.
    nx_eff <- if (shift != 0) max(as.integer(nx), 301L) else as.integer(nx)
    ny_eff <- if (shift != 0) max(as.integer(ny), 301L) else as.integer(ny)

    rangex <- if (is.null(xlim)) {
      if (!is.null(scales$x)) scales$x$dimension() else c(-1, 1)
    } else scales::expand_range(xlim, mul = mul)

    rangey <- if (is.null(ylim)) {
      if (!is.null(scales$y)) scales$y$dimension() else c(-1, 1)
    } else scales::expand_range(ylim, mul = mul)

    if (is.mpoly(poly)) {
      df0 <- poly_to_df(poly, rangex, rangey, nx_eff, ny_eff, shift = 0)
      check_sign_warning(df0$z + shift, shift)
      no_sign_change0 <- !has_strict_sign_change(df0$z)
      if (shift == 0 && no_sign_change0) {
        message("Zero contours were generated")
        return(tibble::tibble())
      }
      df <- variety_paths_with_refinement(
        poly = poly,
        rangex = rangex,
        rangey = rangey,
        nx = nx_eff,
        ny = ny_eff,
        shift = shift,
        group = data$group[1]
      )
      if (shift != 0) {
        if (no_sign_change0) {
          df <- snap_shifted_contours_to_variety(df, poly)
        }
        # Merge near-coincident offset contours produced by shifted surfaces.
        dx <- (rangex[2] - rangex[1]) / max(nx_eff - 1, 1)
        dy <- (rangey[2] - rangey[1]) / max(ny_eff - 1, 1)
        df <- collapse_near_duplicate_contours(df, tol = 0.75 * max(dx, dy))
      }
      df$Polynomial <- as.character(mpoly_to_stan(poly))
      return(df)
    } else if (is.mpolyList(poly)) {
      data_list <- lapply(seq_along(poly), function(i) {
        df0 <- poly_to_df(poly[[i]], rangex, rangey, nx_eff, ny_eff, shift = 0)
        check_sign_warning(df0$z + shift, shift)
        no_sign_change0 <- !has_strict_sign_change(df0$z)
        if (shift == 0 && no_sign_change0) {
          message("Zero contours were generated")
          return(tibble::tibble())
        }
        df <- variety_paths_with_refinement(
          poly = poly[[i]],
          rangex = rangex,
          rangey = rangey,
          nx = nx_eff,
          ny = ny_eff,
          shift = shift,
          group = paste0(data$group[1], "_", i)
        )
        if (shift != 0) {
          if (no_sign_change0) {
            df <- snap_shifted_contours_to_variety(df, poly[[i]])
          }
          # Same de-duplication policy for each polynomial in an mpolyList.
          dx <- (rangex[2] - rangex[1]) / max(nx_eff - 1, 1)
          dy <- (rangey[2] - rangey[1]) / max(ny_eff - 1, 1)
          df <- collapse_near_duplicate_contours(df, tol = 0.75 * max(dx, dy))
        }
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
    vary_colour = FALSE,
    shift = 0,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE
) {
  # Default to linetype differences; colour mapping is opt-in via vary_colour.
  if (is.null(data)) {
    data <- ensure_nonempty_data
  }

  default_mapping <- aes(group = after_stat(group), linetype = after_stat(Polynomial))
  if (isTRUE(vary_colour)) {
    default_mapping <- modifyList(
      default_mapping,
      aes(colour = after_stat(Polynomial))
    )
  }

  mapping <- if (is.null(mapping)) {
    default_mapping
  } else {
    modifyList(mapping, default_mapping)
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

  n_poly <- if (is.mpolyList(poly)) length(poly) else 1L
  out <- list(
    layer_obj,
    ggplot2::scale_linetype_discrete(
      name = NULL,
      labels = function(l) parse(text = l),
      guide = if (isTRUE(vary_colour)) "none" else "legend"
    )
  )
  if (isTRUE(vary_colour)) {
    out <- c(out, list(ggplot2::guides(
      colour = ggplot2::guide_legend(
        title = NULL,
        override.aes = list(linetype = scales::linetype_pal()(max(1L, n_poly)))
      )
    )))
  }
  out
}
