#' 1d Varieties in 2d
#'
#' Draw a variety using ggplot2.
#'
#' @section Aesthetics: [geom_variety()] understands the following aesthetics
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
#' @param geom The geometric object to use display data, defaults to
#'   `GeomVariety`.
#' @param n Number of grid points used in both x and y directions when `nx` and
#'   `ny` are not supplied.
#' @param nx,ny Number of grid points in the x and y directions.
#' @param xlim,ylim Length-2 numeric vectors giving plotting limits.
#' @param poly An `mpoly` or `mpolyList` object describing the implicit
#'   polynomial(s) to plot.
#' @param shift A numeric constant added to the evaluated surface before
#'   contouring, i.e. the plotted level set is `poly + shift = 0`.
#'   Use this when the polynomial does not cross zero on the plotting grid.
#'   If `shift = 0` and all values have one sign, `geom_variety()` prints a
#'   message suggesting a shift value.
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
#' # 3) Folium of Descartes (loop with an asymptote) - variety with singularity
#' p3 <- mp("x^3 + y^3 - 3 x y")
#' ggplot() +
#'   geom_variety(poly = p3, xlim = c(-2, 3), ylim = c(-2, 3)) +
#'   coord_equal()
#'
#' # 4) "Heart" curve (classic implicit heart)
#' p4 <- mp("(x^2 + y^2 - 1)^3 - x^2 y^3")
#' ggplot() +
#'   geom_variety(poly = p4, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   coord_equal() +
#'   theme(legend.position = "top")
#'
#' # 5) A 2-polynomial system (mpolyList): circle ∩ hyperbola-like branch overlay
#' #    Useful example showing multiple contours from one call.
#' ps <- mp( c( "x^2 + y^2 - 1", "x y - 0.25" ) )
#' ggplot() +
#'   geom_variety(poly = ps, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   coord_equal()
#'
#' # 5) Using shift for squared polynomials (same variety, no zero crossing on grid)
#' library(mpoly)
#' p_shift <- mp("x^2 + y^2 - 1")^2
#' ggplot() +
#'   geom_variety(poly = p_shift, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   coord_equal()
#'
#'
#'
#' ## known issues
#' ########################################
#'
#' # 1) singularities Lemniscate of Bernoulli (∞-shaped)
#' p <- mp("(x^2 + y^2)^2 - 2 x^2")
#' ggplot() +
#'   geom_variety(poly = p, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   coord_equal()
#'
#' # 2) non-zero crossing components
#' p <- mp("x^2 + y^2 - 1")
#' ggplot() +
#'   geom_variety(poly = p, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   coord_equal()
#'
#' ggplot() +
#'   geom_variety(poly = p^2, xlim = c(-2, 2), ylim = c(-2, 2)) +
#'   coord_equal()
#'
#' ggplot() +
#'   geom_variety(poly = p^2, xlim = c(-2, 2), ylim = c(-2, 2), shift = -0.00101684) +
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
    poly, n = 101, nx = n, ny = n, xlim = NULL, ylim = NULL, shift = 0, mul = .05
  ) {
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
      check_sign_warning(df0$z, shift)
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
        dx <- (rangex[2] - rangex[1]) / max(nx_eff - 1, 1)
        dy <- (rangey[2] - rangey[1]) / max(ny_eff - 1, 1)
        df <- collapse_near_duplicate_contours(df, tol = 0.75 * max(dx, dy))
      }
      df$Polynomial <- as.character(mpoly_to_stan(poly))
      return(df)
    } else if (is.mpolyList(poly)) {
      data_list <- lapply(seq_along(poly), function(i) {
        df0 <- poly_to_df(poly[[i]], rangex, rangey, nx_eff, ny_eff, shift = 0)
        check_sign_warning(df0$z, shift)
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
    shift = 0,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE
) {
  if (is.null(data)) {
    data <- ensure_nonempty_data
  }

  mapping <- if (is.null(mapping)) {
    aes(group = after_stat(group), linetype = after_stat(Polynomial))
  } else {
    modifyList(mapping, aes(group = after_stat(group), linetype = after_stat(Polynomial)))
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
    ggplot2::scale_linetype_discrete(name = NULL, labels = function(l) parse(text = l))
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

is_fragmented_paths <- function(df) {
  if (nrow(df) == 0) return(FALSE)
  n_per_group <- as.numeric(table(df$group))
  if (length(n_per_group) <= 8) return(FALSE)

  frac_short <- mean(n_per_group < 12)
  frac_short > 0.7
}

mean_nn_distance <- function(df_a, df_b) {
  dx <- outer(df_a$x, df_b$x, "-")
  dy <- outer(df_a$y, df_b$y, "-")
  d2 <- dx * dx + dy * dy
  mean(sqrt(apply(d2, 1, min)))
}

path_proximity <- function(df_a, df_b) {
  max(mean_nn_distance(df_a, df_b), mean_nn_distance(df_b, df_a))
}

collapse_near_duplicate_contours <- function(df, tol) {
  if (nrow(df) == 0) return(df)

  groups <- split(df, df$group)
  if (length(groups) <= 1) return(df)

  g_names <- names(groups)
  centroids <- t(vapply(
    groups,
    function(g) c(mean(g$x), mean(g$y)),
    numeric(2)
  ))
  rownames(centroids) <- g_names

  # Fast path for the common shifted-double case (two near-coincident contours).
  if (length(groups) == 2) {
    g1 <- g_names[1]
    g2 <- g_names[2]
    d_cent <- sqrt(sum((centroids[g1, ] - centroids[g2, ])^2))
    d_path <- path_proximity(groups[[g1]], groups[[g2]])
    if (d_cent <= 3 * tol && d_path <= 3 * tol) {
      k <- if (nrow(groups[[g2]]) > nrow(groups[[g1]])) g2 else g1
      out <- groups[[k]]
      out$group <- factor(out$group)
      out$piece <- as.integer(out$group)
      return(out)
    }
  }

  removed <- stats::setNames(rep(FALSE, length(groups)), g_names)
  keep <- character(0)

  for (i in seq_along(g_names)) {
    g <- g_names[i]
    if (removed[[g]]) next

    cluster <- g
    if (i < length(g_names)) {
      for (j in (i + 1):length(g_names)) {
        h <- g_names[j]
        if (removed[[h]]) next

        d_cent <- sqrt(sum((centroids[g, ] - centroids[h, ])^2))
        if (d_cent > 3 * tol) next

        d_path <- path_proximity(groups[[g]], groups[[h]])
        if (d_path <= 3 * tol) {
          cluster <- c(cluster, h)
        }
      }
    }

    if (length(cluster) == 1) {
      keep <- c(keep, g)
      next
    }

    sizes <- vapply(groups[cluster], nrow, integer(1))
    k <- cluster[which.max(sizes)]
    keep <- c(keep, k)
    removed[setdiff(cluster, k)] <- TRUE
  }

  out <- dplyr::bind_rows(groups[unique(keep)])
  out$group <- factor(out$group)
  out$piece <- as.integer(out$group)
  out
}

variety_paths_with_refinement <- function(poly, rangex, rangey, nx, ny, shift, group) {
  refinement_steps <- c(1L, 2L, 4L)
  best_df <- tibble::tibble()

  for (step in refinement_steps) {
    dfxyz <- poly_to_df(
      poly = poly,
      xlim = rangex,
      ylim = rangey,
      nx = nx * step,
      ny = ny * step,
      shift = shift
    )
    isolines <- xyz_to_isolines(dfxyz, 0)
    df <- iso_to_path(isolines, group)
    best_df <- df

    if (!is_fragmented_paths(df)) {
      return(df)
    }
  }

  best_df
}

check_sign_warning <- function(zvals, shift) {
  signs <- sign(zvals[zvals != 0])
  if (length(signs) == 0) return()

  scale_z <- max(abs(zvals), na.rm = TRUE)
  tol <- sqrt(.Machine$double.eps) * max(1, scale_z)
  near_zero_touch <- min(abs(zvals), na.rm = TRUE) <= tol

  if (all(signs == 1)) {
    q1 <- quantile(zvals, 0.01, na.rm = TRUE)
    if (shift == 0) {
      message(
        "All values are positive on the plotting grid; ",
        "try shift = ", format(-q1, digits = 6), "."
      )
    } else if (near_zero_touch) {
      message("Using shift = ", format(shift, digits = 6), "; duplicate contours merged.")
    } else {
      message("All values positive after applying shift = ", format(shift, digits = 6), ".")
    }
  } else if (all(signs == -1)) {
    q99 <- quantile(zvals, 0.99, na.rm = TRUE)
    if (shift == 0) {
      message(
        "All values are negative on the plotting grid; ",
        "try shift = ", format(-q99, digits = 6), "."
      )
    } else if (near_zero_touch) {
      message("Using shift = ", format(shift, digits = 6), "; duplicate contours merged.")
    } else {
      message("All values negative after applying shift = ", format(shift, digits = 6), ".")
    }
  }
}

mpoly_to_stan <- function(mpoly) {
  p <- get("print.mpoly", asNamespace("mpoly"))
  result <- p(mpoly, stars = TRUE, silent = TRUE, plus_pad = 0L, times_pad = 0L)
  result <- gsub("[*]{2}", "^", result)
  result
}
