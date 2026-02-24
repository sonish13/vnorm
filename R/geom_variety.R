#' One-Dimensional Varieties in Two Dimensions
#'
#' Plot implicit polynomial varieties with `ggplot2`.
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
      check_sign_warning(df0$z, shift)
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
        check_sign_warning(df0$z, shift)
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

poly_to_df <- function(poly, xlim, ylim, nx, ny, shift = 0) {
  # Evaluate polynomial on a regular x/y grid for contour extraction.
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
  # Heuristic: many tiny groups indicates under-resolved contours.
  if (nrow(df) == 0) return(FALSE)
  n_per_group <- as.numeric(table(df$group))
  if (length(n_per_group) <= 8) return(FALSE)

  frac_short <- mean(n_per_group < 12)
  frac_short > 0.7
}

mean_nn_distance <- function(df_a, df_b) {
  # Mean nearest-neighbor distance from A to B.
  dx <- outer(df_a$x, df_b$x, "-")
  dy <- outer(df_a$y, df_b$y, "-")
  d2 <- dx * dx + dy * dy
  mean(sqrt(apply(d2, 1, min)))
}

path_proximity <- function(df_a, df_b) {
  # Symmetric path proximity metric based on NN distances.
  max(mean_nn_distance(df_a, df_b), mean_nn_distance(df_b, df_a))
}

collapse_near_duplicate_contours <- function(df, tol) {
  # Merge contours that are geometric near-duplicates (common after shift).
  if (nrow(df) == 0) return(df)
  dup_mult <- 6

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
    if (d_cent <= dup_mult * tol && d_path <= dup_mult * tol) {
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
        if (d_cent > dup_mult * tol) next

        d_path <- path_proximity(groups[[g]], groups[[h]])
        if (d_path <= dup_mult * tol) {
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

  # Post-merge cleanup: drop tiny fragments that hug a dominant contour.
  out_groups <- split(out, out$group)
  if (length(out_groups) >= 2) {
    sizes <- vapply(out_groups, nrow, integer(1))
    g_main <- names(out_groups)[which.max(sizes)]
    main_df <- out_groups[[g_main]]
    main_n <- nrow(main_df)

    drop_groups <- character(0)
    for (g in names(out_groups)) {
      if (identical(g, g_main)) next
      frag_df <- out_groups[[g]]
      frag_n <- nrow(frag_df)
      if (frag_n > max(25L, as.integer(0.8 * main_n))) next

      # One-sided distance is the right test here: a short fragment can lie on
      # top of the main path while the symmetric distance stays large.
      d_path_frag <- mean_nn_distance(frag_df, main_df)
      if (d_path_frag <= 20 * tol) {
        drop_groups <- c(drop_groups, g)
      }
    }

    if (length(drop_groups) > 0) {
      out <- out[!(out$group %in% drop_groups), , drop = FALSE]
      out$group <- factor(out$group)
      out$piece <- as.integer(out$group)
    }
  }
  out
}

variety_paths_with_refinement <- function(
    poly, rangex, rangey, nx, ny, shift, group
  ) {
  # Retry at higher grid resolution if initial contours are fragmented.
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

snap_shifted_contours_to_variety <- function(df, poly) {
  # For no-sign-change polynomials (e.g., squared factors), shift generates a
  # nearby level set. Snap those points back to poly = 0 to recover topology.
  if (nrow(df) == 0) return(df)
  if (!all(c("x", "y") %in% names(df))) return(df)

  proj <- tryCatch({
    varorder <- c("x", "y")
    gfunc <- as.function(poly, varorder = varorder, silent = TRUE)
    dg <- stats::deriv(poly, var = varorder)
    dgfunc <- as.function(dg, varorder = varorder, silent = TRUE)

    x <- as.matrix(df[, c("x", "y"), drop = FALSE])

    # Typical contour-point spacing for step capping.
    dseg <- sqrt(diff(x[, 1])^2 + diff(x[, 2])^2)
    dseg <- dseg[is.finite(dseg) & dseg > 0]
    base_step <- if (length(dseg) > 0) stats::median(dseg) else 0.05
    if (!is.finite(base_step) || base_step <= 0) base_step <- 0.05
    # Squared-polynomial correction is usually underpowered near singularities;
    # allow larger but still bounded steps.
    max_step <- 12 * base_step
    snap_gains <- c(0.5, 1, 2, 4)

    for (iter in seq_len(12L)) {
      vals <- as.numeric(gfunc(x))
      if (!all(is.finite(vals))) break

      grads <- t(vapply(
        seq_len(nrow(x)),
        function(i) as.numeric(dgfunc(x[i, ])),
        numeric(2)
      ))
      g2 <- rowSums(grads * grads)
      g2_fin <- g2[is.finite(g2)]
      g2_scale <- if (length(g2_fin) > 0) max(1, max(g2_fin)) else 1
      g2_tol <- 1e3 * .Machine$double.eps * g2_scale
      ok <- is.finite(g2) & g2 > g2_tol
      if (!any(ok)) break

      step0 <- matrix(0, nrow = nrow(x), ncol = 2)
      step0[ok, ] <- (vals[ok] / (g2[ok] + g2_tol)) * grads[ok, , drop = FALSE]

      best_x <- x
      best_abs <- abs(vals)
      best_abs[!is.finite(best_abs)] <- Inf

      for (gain in snap_gains) {
        step <- gain * step0

        # Prevent jumping across branches near singular crossings.
        step_norm <- sqrt(rowSums(step * step))
        too_big <- is.finite(step_norm) & step_norm > max_step
        if (any(too_big)) {
          step[too_big, ] <- step[too_big, , drop = FALSE] *
            (max_step / step_norm[too_big])
        }

        x_trial <- x - step
        vals_trial <- as.numeric(gfunc(x_trial))
        abs_trial <- abs(vals_trial)
        abs_trial[!is.finite(abs_trial)] <- Inf

        improve <- ok & (abs_trial < best_abs)
        if (any(improve)) {
          best_x[improve, ] <- x_trial[improve, , drop = FALSE]
          best_abs[improve] <- abs_trial[improve]
        }
      }

      x <- best_x

      best_fin <- best_abs[is.finite(best_abs)]
      if (length(best_fin) == 0) break
      max_abs_new <- max(best_fin)
      if (is.finite(max_abs_new) && max_abs_new <= 1e-8) break
    }

    x
  }, error = function(e) NULL)

  if (is.null(proj)) return(df)
  proj <- as.matrix(proj)
  if (!all(dim(proj) == c(nrow(df), 2L))) return(df)
  bad <- !is.finite(proj[, 1]) | !is.finite(proj[, 2])
  if (any(bad)) {
    proj[bad, ] <- as.matrix(df[bad, c("x", "y"), drop = FALSE])
  }

  df$x <- proj[, 1]
  df$y <- proj[, 2]
  split_large_projected_jumps(df)
}

split_large_projected_jumps <- function(df) {
  # Projection can collapse different shifted components onto the same variety,
  # creating jumps within a path. Split those jumps before plotting.
  if (nrow(df) == 0 || !"group" %in% names(df)) return(df)

  groups <- split(df, df$group)
  if (length(groups) == 0) return(df)

  out <- lapply(seq_along(groups), function(i) {
    g <- groups[[i]]
    if (nrow(g) <= 2) return(g)
    dx <- diff(g$x)
    dy <- diff(g$y)
    d <- sqrt(dx * dx + dy * dy)
    d_pos <- d[is.finite(d) & d > 0]
    if (length(d_pos) == 0) return(g)
    base <- stats::median(d_pos, na.rm = TRUE)
    if (!is.finite(base) || base <= 0) return(g)
    cut_idx <- which(d > 8 * base)
    if (length(cut_idx) == 0) return(g)

    seg_id <- cumsum(c(TRUE, seq_len(nrow(g) - 1) %in% cut_idx))
    g$group <- factor(paste0(as.character(g$group[1]), "_s", seg_id))
    g
  })

  out <- dplyr::bind_rows(out)
  out$group <- factor(out$group)
  out$piece <- as.integer(out$group)
  out
}

has_strict_sign_change <- function(zvals) {
  # True only when the field takes both positive and negative values,
  # ignoring near-zero floating-point noise.
  z <- zvals[is.finite(zvals)]
  if (length(z) == 0) return(FALSE)
  scale_z <- max(abs(z), na.rm = TRUE)
  tol <- 100 * sqrt(.Machine$double.eps) * max(1, scale_z)
  any(z > tol, na.rm = TRUE) && any(z < -tol, na.rm = TRUE)
}

check_sign_warning <- function(zvals, shift) {
  # Emit user guidance when grid values do not cross zero.
  z <- zvals[is.finite(zvals)]
  if (length(z) == 0) return()
  scale_z <- max(abs(z), na.rm = TRUE)
  tol <- sqrt(.Machine$double.eps) * max(1, scale_z)
  signs <- ifelse(z > tol, 1L, ifelse(z < -tol, -1L, 0L))
  signs <- signs[signs != 0]
  if (length(signs) == 0) return()
  near_zero_touch <- min(abs(z), na.rm = TRUE) <= tol

  if (all(signs == 1)) {
    z_pos <- z[z > tol]
    q_small <- if (length(z_pos) > 0) {
      stats::quantile(z_pos, 0.01, na.rm = TRUE)
    } else {
      tol
    }
    q_small <- max(as.numeric(q_small), tol)
    if (shift == 0) {
      message(
        "All values are positive on the plotting grid; ",
        "try shift = ", format(-q_small, digits = 6), "."
      )
    } else if (near_zero_touch) {
      message(
        "Using shift = ",
        format(shift, digits = 6),
        "; duplicate contours merged."
      )
    } else {
      message(
        "All values positive after applying shift = ",
        format(shift, digits = 6),
        "."
      )
    }
  } else if (all(signs == -1)) {
    z_neg <- z[z < -tol]
    q_large <- if (length(z_neg) > 0) {
      stats::quantile(z_neg, 0.99, na.rm = TRUE)
    } else {
      -tol
    }
    q_large <- min(as.numeric(q_large), -tol)
    if (shift == 0) {
      message(
        "All values are negative on the plotting grid; ",
        "try shift = ", format(-q_large, digits = 6), "."
      )
    } else if (near_zero_touch) {
      message(
        "Using shift = ",
        format(shift, digits = 6),
        "; duplicate contours merged."
      )
    } else {
      message(
        "All values negative after applying shift = ",
        format(shift, digits = 6),
        "."
      )
    }
  }
}

mpoly_to_stan <- function(mpoly) {
  # Internal pretty-printer to generate parseable legend labels.
  p <- get("print.mpoly", asNamespace("mpoly"))
  result <- p(mpoly, stars = TRUE, silent = TRUE, plus_pad = 0L, times_pad = 0L)
  result <- gsub("[*]{2}", "^", result)
  result
}
