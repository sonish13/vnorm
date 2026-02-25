## geom_variety-specific helpers
## Keep generic ggplot/contour conversion utilities in helpers-ggplot2.R
## (e.g., xyz_to_isolines(), iso_to_path(), ensure_nonempty_data()).

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

## Sign diagnostics for shifted contour guidance in geom_variety()

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
        "; near-zero values remain on the plotting grid."
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
        "; near-zero values remain on the plotting grid."
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
