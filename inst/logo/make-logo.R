make_vnorm_logo <- function(
    out_dir = "man/figures",
    bg = "#F3F4F6",
    curve_color = "#EAF2FF",
    curve_understroke_color = "#FFFFFF",
    curve_double_stroke = FALSE,
    points_color = "#F59E0B",
    hex_border = "#4FD1C5",
    text_color = "#EAF2FF",
    point_alpha = 0.18,
    point_size = 0.22,
    curve_linewidth = 0.78,
    curve_understroke_linewidth = 1.10,
    curve_n = 701,
    n_points = 2600,
    sample_sd = 0.045,
    sample_w = 1.6,
    dpi = 600,
    show_plot = TRUE,
    text_on_top = TRUE,
    plot_xlim = c(-1.43, 1.43),
    plot_ylim = c(-0.78, 0.78),
    subplot_x = 1,
    subplot_y_top = 0.735,
    subplot_y_bottom = 0.95,
    subplot_width_top = 1.56,
    subplot_width_bottom = 0.96,
    subplot_height_top = 1.06,
    subplot_height_bottom = 0.90,
    text_size = 18
  ) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!requireNamespace("mpoly", quietly = TRUE)) {
    stop("Package 'mpoly' is required.", call. = FALSE)
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  theme_logo <- ggplot2::theme_void() +
    ggplot2::theme(
      # Keep subplot transparent so no rectangular panel bleeds outside the hex.
      plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.position = "none",
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  p <- mpoly::mp("(x^2 + y^2)^2 - 2*(x^2 - y^2)")

  pts <- tryCatch(
    vnorm::rvnorm(
      n = n_points,
      poly = p,
      sd = sample_sd,
      output = "tibble",
      rejection = TRUE,
      w = sample_w,
      verbose = FALSE
    ),
    error = function(e) {
      message("rvnorm point generation failed for logo: ", conditionMessage(e))
      NULL
    }
  )

  layers <- list()
  if (!is.null(pts)) {
    layers <- c(layers, list(
      ggplot2::geom_point(
        data = pts,
        ggplot2::aes(x, y),
        colour = points_color,
        alpha = point_alpha,
        size = point_size,
        shape = 16
      )
    ))
  }

  if (isTRUE(curve_double_stroke)) {
    layers <- c(layers, list(
      vnorm::geom_variety(
        poly = p,
        xlim = plot_xlim,
        ylim = plot_ylim,
        n = curve_n,
        colour = curve_understroke_color,
        linewidth = curve_understroke_linewidth,
        lineend = "round",
        linejoin = "round",
        show.legend = FALSE
      )
    ))
  }

  layers <- c(layers, list(
    vnorm::geom_variety(
      poly = p,
      xlim = plot_xlim,
      ylim = plot_ylim,
      n = curve_n,
      colour = curve_color,
      linewidth = curve_linewidth,
      lineend = "round",
      linejoin = "round",
      show.legend = FALSE
    ),
    ggplot2::coord_equal(),
    theme_logo
  ))

  g <- Reduce(`+`, c(list(ggplot2::ggplot()), layers))

  preview_path <- file.path(out_dir, "logo-preview.png")
  ggplot2::ggsave(
    filename = preview_path,
    plot = g,
    width = 4,
    height = 4.6,
    dpi = dpi,
    bg = bg,
    device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  )

  logo_path <- file.path(out_dir, "logo.png")
  if (requireNamespace("hexSticker", quietly = TRUE)) {
    p_y <- if (isTRUE(text_on_top)) 1.53 else 0.18
    s_y <- if (isTRUE(text_on_top)) subplot_y_top else subplot_y_bottom
    s_h <- if (isTRUE(text_on_top)) subplot_height_top else subplot_height_bottom
    s_w <- if (isTRUE(text_on_top)) subplot_width_top else subplot_width_bottom
    hexSticker::sticker(
      subplot = g,
      package = "vnorm",
      filename = logo_path,
      s_x = subplot_x,
      s_y = s_y,
      s_width = s_w,
      s_height = s_h,
      p_x = 1,
      p_y = p_y,
      p_color = text_color,
      p_family = "sans",
      p_size = if (isTRUE(text_on_top)) text_size else 22,
      h_fill = bg,
      h_color = hex_border,
      dpi = dpi
    )
    message("Saved preview to ", preview_path, " and hex logo to ", logo_path)
  } else {
    message(
      "Saved preview to ", preview_path, ". ",
      "Install 'hexSticker' to generate the final hex logo PNG."
    )
  }

  if (isTRUE(show_plot)) {
    print(g)
  }

  invisible(list(plot = g, preview = preview_path, logo = logo_path))
}

make_vnorm_logo_palette <- function(
    palette = c("teal-gold", "crimson-cyan", "mint-copper"),
    ...
  ) {
  palette <- match.arg(palette)

  args <- switch(
    palette,
    "teal-gold" = list(
      bg = "#2B2F36",
      curve_color = "#EAF2FF",
      points_color = "#F59E0B",
      hex_border = "#4FD1C5",
      text_color = "#EAF2FF"
    ),
    "crimson-cyan" = list(
      bg = "#111827",
      curve_color = "#F9FAFB",
      points_color = "#F43F5E",
      hex_border = "#22D3EE",
      text_color = "#F9FAFB"
    ),
    "mint-copper" = list(
      bg = "#0F172A",
      curve_color = "#E2FDF7",
      points_color = "#C08457",
      hex_border = "#7DD3C7",
      text_color = "#E2FDF7"
    )
  )

  do.call(make_vnorm_logo, c(args, list(...)))
}
