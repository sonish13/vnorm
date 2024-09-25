library(mpoly)
library(tidyverse)

poly_to_df <- function(
    poly,
    xlim = c(-1,1),
    ylim = c(-1,1),
    n = 101,
    nx = n,
    ny = n) {
  if (!is.mpoly(poly)) poly <- mpoly::mp(poly)
  f <- mpoly:::as.function.mpoly(x = poly, varorder = c("x", "y"),
                                 vector = FALSE, silent = TRUE)
  df <- tibble(
    "x" = seq(xlim[1], xlim[2], length.out = nx),
    "y" = seq(ylim[1], ylim[2], length.out = ny)
  ) %>%
    base::expand.grid() %>%
    dplyr::mutate("z" = f(x, y))
  df
}

StatVariety <- ggproto("StatVariety", StatContour,
                       setup_data = function(data, params) {
                         if (!is.mpoly(poly)) poly <- mpoly::mp(poly)
                         f <- mpoly:::as.function.mpoly(x = poly, varorder = c("x", "y"),
                                                        vector = FALSE, silent = TRUE)
                         df <- tibble(
                           "x" = seq(xlim[1], xlim[2], length.out = nx),
                           "y" = seq(ylim[1], ylim[2], length.out = ny)
                         ) %>%
                           base::expand.grid() %>%
                           dplyr::mutate("z" = f(x, y))
                         df
                       }
)

stat_variety <- function(mapping = NULL, data = NULL,
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
                         inherit.aes = TRUE) {
  if (is.null(data)) {
    data <- ggplot2:::ensure_nonempty_data
  }
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

GeomVariety <- ggproto("GoemVariety", GeomContour,

                       default_aes = aes(
                         weight = 1,
                         colour = "#123456",
                         linewidth = 0.5,
                         linetype = 1,
                         alpha = NA,
                         breaks = 0
                       )
)

geom_variety <- function(mapping = NULL, data = NULL,
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
                         inherit.aes = TRUE) {
  if (is.null(data)) {
    data <- ggplot2:::ensure_nonempty_data
  }
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

p <- mp("x^2 + y^2 - 1")
ggplot() +
  geom_variety(poly = p)

dft <- poly_to_df(p)
dft %>% ggplot(aes(x , y , z = z )) +
  geom_contour(breaks = 0)
