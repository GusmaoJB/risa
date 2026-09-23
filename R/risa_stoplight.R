#' Add stoplight data-quality indicators to facet panels
#'
#' `risa_stoplight()` adds a small "traffic light" housing to a corner of each
#' panel of a ggplot. Each light is a filled circle whose colour is taken from
#' the data, so it can show the quality or status of several data sources for
#' the panel (for example, "good", "caution" or "poor" data). The layer is
#' drawn in panel-relative coordinates, so it does not affect the axis scales
#' or the position of any other layer.
#'
#' @details
#' The layer is drawn once per panel. `data` is split across panels through
#' any faceting columns it shares with the plot (for example `species` and
#' `stressor` when using `facet_grid(stressor ~ species)`), and each panel
#' gets its own set of lights. A data frame without faceting columns draws
#' the same lights in every panel.
#'
#' Lights are matched to `light_order` by name, ignoring case. If a level in
#' `light_order` has no matching row in a panel, that circle is filled with
#' `na_col`. If a panel has several rows for the same light, only the first
#' is used. Lights are stacked top to bottom in the order of `light_order`.
#'
#' The `location` argument sets the corner where the housing is anchored, and
#' `pad_x` and `pad_y` set its distance from the panel edges. Labels (when
#' `show_labels = TRUE`) are placed to the right of the housing, so they
#' are best used with left-hand locations (`"tl"` and `"bl"`). With `"tr"` or
#' `"br"` the labels may extend beyond the panel edge.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()]. It
#'   must map `light` and `value`. If `NULL` (the default), it uses
#'   `aes(light = .data$light, value = .data$value)`.
#' @param data A data frame with one row per light and panel. It must contain
#'   a character column `light` (the name of each light) and a character
#'   column `value` (a valid R colour, such as `"green"`, `"yellow"` or
#'   `"#FF0000"`, used as the fill of that light). Include any faceting
#'   columns needed to place lights in specific panels. This argument is
#'   required.
#' @param light_order Character vector of unique light names, giving the order
#'   in which lights are drawn from top to bottom. If `NULL` (the default),
#'   the order in which lights first appear in `data` is used, and a message
#'   is printed.
#' @param location Corner of each panel where the housing is placed: `"tl"`
#'   (top left), `"tr"` (top right), `"bl"` (bottom left) or `"br"` (bottom
#'   right, the default).
#' @param circle_diameter A [grid::unit()] of length 1 giving the diameter of
#'   each light. Defaults to `grid::unit(0.5, "cm")`.
#' @param spacing A [grid::unit()] of length 1 giving the vertical gap between
#'   adjacent lights. Defaults to `grid::unit(0.15, "cm")`.
#' @param housing_pad A [grid::unit()] of length 1 giving the padding between
#'   the housing edge and the circles. Defaults to `grid::unit(0.2, "cm")`.
#' @param pad_x,pad_y [grid::unit()]s of length 1 giving the horizontal and
#'   vertical distance between the housing and the nearest panel edges.
#'   Both default to `grid::unit(0.25, "cm")`.
#' @param housing_fill Fill colour of the housing. Defaults to `"grey20"`.
#' @param housing_col Outline colour of the housing. Defaults to `"black"`.
#' @param circle_col Outline colour of the circles. Defaults to `"black"`.
#' @param line_width Line width of the housing outline. The circle outlines use
#'   0.75 times this value. Defaults to `1`.
#' @param na_col Fill colour for lights in `light_order` that have no matching
#'   row in a panel's data. Defaults to `"grey90"`.
#' @param show_labels Logical. If `TRUE`, the light names are drawn as text
#'   labels next to each circle. Defaults to `FALSE`.
#' @param text_col Colour of the labels. Defaults to `"black"`.
#' @param text_cex Numeric character expansion factor of the labels. Defaults
#'   to `0.7`.
#' @param text_face Font face of the labels (for example `"bold"` or
#'   `"italic"`). If `NULL` (the default), the graphics device default is used.
#' @param text_family Font family of the labels. The default `""` uses the
#'   graphics device default.
#' @param text_pad A [grid::unit()] of length 1 giving the gap between the
#'   housing and the labels. Defaults to `grid::unit(0.15, "cm")`.
#'
#' @return A ggplot2 layer that can be added to a plot with `+`.
#'
#' @section Aesthetics:
#' `risa_stoplight()` understands the following aesthetics (required
#' aesthetics are in bold):
#'
#' * **`light`**: name of the light (for example, `"Species data"`).
#' * **`value`**: fill colour of the light.
#' * `location`: corner of the panel. It is set through the `location`
#'   argument and should not be mapped to a column of the data.
#'
#' @seealso [ggplot2::facet_grid()], [ggplot2::facet_wrap()],
#'   [grid::unit()]
#'
#' @examples
#' library(ggplot2)
#'
#' # Creating a mock data-quality data frame
#' stoplight_df <- data.frame(
#'   species = rep(c("Dolphin", "Turtle"), each = 6),
#'   stressor = rep(rep(c("Trawling", "Gillnet"), each=3), 2),
#'   light = rep(c("Species data", "Stressor data", "Bycatch data"), 4),
#'   value = c("green", "green", "yellow",
#'             "red", "green", "red",
#'             "green", "red", "yellow",
#'             "yellow", "yellow", "yellow"),
#'   stringsAsFactors = FALSE)
#'
#' # Creating a mock data.frame
#' base_df <- data.frame(
#'   species = rep(c("Dolphin", "Turtle"), each = 20),
#'   stressor = rep(c("Trawling", "Gillnet"), 20),
#'   x = c(rnorm(20, 0, 1), rnorm(20, 5, 1)),
#'   y = c(rnorm(20, 0, 1), rnorm(20, 5, 1))
#' )
#'
#' # Plot with stoplight layer
#' ggplot(base_df, aes(x, y)) +
#'   geom_point(alpha = 0.5) +
#'   facet_grid(stressor ~ species) +
#'   risa_stoplight(
#'     location = "tl",
#'     show_labels = TRUE,
#'     data = stoplight_df,
#'     light_order = c("Species data", "Stressor data", "Bycatch data"),
#'     circle_diameter = grid::unit(0.45, "cm")
#'   ) +
#'   theme_bw()
#'
#' @export
risa_stoplight <- function(mapping = NULL,
                           data = NULL,
                           light_order = NULL,
                           location = "br",
                           circle_diameter = grid::unit(0.5, "cm"),
                           spacing = grid::unit(0.15, "cm"),
                           housing_pad = grid::unit(0.2, "cm"),
                           pad_x = grid::unit(0.25, "cm"),
                           pad_y = grid::unit(0.25, "cm"),
                           housing_fill = "grey20",
                           housing_col = "black",
                           circle_col = "black",
                           line_width = 1,
                           na_col = "grey90",
                           show_labels = FALSE,
                           text_col = "black",
                           text_cex = 0.7,
                           text_face = NULL,
                           text_family = "",
                           text_pad = grid::unit(0.15, "cm")) {
  
  if (is.null(data)) {
    stop(
      "`risa_stoplight()` requires a `data` data.frame with columns ",
      "`light` and `value` (plus any faceting column).",
      call. = FALSE
    )
  }
  
  if (!location %in% c("tl", "tr", "bl", "br")) {
    stop(
      "`location` must be one of 'tl', 'tr', 'bl', or 'br'.",
      call. = FALSE
    )
  }
  
  if (!all(c("light", "value") %in% names(data))) {
    stop("`data` must contain columns named `light` and `value`.", call. = FALSE)
  }
  
  if (is.null(mapping)) {
    mapping <- ggplot2::aes(
      light = .data$light,
      value = .data$value
    )
  }
  
  mapping$location <- location
  
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = ggplot2::StatIdentity,
    geom = GeomStoplight,
    position = ggplot2::PositionIdentity,
    show.legend = FALSE,
    inherit.aes = FALSE,
    params = list(
      light_order = light_order,
      circle_diameter = circle_diameter,
      spacing = spacing,
      housing_pad = housing_pad,
      pad_x = pad_x,
      pad_y = pad_y,
      housing_fill = housing_fill,
      housing_col = housing_col,
      circle_col = circle_col,
      line_width = line_width,
      na_col = na_col,
      show_labels = show_labels,
      text_col = text_col,
      text_cex = text_cex,
      text_face = text_face,
      text_family = text_family,
      text_pad = text_pad
    )
  )
}

#' @rdname risa_stoplight
#' @format NULL
#' @usage NULL
#' @export
GeomStoplight <- ggplot2::ggproto(
  "GeomStoplight",
  ggplot2::Geom,
  
  extra_params = "",
  
  required_aes = c("light", "value"),
  
  default_aes = ggplot2::aes(location = "br"),
  
  handle_na = function(data, params) {
    data
  },
  
  draw_key = ggplot2::draw_key_blank,
  
  draw_panel = function(data, panel_params, coordinates,
                        light_order = NULL,
                        circle_diameter = grid::unit(0.5, "cm"),
                        spacing = grid::unit(0.15, "cm"),
                        housing_pad = grid::unit(0.2, "cm"),
                        pad_x = grid::unit(0.25, "cm"),
                        pad_y = grid::unit(0.25, "cm"),
                        housing_fill = "grey20",
                        housing_col = "black",
                        circle_col = "black",
                        line_width = 1,
                        na_col = "grey90",
                        show_labels = FALSE,
                        text_col = "black",
                        text_cex = 0.7,
                        text_face = NULL,
                        text_family = "",
                        text_pad = grid::unit(0.15, "cm")) {
    
    if (nrow(data) == 0) {
      return(grid::nullGrob())
    }
    
    if (is.null(light_order)) {
      message("User did not inform light order. Assuming light orders from data input.")
      light_order <- unique(data$light)
    } else if (is.character(light_order) &
               length(light_order) >= 1) {
      message("Using user-informed light order.")
      if (anyDuplicated(light_order)) {
        stop("`light_order` must contain unique light names.", call. = FALSE)
      }
    } else {
      stop("Error: Light order input must be a character vector.")
    }
    
    stopifnot(
      grid::is.unit(circle_diameter), length(circle_diameter) == 1,
      grid::is.unit(spacing), length(spacing) == 1,
      grid::is.unit(housing_pad), length(housing_pad) == 1,
      grid::is.unit(pad_x), length(pad_x) == 1,
      grid::is.unit(pad_y), length(pad_y) == 1
    )
    
    location <- data$location[1]
    adj_x <- as.numeric(grepl("r", location))
    adj_y <- as.numeric(grepl("t", location))
    
    # match a fill color to each level in light_order (case-insensitive),
    # falling back to na_col for any level with no matching row in this panel
    light_vals <- tolower(as.character(data$light))
    fills <- vapply(tolower(light_order), function(lvl) {
      idx <- which(light_vals == lvl)
      if (length(idx) == 0) na_col else as.character(data$value[idx[1]])
    }, character(1))
    
    n <- length(light_order)
    gap <- if (n > 1) spacing else grid::unit(0, "cm")
    
    box_width  <- circle_diameter + housing_pad * 2
    box_height <- circle_diameter * n + gap * (n - 1) + housing_pad * 2
    
    origin_x <- grid::unit(adj_x, "npc") - adj_x * box_width + (0.5 - adj_x) * 2 * pad_x
    origin_y <- grid::unit(adj_y, "npc") - adj_y * box_height + (0.5 - adj_y) * 2 * pad_y
    
    housing_grob <- grid::roundrectGrob(
      x = origin_x, y = origin_y,
      width = box_width, height = box_height,
      just = c("left", "bottom"),
      r = grid::unit(0.15, "snpc"),
      gp = grid::gpar(fill = housing_fill, col = housing_col, lwd = line_width)
    )
    
    idx <- seq_len(n)
    cx <- rep(origin_x + housing_pad + circle_diameter / 2, n)
    cy <- origin_y + box_height - housing_pad - circle_diameter / 2 -
      (idx - 1) * (circle_diameter + gap)
    
    circle_grobs <- grid::circleGrob(
      x = cx, y = cy, r = circle_diameter / 2,
      gp = grid::gpar(fill = fills, col = circle_col, lwd = line_width * 0.75)
    )
    
    label_grob <- NULL
    if (show_labels) {
      lx <- rep(origin_x + housing_pad * 2 + circle_diameter + text_pad, n)
      ly <- cy
      just <- c("left", "center")
      label_grob <- grid::textGrob(
        label = light_order,
        x = lx, y = ly,
        just = just,
        gp = grid::gpar(
          col = text_col, cex = text_cex,
          fontface = text_face, fontfamily = text_family
        )
      )
    }
    
    grid::gList(housing_grob, circle_grobs, label_grob)
  }
)