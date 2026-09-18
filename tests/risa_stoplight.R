#' Add a stoplight (traffic light) annotation to a ggplot, one per facet panel
#'
#' Draws a small "housing" with a stack of colored circles ("lights"), sized
#' and positioned in plot (npc) coordinates so it behaves like a scale bar or
#' north arrow: it stays a fixed size and corner position regardless of the
#' underlying data coordinate system (works fine with sf/coord_sf maps).
#'
#' Unlike `annotation_scale()`, which draws the *same* thing in every panel
#' from a dummy one-row data.frame, this geom is designed to draw *different*
#' colors in each facet panel. To do that, `data` must contain:
#'   - `light` : the light's position/level, matched against `light_order`
#'   - `value` : the fill color for that light (any valid R color spec)
#'   - one column per faceting variable you use (e.g. `city`), so that
#'     ggplot2 can route each row to the correct panel. This routing happens
#'     automatically because ggplot2 matches facet variables against any
#'     column present in a layer's raw `data`, whether or not it is mapped
#'     via `aes()`.
#'
#' @param mapping Defaults to `aes(light = light, value = value)`, which is
#'   almost always what you want. Override only if your columns are named
#'   differently.
#' @param data A data.frame with (at least) `light` and `value` columns, plus
#'   any faceting column(s) used elsewhere in the plot.
#' @param light_order Character vector giving the levels of `light` from
#'   top/first to bottom/last (vertical) or left to right (horizontal).
#'   Levels present in `light_order` but missing from `data` for a given
#'   panel are drawn using `na_col`.
#' @param orientation "vertical" (default, classic traffic-light column) or
#'   "horizontal".
#' @param circle_diameter,spacing,housing_pad Sizes (grid units) of each
#'   light circle, the gap between circles, and the padding between the
#'   circles and the housing's edge.
#' @param pad_x,pad_y Padding (grid units) between the housing and the panel
#'   edge.
#' @param housing_fill,housing_col Fill/outline color of the housing box.
#' @param circle_col Outline color of each light circle.
#' @param line_width Line width used for both housing and circle outlines.
#' @param na_col Fill color used for a light level with no matching row.
#' @param show_labels Whether to draw a text label (from `light_order`) next
#'   to each circle.
#' @param text_col,text_cex,text_face,text_family,text_pad Label styling.
#'
#' @export
risa_stoplight <- function(mapping = NULL,
                                 data = NULL,
                                 ...,
                                 light_order = c("high", "medium", "low"),
                                 orientation = c("vertical", "horizontal"),
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
      "`light` and `value` (plus any faceting column, e.g. `city`).",
      call. = FALSE
    )
  }
  if (!all(c("light", "value") %in% names(data))) {
    stop("`data` must contain columns named `light` and `value`.", call. = FALSE)
  }
  
  # Default mapping: pull `light`/`value` through as aesthetics so they
  # survive into the built layer data that draw_panel() receives. Facet
  # routing (e.g. by `city`) does NOT depend on this -- it happens against
  # the raw `data` argument regardless of what's mapped in `aes()`.
  if (is.null(mapping)) {
    mapping <- ggplot2::aes(light = .data$light, value = .data$value)
  }
  
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = ggplot2::StatIdentity,
    geom = GeomStoplight,
    position = ggplot2::PositionIdentity,
    show.legend = FALSE,
    inherit.aes = FALSE,
    params = list(
      ...,
      light_order = light_order,
      orientation = match.arg(orientation),
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
#' @export
GeomStoplight <- ggplot2::ggproto(
  "GeomStoplight",
  ggplot2::Geom,
  
  extra_params = "",
  
  # `location` works the same way as in annotation_scale(): it's not a real
  # column in the input data, it's supplied purely via default_aes so each
  # panel/row gets "br" unless the user maps/supplies something else.
  default_aes = ggplot2::aes(location = "br"),
  
  handle_na = function(data, params) {
    data
  },
  
  draw_key = ggplot2::draw_key_blank,
  
  draw_panel = function(data, panel_params, coordinates,
                        light_order = c("high", "medium", "low"),
                        orientation = "vertical",
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
                        show_labels = TRUE,
                        text_col = "black",
                        text_cex = 0.7,
                        text_face = NULL,
                        text_family = "",
                        text_pad = grid::unit(0.15, "cm")) {
    
    if (nrow(data) == 0) {
      return(grid::nullGrob())
    }
    
    stopifnot(
      is.character(light_order), length(light_order) >= 1,
      orientation %in% c("vertical", "horizontal"),
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
    
    if (orientation == "vertical") {
      box_width  <- circle_diameter + housing_pad * 2
      box_height <- circle_diameter * n + gap * (n - 1) + housing_pad * 2
    } else {
      box_width  <- circle_diameter * n + gap * (n - 1) + housing_pad * 2
      box_height <- circle_diameter + housing_pad * 2
    }
    
    # same corner-anchoring trick as annotation_scale(): pick which corner
    # of the housing sits at the npc origin, then pad inward from there
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
    if (orientation == "vertical") {
      cx <- rep(origin_x + housing_pad + circle_diameter / 2, n)
      cy <- origin_y + box_height - housing_pad - circle_diameter / 2 -
        (idx - 1) * (circle_diameter + gap)
    } else {
      cx <- origin_x + housing_pad + circle_diameter / 2 +
        (idx - 1) * (circle_diameter + gap)
      cy <- rep(origin_y + housing_pad + circle_diameter / 2, n)
    }
    
    circle_grobs <- grid::circleGrob(
      x = cx, y = cy, r = circle_diameter / 2,
      gp = grid::gpar(fill = fills, col = circle_col, lwd = line_width * 0.75)
    )
    
    label_grob <- NULL
    if (show_labels) {
      if (orientation == "vertical") {
        lx <- rep(origin_x + housing_pad * 2 + circle_diameter + text_pad, n)
        ly <- cy
        just <- c("left", "center")
      } else {
        lx <- cx
        ly <- rep(origin_y - text_pad, n)
        just <- c("center", "top")
      }
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


# ---------------------------------------------------------------------------
# Example usage, matching the Santos / Ubatuba case from the prompt
# ---------------------------------------------------------------------------
library(ggplot2)

stoplight_df <- data.frame(
  city  = c("Santos", "Santos", "Santos", "Santos", "Ubatuba", "Ubatuba", "Ubatuba", "Ubatuba"),
  light = c("low", "medium", "high", "very high", "low", "medium", "high", "very high"),
  value = c("green", "green", "yellow", "red", "green", "red", "yellow", "red"),
  stringsAsFactors = FALSE
)

# any per-panel base plot -- works the same with sf/coord_sf layers,
# since the stoplight is drawn in npc space, not data space
base_df <- data.frame(
  city = rep(c("Santos", "Ubatuba"), each = 20),
  x = c(rnorm(20, 0, 1), rnorm(20, 5, 1)),
  y = c(rnorm(20, 0, 1), rnorm(20, 5, 1))
)

ggplot(base_df, aes(x, y)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~city) +
  risa_stoplight(
    data = stoplight_df,
    light_order = c("very high", "high", "medium", "low"),
    circle_diameter = grid::unit(0.45, "cm")
  ) +
  theme_bw()


if (interactive()) {
  
  
}
