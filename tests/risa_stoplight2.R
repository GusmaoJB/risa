risa_stoplight2 <- function(mapping = NULL,
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
      "`light` and `value` (plus any faceting column, e.g. `city`).",
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


# Example
library(ggplot2)

# Creating a mock data-quality data frame
stoplight_df <- data.frame(
  species = rep(c("Dolphin", "Turtle"), each = 6),
  stressor = rep(rep(c("Trawling", "Gillnet"), each=3), 2),
  light = rep(c("Species data", "Stressor data", "Bycatch data"), 4),
  value = c("green", "green", "yellow",
            "red", "green", "red",
            "green", "red", "yellow",
            "yellow", "yellow", "yellow"),
  stringsAsFactors = FALSE)

# Creating a mock data.frame
base_df <- data.frame(
  species = rep(c("Dolphin", "Turtle"), each = 20),
  stressor = rep(c("Trawling", "Gillnet"), 20),
  x = c(rnorm(20, 0, 1), rnorm(20, 5, 1)),
  y = c(rnorm(20, 0, 1), rnorm(20, 5, 1))
)

ggplot(base_df, aes(x, y)) +
  geom_point(alpha = 0.5) +
  facet_grid(stressor ~ species) +
  risa_stoplight2(
    location = "tl",
    show_labels = TRUE,
    data = stoplight_df,
    light_order = c("Species data", "Stressor data", "Bycatch data"),
    circle_diameter = grid::unit(0.45, "cm")
  ) +
  theme_bw()
