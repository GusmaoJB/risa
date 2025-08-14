#' Create polygon of the area of interest
#'
#' Generate an area polygon (convex hull or bounding box) around input points, with optional scaling.
#' The scaling is an affine expansion from the centroid by a factor of (1 + buffer_frac).
#'
#' @param x An `sf` object or data frame with longitude/latitude in the first two columns.
#' @param crs Integer or string; EPSG or proj string for the metric transform. If `NULL`, UTM is chosen automatically.
#' @param area_type One of `"convex_hull"` or `"bbox"`. Default `"convex_hull"`.
#' @param buffer_frac Numeric ≥ 0; fractional expansion relative to the centroid (e.g., 0.5 = +50 percent).
#' @param quiet Logical; suppress informative messages. Default `TRUE`.
#' @return An `sf` POLYGON.
#' @importFrom sf st_crs st_bbox st_sfc st_polygon st_as_sf
#' @examples
#' # Create test data
#' coords <- data.frame(long = c(1,2,2,4), lat = c(4,4,2,2))
#' coords_vec <- df_to_shp(coords)
#' coords_vec_metric <- transform_to_metric(coords_vec)$shape
#'
#' # Create area of interest (default: convex hull)
#' aoi <- create_area(coords)
#'
#' # Plot results
#' plot(sf::st_geometry(aoi), border = "blue", axes=TRUE)
#' plot(sf::st_geometry(coords_vec_metric), add = TRUE, col = "red")
#' @export
create_area <- function(x,
                        crs = NULL,
                        area_type = c("convex_hull", "bbox"),
                        buffer_frac = 0.5,
                        quiet = TRUE) {
  area_type <- match.arg(area_type)
  if (!is.numeric(buffer_frac) || length(buffer_frac) != 1L || buffer_frac < 0)
    stop("`buffer_frac` must be a single non-negative number.")

  # Accept sf or data.frame
  if (inherits(x, "sf")) {
    x_sf <- x
  } else if (inherits(x, "data.frame")) {
    x_sf <- df_to_shp(x)
  } else {
    stop("`x` must be an sf or a data.frame.")
  }

  # Project to metric
  x_m <- transform_to_metric(x_sf, metric_crs = crs, quiet = quiet)
  x_sf     <- x_m$shape
  coords   <- x_m$coordinates
  crs_proj <- sf::st_crs(x_sf)

  # Use XY only (avoid L1/L2/Z/M columns)
  xy <- coords[, 1:2, drop = FALSE]

  # Helper: build polygon from a set of corner points (close the ring)
  make_poly <- function(mat_xy) {
    mat_xy <- as.matrix(mat_xy)
    if (!all(mat_xy[1, ] == mat_xy[nrow(mat_xy), ])) {
      mat_xy <- rbind(mat_xy, mat_xy[1, , drop = FALSE])
    }
    sf::st_as_sf(sf::st_sfc(sf::st_polygon(list(mat_xy)), crs = crs_proj))
  }

  # Helper: scale a set of points about its centroid
  scale_about_centroid <- function(mat_xy, s) {
    cent <- colMeans(mat_xy)
    sweep(mat_xy, 2, cent, function(v, c) (v - c) * s + c)
  }

  scale_factor <- 1 + buffer_frac

  if (area_type == "convex_hull") {
    # Handle degenerate cases (need ≥ 3 unique points)
    uniq <- unique.data.frame(as.data.frame(xy))
    if (nrow(uniq) < 3L) {
      # fall back to a scaled bbox
      if (!quiet) message("Fewer than 3 unique points; using scaled bounding box instead.")
      bb <- sf::st_bbox(x_sf)
      corners <- matrix(
        c(bb$xmin, bb$ymin,
          bb$xmin, bb$ymax,
          bb$xmax, bb$ymax,
          bb$xmax, bb$ymin),
        ncol = 2, byrow = TRUE
      )
      corners <- scale_about_centroid(corners, scale_factor)
      return(make_poly(corners))
    }

    # Compute convex hull on points, then scale the hull vertices
    hull_ix <- chull(uniq[, 1], uniq[, 2])
    hull_xy <- as.matrix(uniq[hull_ix, , drop = FALSE])
    hull_xy <- scale_about_centroid(hull_xy, scale_factor)
    return(make_poly(hull_xy))

  } else { # area_type == "bbox"
    bb <- sf::st_bbox(x_sf)
    corners <- matrix(
      c(bb$xmin, bb$ymin,
        bb$xmin, bb$ymax,
        bb$xmax, bb$ymax,
        bb$xmax, bb$ymin),
      ncol = 2, byrow = TRUE
    )
    corners <- scale_about_centroid(corners, scale_factor)
    return(make_poly(corners))
  }
}
