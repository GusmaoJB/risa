#' Create polygon of the area of interest
#'
#' Generate an area polygon (convex hull or bounding box) around input points, with optional buffering.
#'
#' @param x An `sf` object or data frame with longitude/latitude in the first two columns.
#' @param crs Integer or string; EPSG code or proj4 string for metric transformation. If `NULL`, select UTM zone automatically.
#' @param area_type Character; one of \code{"convex_hull"} or \code{"bbox"}. Default is \code{"convex_hull"}.
#' @param buffer_frac Numeric; fraction by which to buffer the hull or bounding box. (e.g. 0.5 for 50\%.)
#' @return An `sf` object representing the area polygon.
#' @importFrom sf st_crs st_coordinates st_centroid st_transform st_convex_hull st_union st_polygon st_sfc st_geometry
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
#' plot(st_geometry(aoi), border = "blue", axes=TRUE)
#' plot(st_geometry(coords_vec_metric), add = TRUE, col = "red")
#' @export
create_area <- function(x, crs = NULL, area_type = c("convex_hull", "bbox"), buffer_frac = 0.5) {
  area_type <- match.arg(area_type)

  # Check input
  if (inherits(x, "sf")) {
    x_sf <- x
  } else if (inherits(x, "data.frame")) {
    x_sf <- df_to_shp(x)
  } else {
    stop("`x` must be an sf or data.frame")
  }

  # Project to metric
  x_m <- transform_to_metric(x_sf, metric_crs = crs)
  x_sf <- x_m$shape
  coords <- x_m$coordinates
  crs_proj <- st_crs(x_sf)

  # Convex hull case (default)
  if (area_type == "convex_hull") {
    # Original centroid of points
    cent <- colMeans(coords)
    # Scale factor
    scale_factor <- 1 + buffer_frac
    # Scale coords away from centroid, then re-center
    scaled_coords <- sweep(coords, 2, cent, FUN = function(pt, c) (pt - c) * scale_factor + c)
    # Find convex hull on the scaled points
    hull_ix      <- chull(scaled_coords)
    hull_coords  <- scaled_coords[c(hull_ix, hull_ix[1]), , drop = FALSE]
    # Build sf polygon and return
    poly     <- st_polygon(list(hull_coords))
    out_sfc  <- st_sfc(poly, crs = crs_proj)
    return(st_as_sf(out_sfc))

  } else {
    # Bounding box case
    # Compute original centroid
    cent <- colMeans(coords)
    # Scale factor
    scale_factor <- 1 + buffer_frac
    # Extract original bbox
    bb <- st_bbox(x_sf)
    corners <- matrix(c(
      bb$xmin, bb$ymin,
      bb$xmin, bb$ymax,
      bb$xmax, bb$ymax,
      bb$xmax, bb$ymin
    ), ncol = 2, byrow = TRUE)
    # Scale corners from centroid and re-center
    scaled_corners <- sweep(corners, 2, cent, FUN = function(pt, c) (pt - c) * scale_factor + c)
    # Close the polygon ring
    poly_coords <- rbind(scaled_corners, scaled_corners[1, ])
    poly <- st_polygon(list(poly_coords))
    out_sfc <- st_sfc(poly, crs = crs_proj)
    return(st_as_sf(out_sfc))
  }
}
