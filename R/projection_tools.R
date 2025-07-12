#' Guess coordinate reference system from geographic bounds
#'
#' If coordinates fall within longitude [-180,180] and latitude [-90,90], returns EPSG:4326; else `NA`.
#'
#' @param shp An `sf` object or data frame with longitude/latitude in first two columns.
#' @return Integer EPSG code (4326) or `NA` if bounds are invalid.
#' @importFrom sf st_coordinates
#' @examples
#' # Create test data
#' coords <- data.frame(long = c(1,2,2,4), lat = c(4,4,2,2))
#'
#' # Gussing coordinates
#' guess_crs(coords)
#'
#' @export
guess_crs <- function(shp) {
  # Extract coordinates depending on input type
  if (inherits(shp, "sf")) {
    coords <- sf::st_coordinates(shp)
  } else if (is.data.frame(shp)) {
    if (ncol(shp) < 2) {
      stop("Data frame must have at least two columns: longitude (1st) and latitude (2nd).")
    }
    coords <- as.matrix(shp[, 1:2])
    if (!is.numeric(coords)) {
      stop("The first two columns of the data frame must be numeric (lon/lat).")
    }
  } else {
    stop("Input must be either an sf object or a data.frame with lon/lat in the first two columns.")
  }

  # Check whether they fall within geographic bounds
  lon_ok <- all(coords[,1] >= -180 & coords[,1] <= 180, na.rm = TRUE)
  lat_ok <- all(coords[,2] >=  -90 & coords[,2] <=  90, na.rm = TRUE)

  if (lon_ok && lat_ok) {
    message("Guessed CRS as EPSG:4326 (WGS84).")
    return(4326)
  } else {
    message("Could not guess CRS: coordinates out of geographic bounds.")
    return(NA)
  }
}


#' Transform to a metric CRS (e.g., UTM)
#'
#' Transforms an `sf` object to a metric CRS based on its centroid (UTM) or a user-specified CRS.
#'
#' @param shp An `sf` object with a defined geographic CRS or guessed via `guess_crs()`.
#' @param metric_crs Integer EPSG code for target metric CRS. If `NULL`, selected automatically by UTM zone.
#' @return A list with: `shape` (transformed sf), `coordinates` (matrix of xy), and `crs` (EPSG code).
#' @importFrom sf st_crs st_set_crs st_coordinates st_centroid st_transform
#' @examples
#' # Create test data
#' coords <- data.frame(long = c(1,2,2,4), lat = c(4,4,2,2))
#' coords_vec <- df_to_shp(coords)
#'
#' # Transform vector to a metric projection
#' transform_to_metric(coords_vec)
#'
#' @export
transform_to_metric <- function(shp, metric_crs = NULL) {
  # Ensure CRS is defined
  current_epsg <- st_crs(shp)$epsg
  if (is.na(current_epsg)) {
    message("CRS is undefined. Attempting to guess...")
    epsg <- guess_crs(shp)
    if (is.na(epsg)) stop("CRS couldn't be guessed. Use st_set_crs().")
    shp <- st_set_crs(shp, epsg)
    current_epsg <- epsg
  }

  # If already in WGS84 UTM (metric), return as‐is
  if ((current_epsg >= 32601 && current_epsg <= 32660) ||
      (current_epsg >= 32701 && current_epsg <= 32760)) {
    message(sprintf("Input is already in WGS84 UTM zone (EPSG:%d). Returning original.", current_epsg))
    coords <- st_coordinates(shp)
    return(list(shape = shp, coordinates = coords, crs = current_epsg))
  }

  # Otherwise, determine metric CRS if not provided
  if (is.null(metric_crs)) {
    # get centroid in lon/lat
    centroid <- st_coordinates(st_centroid(st_transform(shp, 4326)))
    lon <- centroid[1]; lat <- centroid[2]
    utm_zone <- floor((lon + 180) / 6) + 1
    metric_crs <- if (lat >= 0) 32600 + utm_zone else 32700 + utm_zone
    message(sprintf("Auto-selected UTM zone %d → EPSG:%d.", utm_zone, metric_crs))
  }

  # Transform and return
  shp_m <- st_transform(shp, crs = metric_crs)
  coords <- st_coordinates(shp_m)
  return(list(shape = shp_m, coordinates = coords, crs = metric_crs))
}


#' Convert spatial object to decimal degrees (EPSG:4326)
#'
#' Projects a `SpatRaster` or `sf` object into geographic coordinates.
#'
#' @param obj A `SpatRaster` or `sf` object with defined CRS.
#' @return Transformed object in EPSG:4326.
#' @importFrom terra crs project
#' @importFrom sf st_crs st_transform
#' @examples
#' # Create test data
#' df <- data.frame(long = c(277000,389000,389000,611000),
#'                  lat = c(442000,442000,221000,221000))
#' sf_obj <- st_as_sf(df, coords = names(df), crs = 32631, remove = FALSE)
#'
#' # Convert to decimal degrees
#' convert_to_decimal_degrees(sf_obj)
#' @export
convert_to_decimal_degrees <- function(obj) {
  target_crs <- "EPSG:4326"

  if (inherits(obj, "SpatRaster")) {
    message("Current raster CRS: ", terra::crs(obj))
    if (is.na(terra::crs(obj))) stop("Input raster has no CRS defined.")
    obj_out <- terra::project(obj, target_crs)
    message("Transformed raster CRS: ", terra::crs(obj_out))
    return(obj_out)
  }

  if (inherits(obj, "sf")) {
    message("Current sf CRS: ", sf::st_crs(obj)$input)
    if (is.na(sf::st_crs(obj))) stop("Input sf object has no CRS defined.")
    obj_out <- sf::st_transform(obj, crs = target_crs)
    message("Transformed sf CRS: ", sf::st_crs(obj_out)$input)
    return(obj_out)
  }

  stop("Input must be a 'SpatRaster' or 'sf' object.")
}




