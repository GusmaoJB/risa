#' Guess coordinate reference system from geographic bounds
#'
#' If coordinates look like degrees (lon in [-180, 180], lat in [-90, 90]),
#' returns EPSG:4326; else `NA`. If `shp` is `sf` with a geographic CRS set,
#' returns its EPSG when available.
#'
#' @param shp An `sf` object or a data frame with longitude/latitude.
#' @param lon,lat Optional column names or indices for data frames; if omitted,
#'   the first two numeric columns are used.
#' @param quiet Logical; suppress messages (default `TRUE`).
#' @return Integer EPSG code (e.g., `4326L`) or `NA_integer_` if not guessed.
#' @importFrom sf st_coordinates st_crs st_is_longlat st_geometry
#' @examples
#' # Create test data
#' coords <- data.frame(long = c(1,2,2,4), lat = c(4,4,2,2))
#' # Gussing coordinates
#' guess_crs(coords)
#' @export
guess_crs <- function(shp, lon = NULL, lat = NULL, quiet = TRUE) {
  # Fast path: sf with a known geographic CRS
  if (inherits(shp, "sf")) {
    crs <- sf::st_crs(shp)
    if (!is.na(crs) && isTRUE(sf::st_is_longlat(sf::st_geometry(shp)))) {
      if (!is.na(crs$epsg)) return(as.integer(crs$epsg))
      # fall through to bounds heuristic if EPSG unknown
    }
    coords <- sf::st_coordinates(shp)[, 1:2, drop = FALSE]
  } else if (is.data.frame(shp)) {
    nms <- names(shp)
    if (is.null(lon) || is.null(lat)) {
      num_idx <- which(vapply(shp, is.numeric, logical(1)))
      if (length(num_idx) < 2L) stop("Data frame needs two numeric columns for lon/lat (or specify `lon` and `lat`).")
      lon <- nms[num_idx[1]]; lat <- nms[num_idx[2]]
    } else {
      if (is.numeric(lon)) lon <- nms[lon]
      if (is.numeric(lat)) lat <- nms[lat]
    }
    if (!all(c(lon, lat) %in% nms)) stop("`lon`/`lat` columns not found in data frame.")
    if (!is.numeric(shp[[lon]]) || !is.numeric(shp[[lat]])) stop("`lon` and `lat` must be numeric.")
    coords <- cbind(shp[[lon]], shp[[lat]])
  } else {
    stop("Input must be an `sf` object or a data frame with lon/lat.")
  }

  # Require at least one finite pair
  keep <- is.finite(coords[, 1]) & is.finite(coords[, 2])
  if (!any(keep)) {
    if (!quiet) message("Could not guess CRS: no finite coordinates found.")
    return(NA_integer_)
  }

  xr <- range(coords[keep, 1], na.rm = TRUE)
  yr <- range(coords[keep, 2], na.rm = TRUE)
  lon_ok <- xr[1] >= -180 && xr[2] <= 180
  lat_ok <- yr[1] >=  -90 && yr[2] <=  90

  if (lon_ok && lat_ok) {
    if (!quiet) message("Guessed CRS as EPSG:4326 (WGS84).")
    return(4326L)
  } else {
    if (!quiet) message("Could not guess CRS as geographic; returning NA.")
    return(NA_integer_)
  }
}


#' Transform to a metric CRS (UTM or Polar Stereographic fallback)
#'
#' Transforms an `sf` layer to a metric CRS. If `metric_crs` is not supplied,
#' chooses a WGS84 UTM zone from the layer's geographic center. If the layer
#' extends beyond UTM's latitude limits (default ≥84°N or ≤80°S), it falls back
#' to Polar Stereographic (WGS84: EPSG 3995 for Arctic, 3031 for Antarctic).
#'
#' @param shp An `sf` object with a defined CRS. If missing, tries `guess_crs()` (if available).
#' @param metric_crs Target CRS (EPSG integer, WKT, or `sf::crs`). If `NULL`, it will be selected automatically.
#' @param keep_if_projected Logical; if `TRUE` and input is already projected (not lon/lat),
#'   return unchanged. Default `FALSE`.
#' @param polar One of `"auto"`, `"never"`, `"north"`, `"south"`. `"auto"` (default)
#'   uses Polar Stereographic when latitudes exceed thresholds. `"never"` forces UTM
#'   selection (even near poles). `"north"`/`"south"` force the respective polar CRS.
#' @param polar_threshold_north,polar_threshold_south Latitude thresholds (degrees) to trigger
#'   polar fallback. Defaults: `84` (north), `-80` (south).
#' @param polar_crs_north,polar_crs_south EPSG codes for polar projections.
#'   Defaults: `3995` (Arctic), `3031` (Antarctic).
#' @param quiet Logical; suppress informative messages. Default `TRUE`.
#' @return A list with: `shape` (transformed `sf`), `coordinates` (matrix from `st_coordinates()`), and `crs` (EPSG or `NA`).
#' @importFrom sf st_crs st_set_crs st_coordinates st_transform st_geometry st_bbox st_sfc st_point st_is_longlat st_as_sfc
#' @examples
#' # Create test data
#' coords <- data.frame(long = c(1,2,2,4), lat = c(4,4,2,2))
#' coords_vec <- df_to_shp(coords)
#'
#' # Transform vector to a metric projection
#' transform_to_metric(coords_vec)
#' @export
transform_to_metric <- function(
    shp,
    metric_crs = NULL,
    keep_if_projected = FALSE,
    polar = c("auto", "never", "north", "south"),
    polar_threshold_north = 84,
    polar_threshold_south = -80,
    polar_crs_north = 3995,  # WGS 84 / Arctic Polar Stereographic
    polar_crs_south = 3031,  # WGS 84 / Antarctic Polar Stereographic
    quiet = TRUE
) {
  polar <- match.arg(polar)
  if (!inherits(shp, "sf")) stop("`shp` must be an sf object.")

  crs <- sf::st_crs(shp)

  # If CRS missing, attempt to guess only for long/lat-looking coords
  if (is.na(crs)) {
    epsg <- tryCatch({
      if (exists("guess_crs", mode = "function")) guess_crs(shp, quiet = TRUE) else NA_integer_
    }, error = function(e) NA_integer_)
    if (!is.na(epsg)) {
      shp <- sf::st_set_crs(shp, epsg)
      crs <- sf::st_crs(epsg)
      if (!quiet) message(sprintf("Set missing CRS to EPSG:%d via guess.", epsg))
    } else {
      stop("CRS is undefined and could not be guessed. Set it explicitly with `sf::st_set_crs()`.")
    }
  }

  # If already projected and user wants to keep it, return early
  is_geographic <- sf::st_is_longlat(sf::st_geometry(shp))
  if (!is_geographic && isTRUE(keep_if_projected)) {
    coords <- sf::st_coordinates(shp)
    return(list(shape = shp, coordinates = coords, crs = crs$epsg %||% NA_integer_))
  }

  # If already WGS84 UTM, return as-is
  if (!is.na(crs$epsg) && ((crs$epsg >= 32601 && crs$epsg <= 32660) ||
                           (crs$epsg >= 32701 && crs$epsg <= 32760))) {
    if (!quiet) message(sprintf("Input already in WGS84 UTM (EPSG:%d); returning unchanged.", crs$epsg))
    coords <- sf::st_coordinates(shp)
    return(list(shape = shp, coordinates = coords, crs = crs$epsg))
  }

  # Choose target CRS if not provided
  if (is.null(metric_crs)) {

    # Compute bbox and its center in current CRS
    bb <- sf::st_bbox(shp)
    cx <- (bb["xmin"] + bb["xmax"]) / 2
    cy <- (bb["ymin"] + bb["ymax"]) / 2

    # Center point to lon/lat for UTM decision
    center_pt <- sf::st_sfc(sf::st_point(c(cx, cy)), crs = crs)
    center_ll <- sf::st_transform(center_pt, 4326)
    ll <- sf::st_coordinates(center_ll)[1, ]
    lon_center <- ((ll[1] + 180) %% 360) - 180  # wrap to [-180,180]
    lat_center <- ll[2]

    # Bbox in lon/lat to check latitude span (cheap and robust)
    bb_sfc   <- sf::st_as_sfc(bb, crs = crs)
    bb_ll    <- sf::st_transform(bb_sfc, 4326)
    bb_ll_bb <- sf::st_bbox(bb_ll)
    lat_min  <- bb_ll_bb["ymin"]
    lat_max  <- bb_ll_bb["ymax"]

    # Decide on polar fallback
    use_polar <- switch(
      polar,
      "never" = FALSE,
      "north" = TRUE,
      "south" = TRUE,
      "auto"  = (lat_max >= polar_threshold_north) || (lat_min <= polar_threshold_south)
    )

    if (use_polar) {
      # If forced north/south, obey; else pick by center latitude sign
      if (polar == "north" || (polar == "auto" && lat_center >= 0)) {
        metric_crs <- polar_crs_north
        if (!quiet) message(sprintf("Auto-selected Arctic Polar Stereographic → EPSG:%d.", metric_crs))
      } else {
        metric_crs <- polar_crs_south
        if (!quiet) message(sprintf("Auto-selected Antarctic Polar Stereographic → EPSG:%d.", metric_crs))
      }
    } else {
      # Standard UTM
      utm_zone <- max(1, min(60, floor((lon_center + 180) / 6) + 1))
      metric_crs <- if (lat_center >= 0) 32600 + utm_zone else 32700 + utm_zone
      if (!quiet) message(sprintf("Auto-selected UTM zone %d → EPSG:%d.", utm_zone, metric_crs))
    }
  }

  # Transform and return
  target <- sf::st_crs(metric_crs)
  shp_m  <- sf::st_transform(shp, crs = target)
  coords <- sf::st_coordinates(shp_m)
  list(shape = shp_m, coordinates = coords, crs = sf::st_crs(shp_m)$epsg %||% NA_integer_)
}

# helper
`%||%` <- function(x, y) if (is.null(x) || (is.atomic(x) && length(x) == 1L && is.na(x))) y else x


#' Convert spatial object to decimal degrees (EPSG:4326)
#'
#' Reprojects a `SpatRaster` or `sf` object to geographic coordinates (EPSG:4326).
#' Note: for `sf`, only the geometry is transformed; any original numeric
#' columns you kept (e.g., Easting/Northing) are not altered.
#'
#' @param obj A `terra::SpatRaster` or `sf` object with a defined CRS.
#' @param method Resampling method used for rasters. One of
#'   `"near"`, `"bilinear"`, `"cubic"`, `"cubicspline"`, `"lanczos"`. Ignored for `sf`.
#' @param quiet Logical; suppress messages. Default `TRUE`.
#' @return Object reprojected to EPSG:4326 (same class as input).
#' @importFrom terra crs project same.crs
#' @importFrom sf st_crs st_transform
#' @examples
#' df <- data.frame(long = c(277000,389000,389000,611000),
#'                  lat  = c(442000,442000,221000,221000))
#' sf_obj <- sf::st_as_sf(df, coords = names(df), crs = 32631, remove = FALSE)
#' convert_to_decimal_degrees(sf_obj)
#' @export
convert_to_decimal_degrees <- function(obj, method = "near", quiet = TRUE) {
  target_crs <- "EPSG:4326"

  # If SpatRaster
  if (inherits(obj, "SpatRaster")) {
    crs_in <- terra::crs(obj)
    if (is.na(crs_in) || crs_in == "") stop("Input raster has no CRS defined.")

    # No-op if already 4326
    if (isTRUE(try(terra::same.crs(obj, target_crs), silent = TRUE))) {
      if (!quiet) message("Raster already in EPSG:4326; returning unchanged.")
      return(obj)
    }

    if (!quiet) message(sprintf("Reprojecting raster → %s using '%s' resampling.", target_crs, method))
    return(terra::project(obj, target_crs, method = method))
  }

  # If sf
  if (inherits(obj, "sf")) {
    crs_in <- sf::st_crs(obj)
    if (is.na(crs_in)) stop("Input sf object has no CRS defined.")

    # No-op if already 4326
    if (!is.na(crs_in$epsg) && crs_in$epsg == 4326L) {
      if (!quiet) message("sf already in EPSG:4326; returning unchanged.")
      return(obj)
    }

    if (!quiet) message("Reprojecting sf → EPSG:4326.")
    return(sf::st_transform(obj, crs = 4326))
  }

  stop("Input must be a 'SpatRaster' or 'sf' object.")
}
