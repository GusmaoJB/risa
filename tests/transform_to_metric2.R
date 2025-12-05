#' Transform to a metric CRS (UTM or Polar Stereographic fallback)
#'
#' Transforms an `sf` or `SpatRaster` object to a metric CRS. If `metric_crs`
#' is not supplied, chooses a WGS84 UTM zone from the layer's geographic center.
#' If the layer extends beyond UTM's latitude limits (default ≥84°N or ≤80°S),
#' it falls back to Polar Stereographic (WGS84: EPSG 3995 for Arctic, 3031 for
#' Antarctic).
#'
#' @param input_map An `sf` or `terra::SpatRaster` object with a defined CRS.
#'   If `sf` and CRS is missing, tries `guess_crs()` (if available).
#' @param metric_crs Target CRS (EPSG integer, WKT, proj4, or `sf::crs`).
#'   If `NULL`, it will be selected automatically.
#' @param keep_if_projected Logical; if `TRUE` and input is already projected
#'   (not lon/lat), return unchanged. Default `TRUE`.
#' @param polar One of `"auto"`, `"never"`, `"north"`, `"south"`. `"auto"` (default)
#'   uses Polar Stereographic when latitudes exceed thresholds. `"never"` forces UTM
#'   selection (even near poles). `"north"`/`"south"` force the respective polar CRS.
#' @param polar_threshold_north,polar_threshold_south Latitude thresholds (degrees)
#'   to trigger polar fallback. Defaults: `84` (north), `-80` (south).
#' @param polar_crs_north,polar_crs_south EPSG codes for polar projections.
#'   Defaults: `3995` (Arctic), `3031` (Antarctic).
#' @param quiet Logical; suppress informative messages. Default `TRUE`.
#'
#' @return A list with:
#'   - `shape`: transformed `sf` or `SpatRaster`
#'   - `coordinates`: matrix from `st_coordinates()` for `sf`, or `NULL` for rasters
#'   - `crs`: EPSG code (integer) if available, otherwise `NA_integer_`
#'
#' @importFrom sf st_crs st_set_crs st_coordinates st_transform st_geometry
#'   st_bbox st_sfc st_point st_is_longlat st_as_sfc
#' @importFrom terra crs ext is.lonlat project xmin xmax ymin ymax is.int
#' @export
transform_to_metric2 <- function(
    input_map,
    metric_crs = NULL,
    keep_if_projected = TRUE,
    polar = c("auto", "never", "north", "south"),
    polar_threshold_north = 84,
    polar_threshold_south = -80,
    polar_crs_north = 3995,  # WGS 84 / Arctic Polar Stereographic
    polar_crs_south = 3031,  # WGS 84 / Antarctic Polar Stereographic
    quiet = TRUE) {
  
  polar <- match.arg(polar)
  is_sf  <- inherits(input_map, "sf")
  is_rst <- inherits(input_map, "SpatRaster")
  
  if (!is_sf && !is_rst) {
    stop("`input_map` must be either an sf object or a SpatRaster.")
  }
  
  # helper: a cool custom operator that implements if statements
  `%||%` <- function(x, y) {
    if (is.null(x) || (is.atomic(x) && length(x) == 1L && is.na(x))) y else x
  }
  
  ## Get CRS as sf::crs object
  if (is_sf) {
    crs <- sf::st_crs(input_map)
  } else {
    crs <- sf::st_crs(terra::crs(input_map))  # parse WKT/PROJ string into sf::crs
  }
  
  ## --- If CRS missing: try guessing for sf; fail for raster ---------------
  if (is.na(crs)) {
    if (is_sf) {
      epsg <- guess_crs(input_map, quiet = TRUE)
      if (!is.na(epsg) && !is.null(epsg)) {
        input_map <- sf::st_set_crs(input_map, epsg)
        crs <- sf::st_crs(epsg)
        if (!quiet) message(sprintf("Set missing CRS to EPSG:%d via guess_crs().", epsg))
      } else {
        stop("CRS is undefined and could not be guessed.")
      }
    } else {
      stop("CRS is undefined for the SpatRaster.")
    }
  }
  
  ## Geographic vs projected? 
  if (is_sf) {
    is_geographic <- sf::st_is_longlat(sf::st_geometry(input_map))
  } else {
    is_geographic <- terra::is.lonlat(input_map)
  }
  
  ## If already projected and user wants to keep it
  if (!is_geographic && isTRUE(keep_if_projected)) {
    coords <- if (is_sf) sf::st_coordinates(input_map) else NULL
    epsg   <- crs$epsg %||% NA_integer_
    return(list(shape = input_map, coordinates = coords, crs = epsg))
  }
  
  ## If already WGS84 UTM, return as-is
  if (!is.na(crs$epsg) &&
      ((crs$epsg >= 32601 && crs$epsg <= 32660) ||
       (crs$epsg >= 32701 && crs$epsg <= 32760))) {
    if (!quiet) {
      message(sprintf("Input already in WGS84 UTM (EPSG:%d); returning unchanged.", crs$epsg))
    }
    coords <- if (is_sf) sf::st_coordinates(input_map) else NULL
    return(list(shape = input_map, coordinates = coords, crs = crs$epsg))
  }
  
  ## Choose target CRS if not provided 
  if (is.null(metric_crs)) {
    # bbox in native CRS
    if (is_sf) {
      bb <- sf::st_bbox(input_map)
    } else {
      bb <- sf::st_bbox(
        c(
          xmin = terra::xmin(input_map),
          ymin = terra::ymin(input_map),
          xmax = terra::xmax(input_map),
          ymax = terra::ymax(input_map)
        ),
        crs = crs
      )
    }
    
    cx <- (bb["xmin"] + bb["xmax"]) / 2
    cy <- (bb["ymin"] + bb["ymax"]) / 2
    
    # center -> lon/lat for UTM decision
    center_pt <- sf::st_sfc(sf::st_point(c(cx, cy)), crs = crs)
    center_ll <- sf::st_transform(center_pt, 4326)
    ll        <- sf::st_coordinates(center_ll)[1, ]
    lon_center <- ((ll[1] + 180) %% 360) - 180  # wrap to [-180, 180]
    lat_center <- ll[2]
    
    # bbox in lon/lat to check latitude span (robust with fallback)
    lat_span <- tryCatch({
      bb_sfc   <- sf::st_as_sfc(bb)
      bb_ll    <- sf::st_transform(bb_sfc, 4326)
      bb_ll_bb <- sf::st_bbox(bb_ll)
      c(lat_min = bb_ll_bb["ymin"], lat_max = bb_ll_bb["ymax"])
    }, error = function(e) {
      if (!quiet) {
        message("Warning: could not transform bbox to lon/lat; using center latitude only for polar/UTM decision.")
      }
      c(lat_min = lat_center, lat_max = lat_center)
    })
    
    lat_min <- lat_span["lat_min"]
    lat_max <- lat_span["lat_max"]
    
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
      utm_zone  <- max(1, min(60, floor((lon_center + 180) / 6) + 1))
      metric_crs <- if (lat_center >= 0) 32600 + utm_zone else 32700 + utm_zone
      if (!quiet) message(sprintf("Auto-selected UTM zone %d → EPSG:%d.", utm_zone, metric_crs))
    }
  }
  
  ## Transform and return
  target_crs <- sf::st_crs(metric_crs)
  
  if (is_sf) {
    input_map_m <- sf::st_transform(input_map, crs = target_crs)
    coords   <- sf::st_coordinates(input_map_m)
    epsg_out <- sf::st_crs(input_map_m)$epsg %||% NA_integer_
  } else {
    # Decide interpolation method based on integerness of the raster
    integer_like <- all(terra::is.int(input_map))
    method       <- if (integer_like) "near" else "bilinear"
    
    if (!quiet) {
      msg <- if (integer_like) "Integer-like raster; using 'near' interpolation." 
      else "Non-integer raster; using 'bilinear' interpolation."
      message(msg)
    }
    
    terra_crs <- target_crs$wkt %||% target_crs$proj4string %||% as.character(metric_crs)
    input_map_m <- terra::project(input_map, terra_crs, method = method)
    coords <- NULL
    epsg_out <- sf::st_crs(terra::crs(input_map_m))$epsg %||% NA_integer_
  }
  
  list(
    shape       = input_map_m,
    coordinates = coords,
    crs         = epsg_out
  )
}

