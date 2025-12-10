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
transform_to_metric <- function(
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
  is_sf <- inherits(input_map, "sf")
  is_rst <- inherits(input_map, "SpatRaster")

  if (!is_sf && !is_rst) {
    stop("input_map must be either an 'sf' object or a 'SpatRaster'.")
  }

  # helper: a cool custom operator that implements if statements
  `%||%` <- function(x, y) {
    if (is.null(x) || (is.atomic(x) && length(x) == 1L && is.na(x))) y else x
  }

  # Get CRS as sf::crs object
  if (is_sf) {
    crs <- sf::st_crs(input_map)
  } else {
    crs <- sf::st_crs(terra::crs(input_map))  # parse WKT/PROJ string into sf::crs
  }

  # If CRS missing: try guessing for sf; fail for raster
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

  # Geographic vs projected?
  if (is_sf) {
    is_geographic <- sf::st_is_longlat(sf::st_geometry(input_map))
  } else {
    is_geographic <- terra::is.lonlat(input_map)
  }

  # If already projected and user wants to keep it
  if (!is_geographic && isTRUE(keep_if_projected)) {
    coords <- if (is_sf) sf::st_coordinates(input_map) else NULL
    epsg   <- crs$epsg %||% NA_integer_
    return(list(shape = input_map, coordinates = coords, crs = epsg))
  }

  # If already WGS84 UTM, return as-is
  if (!is.na(crs$epsg) &&
      ((crs$epsg >= 32601 && crs$epsg <= 32660) ||
       (crs$epsg >= 32701 && crs$epsg <= 32760))) {
    if (!quiet) {
      message(sprintf("Input already in WGS84 UTM (EPSG:%d); returning unchanged.", crs$epsg))
    }
    coords <- if (is_sf) sf::st_coordinates(input_map) else NULL
    return(list(shape = input_map, coordinates = coords, crs = crs$epsg))
  }

  # Choose target CRS if not provided
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
    ll <- sf::st_coordinates(center_ll)[1, ]
    lon_center <- ((ll[1] + 180) %% 360) - 180  # wrap to [-180, 180]
    lat_center <- ll[2]

    # bbox in lon/lat to check latitude span (robust with fallback)
    lat_span <- tryCatch({
      bb_sfc <- sf::st_as_sfc(bb)
      bb_ll <- sf::st_transform(bb_sfc, 4326)
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
      "auto" = {
        if (is.na(lat_max) || is.na(lat_min)) {
          if (!quiet) {
            message(
              "Could not evaluate latitude span; using center latitude for polar/UTM decision."
            )
          }
          if (is.na(lat_center)) {
            NA
          } else {
            (lat_center >= polar_threshold_north) || (lat_center <= polar_threshold_south)
          }
        } else {
          (lat_max >= polar_threshold_north) || (lat_min <= polar_threshold_south)
        }
      }
    )

    if (is.na(use_polar)) {
      use_polar <- FALSE
      if (!quiet) {
        message("Could not determine whether to use polar CRS; defaulting to UTM.")
      }
    }

    if (use_polar) {
      # If forced north/south, obey; else pick by center latitude sign
      if (polar == "north" || (polar == "auto" && lat_center >= 0)) {
        metric_crs <- polar_crs_north
        if (!quiet) message(sprintf("Auto-selected Arctic Polar Stereographic: EPSG:%d.", metric_crs))
      } else {
        metric_crs <- polar_crs_south
        if (!quiet) message(sprintf("Auto-selected Antarctic Polar Stereographic: EPSG:%d.", metric_crs))
      }
    } else {
      # Standard UTM
      utm_zone <- max(1, min(60, floor((lon_center + 180) / 6) + 1))
      metric_crs <- if (lat_center >= 0) 32600 + utm_zone else 32700 + utm_zone
      if (!quiet) message(sprintf("Auto-selected UTM zone %d: EPSG:%d.", utm_zone, metric_crs))
    }
  }

  # Transform and return
  target_crs <- sf::st_crs(metric_crs)

  if (is_sf) {
    input_map_m <- sf::st_transform(input_map, crs = target_crs)
    coords <- sf::st_coordinates(input_map_m)
    epsg_out <- sf::st_crs(input_map_m)$epsg %||% NA_integer_
  } else {
    # Decide interpolation method based on integerness of the raster
    integer_like <- all(terra::is.int(input_map))
    method <- if (integer_like) "near" else "bilinear"

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

  # Output
  if (is_sf){
    list(
      shape = input_map_m,
      coordinates = coords,
      crs = epsg_out
    )
  } else {
    input_map_m
  }
}


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


#' Align a SpatRaster or sf object to a template CRS or geometry
#'
#' Projects a `SpatRaster` or `sf` object (`src`) to match either the full geometry
#' (CRS, extent, resolution) of a `SpatRaster` template, or only the CRS
#' of an `sf` template.
#'
#' @param src A `SpatRaster` or `sf` object to be projected / aligned.
#' @param template Either a `SpatRaster` whose geometry should be matched,
#'   or an `sf` object whose CRS should be matched.
#' @param categorical Logical. If `TRUE`, use nearest-neighbor resampling;
#'   if `FALSE`, use bilinear interpolation.
#' @param pixel_size If `src` is `sf`, this sets the raster resolution (in map units).
#'   If NULL and `template` is a `SpatRaster`, resolution is inferred from template.
#'   Default is `NULL`.
#'
#' @return A `SpatRaster` aligned to the template.
#' @importFrom terra vect crs project rast rasterize ext compareGeom
#' @importFrom sf st_crs
#' @export
align_to <- function(src, template, categorical = TRUE, pixel_size = NULL) {
  # small helper
  `%||%` <- function(x, y) {
    if (is.null(x) || (is.atomic(x) && length(x) == 1L && is.na(x))) y else x
  }

  # interpolation method
  m <- if (categorical) "near" else "bilinear"

  # Convert sf/sfc source to SpatRaster (if needed)
  if (inherits(src, c("sf", "sfc"))) {
    src_vect <- terra::vect(src)
    m <- "near"  # sf: always categorical

    # Determine target CRS from template
    if (inherits(template, "SpatRaster")) {
      target_crs <- terra::crs(template, proj = TRUE)
    } else if (inherits(template, c("sf", "sfc"))) {
      crs_t <- sf::st_crs(template)
      if (is.na(crs_t)) stop("Template sf object has no CRS.")
      target_crs <- crs_t$wkt %||% as.character(crs_t)
    } else {
      stop("Template must be a SpatRaster or sf/sfc object.")
    }

    # Project source vector to target CRS
    src_vect <- terra::project(src_vect, target_crs)

    # Create raster template
    if (inherits(template, "SpatRaster")) {
      r_template <- terra::rast(template)
    } else {
      # sf template: create raster from extent + pixel_size
      if (is.null(pixel_size)) {
        stop("pixel_size must be provided when template is 'sf' and src is 'sf'.")
      }
      r_template <- terra::rast(
        ext = terra::ext(src_vect),
        resolution = pixel_size,
        crs = target_crs
      )
    }

    # Use the first attribute field by default
    fields <- names(src)
    if (length(fields) == 0) stop("No attribute fields found in 'sf' object.")
    field <- fields[1]

    # Rasterize sf object
    src <- terra::rasterize(src_vect, r_template, field = field)
  }

  # At this point, src *must* be a SpatRaster
  if (!inherits(src, "SpatRaster")) {
    stop("align_to(): 'src' must be a 'SpatRaster' after conversion.")
  }

  # Handle template
  # Case 1: template is SpatRaster
  if (inherits(template, "SpatRaster")) {
    # Only compare geometry when both are SpatRaster (they are)
    if (terra::compareGeom(src, template, stopOnError = FALSE)) {
      return(src)
    }
    return(terra::project(src, template, method = m))
  }

  # Case 2: template is sf/sfc
  if (inherits(template, c("sf", "sfc"))) {
    crs_t <- sf::st_crs(template)
    if (is.na(crs_t)) stop("Template sf object has no CRS.")
    target_crs <- crs_t$wkt %||% as.character(crs_t)

    # Check if CRS already matches
    if (terra::crs(src, proj = TRUE) == target_crs) {
      return(src)
    }

    return(terra::project(src, target_crs, method = m))
  }

  stop("Template must be a 'SpatRaster' or an 'sf'/'sfc' object.")
}


