################################################################################
###############################input_tools######################################
################################################################################

#' Merge a list of 'sf' objects into a data frame of coordinates
#'
#' Combines multiple vector objects of class `sf` into a single data frame with X, Y, and group labels.
#' Coordinates are the vertices of each geometry (POINT/LINESTRING/POLYGON/MULTI*).
#' For polygons, ring-closing vertices are included as returned by `sf::st_coordinates()`.
#'
#' @param shp_list A list of `sf` objects.
#' @param group_size Optional name of a column in each `sf` to repeat per vertex and include in the output.
#' @returns A data.frame with columns `X`, `Y`, `group`, and optionally `group_size`.
#' @export
merge_shp <- function(shp_list, group_size = NULL) {
  if (!is.list(shp_list) || length(shp_list) == 0L) {
    stop("`shp_list` must be a non-empty list of sf objects.")
  }
  if (!all(vapply(shp_list, inherits, logical(1), what = "sf"))) {
    stop("All elements of `shp_list` must be sf objects.")
  }

  # Warn if CRSs differ; coordinates will be mixed as-is
  crs_keys <- vapply(shp_list, function(s) {
    crs <- sf::st_crs(s)
    if (!is.null(crs$epsg)) paste0("epsg:", crs$epsg) else crs$wkt %||% "NA"
  }, character(1))
  if (length(unique(crs_keys)) > 1L) {
    warning("Input layers have different CRS; coordinates are kept as-is. Consider transforming beforehand.")
  }

  out_list <- vector("list", length(shp_list))

  for (i in seq_along(shp_list)) {
    shp <- shp_list[[i]]

    group_name <-
      if (!is.null(names(shp_list)) && nzchar(names(shp_list)[i])) names(shp_list)[i] else as.character(i)

    # One row per vertex
    coords <- sf::st_coordinates(shp)

    # Base columns
    data_out <- data.frame(
      X = coords[, 1],
      Y = coords[, 2],
      group = rep(group_name, nrow(coords)),
      stringsAsFactors = FALSE
    )

    # Repeat attribute per vertex if requested and present
    if (!is.null(group_size) && group_size %in% names(shp)) {
      # Number of vertices per feature
      n_per_feat <- sf::st_npoints(sf::st_geometry(shp), by_feature = TRUE)
      data_out[[group_size]] <- rep(shp[[group_size]], times = n_per_feat)
    }

    out_list[[i]] <- data_out
  }

  out <- do.call(rbind, out_list)
  row.names(out) <- NULL
  out
}

# A helper function to deal with null values
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Convert a data.frame to an sf object
#'
#' Builds an `sf` POINT layer from longitude/latitude columns of a data frame.
#' If `crs` is not supplied, tries `guess_crs()`. If guessing fails,
#' falls back to EPSG:4326 only when coordinates look like degrees; otherwise uses `NA` CRS.
#'
#' @param df A data.frame with at least two numeric columns for coordinates, or an `sf` (returned as-is).
#' @param lon,lat Optional column names or indices for longitude and latitude.
#'   If omitted, will try common names (lon/long/longitude/x, lat/latitude/y), else the first two numeric columns.
#' @param crs Optional CRS (e.g., EPSG code like 4326, WKT, or an `sf::crs` object).
#' @param guess Logical; if `TRUE`, try `guess_crs()` when `crs` is not provided. Default `TRUE`.
#' @param keep_coords Logical; keep original lon/lat columns (`TRUE`, default) or drop them (`FALSE`).
#' @param drop_na Logical; drop rows with non-finite lon/lat (`TRUE`, default). If `FALSE`, will error on NA/Inf.
#' @param quiet Logical; suppress messages/warnings (default `TRUE`). Set to `FALSE` for informative messages.
#' @return An `sf` object with geometry column and original attributes.
#' @importFrom sf st_as_sf st_crs
#' @examples
#' df <- data.frame(long = c(1,2,2,4), lat = c(4,4,2,2))
#' vec <- df_to_shp(df)             # silent
#' vec2 <- df_to_shp(df, quiet=FALSE) # verbose
#' @export
df_to_shp <- function(df,
                      lon = NULL, lat = NULL,
                      crs = NULL,
                      guess = TRUE,
                      keep_coords = TRUE,
                      drop_na = TRUE,
                      quiet = TRUE) {

  # 1) Early return for sf ------------------------------------------------------
  if (inherits(df, "sf")) {
    if (!quiet) message("Input is already an 'sf' object; returning as-is.")
    return(df)
  }
  if (!is.data.frame(df)) stop("`df` must be a data.frame.")

  # 2) Choose lon/lat columns ---------------------------------------------------
  nms <- names(df)
  if (is.null(lon) || is.null(lat)) {
    canon  <- tolower(nms)
    lon_ix <- which(canon %in% c("lon","long","longitude","x","longitud"))
    lat_ix <- which(canon %in% c("lat","latitude","y","latitud"))
    if (length(lon_ix) >= 1 && length(lat_ix) >= 1) {
      lon <- nms[lon_ix[1]]; lat <- nms[lat_ix[1]]
    } else {
      num_idx <- which(vapply(df, is.numeric, logical(1)))
      if (length(num_idx) < 2) stop("Could not find lon/lat: supply `lon`/`lat` or add two numeric columns.")
      lon <- nms[num_idx[1]]; lat <- nms[num_idx[2]]
      if (!quiet) message(sprintf("Using numeric columns '%s' (lon) and '%s' (lat); supply `lon`/`lat` to control.", lon, lat))
    }
  } else {
    if (is.numeric(lon)) lon <- nms[lon]
    if (is.numeric(lat)) lat <- nms[lat]
  }
  if (!all(c(lon, lat) %in% nms)) stop("`lon`/`lat` columns not found in `df`.")
  if (!is.numeric(df[[lon]]) || !is.numeric(df[[lat]])) stop("`lon` and `lat` must be numeric.")

  # 3) Handle NA/Inf ------------------------------------------------------------
  if (drop_na) {
    keep <- is.finite(df[[lon]]) & is.finite(df[[lat]])
    if (!all(keep)) df <- df[keep, , drop = FALSE]
  } else if (any(!is.finite(df[[lon]]) | !is.finite(df[[lat]]))) {
    stop("`lon`/`lat` contain NA/Inf; set `drop_na = TRUE` to drop them.")
  }

  # 4) Determine CRS ------------------------------------------------------------
  target_crs <- NULL
  if (!is.null(crs)) {
    target_crs <- sf::st_crs(crs)
  } else if (isTRUE(guess) && exists("guess_crs", mode = "function")) {
    guessed <- tryCatch(
      if (quiet) suppressMessages(guess_crs(df)) else guess_crs(df),
      error = function(e) NA
    )
    target_crs <- tryCatch(sf::st_crs(guessed), error = function(e) NA)
  }
  if (is.null(target_crs) || is.na(target_crs)) {
    rng_lon <- range(df[[lon]], na.rm = TRUE)
    rng_lat <- range(df[[lat]], na.rm = TRUE)
    if (all(rng_lon >= -180 & rng_lon <= 180) && all(rng_lat >= -90 & rng_lat <= 90)) {
      target_crs <- sf::st_crs(4326)
    } else {
      if (!quiet) message("CRS could not be determined; creating geometry with NA CRS.")
      target_crs <- NA
    }
  }

  # 5) Build sf -----------------------------------------------------------------
  sf::st_as_sf(
    df,
    coords = c(lon, lat),
    crs = target_crs,
    remove = !keep_coords
  )
}


#' Split a data.frame or sf into a list of sf by group
#'
#' Splits rows into separate `sf` objects based on a grouping column.
#' If `df` is a data.frame, it's converted once with `df_to_shp()` and then split.
#'
#' @param df A data.frame (with lon/lat columns) or an `sf` object.
#' @param group Optional grouping column (name or index). If omitted:
#'   - for data.frame with ≥3 columns, uses the 3rd column;
#'   - for sf, returns a single-element list (no split).
#' @param drop_na_group Logical; if `FALSE`, rows with NA in `group` are returned
#'   as a `<NA>` element. Default `TRUE`.
#' @param ... Passed to `df_to_shp()` when `df` is a data.frame (e.g., `lon`, `lat`, `crs`).
#' @return A named list of `sf` objects, one per group level (or a single item if no grouping).
#' @examples
#' df <- data.frame(
#'   long = c(1,2,2,4, 2,5,4,6),
#'   lat  = c(4,4,2,2, 4,4,2,2),
#'   species = rep(c("sp1","sp2"), each = 4)
#' )
#' vec_list <- df_to_list(df, group = "species")
#' @export
df_to_list <- function(df, group = NULL, drop_na_group = TRUE, ...) {
  nm <- deparse(substitute(df))

  # Normalize to sf once --------------------------------------------------------
  if (inherits(df, "sf")) {
    sf_obj <- df
  } else if (is.data.frame(df)) {
    # Default grouping for data.frame: 3rd column if not provided
    if (is.null(group) && ncol(df) >= 3L) group <- names(df)[3L]
    sf_obj <- df_to_shp(df, ...)  # uses your improved df_to_shp
  } else {
    stop("`df` must be a data.frame or an sf object.")
  }

  # If no grouping requested/available, return single item ---------------------
  if (is.null(group)) {
    lst <- list(sf_obj)
    names(lst) <- nm
    return(lst)
  }

  # Resolve group column name safely -------------------------------------------
  grp_name <- if (is.numeric(group)) names(sf_obj)[group] else group
  sf_col   <- attr(sf_obj, "sf_column")
  if (!grp_name %in% names(sf_obj) || grp_name == sf_col) {
    stop("`group` must refer to a non-geometry column present in `df`/`sf`.")
  }

  g <- sf_obj[[grp_name]]

  # Split efficiently; control NA handling -------------------------------------
  out <- split(sf_obj, g, drop = TRUE)
  if (!drop_na_group && anyNA(g)) {
    out[["<NA>"]] <- sf_obj[is.na(g), , drop = FALSE]
  }

  # If nothing to split (single level), still return a named list
  if (length(out) == 0L) {
    out <- list(sf_obj)
    names(out) <- nm
  }

  out
}


#' Compute the nesting depth of a list
#'
#' Determines the maximum number of nested list layers in an R object.
#'
#' @param x An R object. If it is a list (possibly nested), its maximum nesting depth is returned; otherwise 0.
#' @return An integer scalar: the depth of list nesting (0 for non‐lists, 1 for an empty list, etc.).
#' @examples
#' # Not a list
#' list_depth_base(42)
#'
#' # Empty list
#' list_depth_base(list())
#'
#' # Nested list
#' nested <- list(a = list(b = list(c = 1)))
#' list_depth_base(nested)
#'
#' @export
list_depth_base <- function(x) {
  if (!is.list(x) || length(x) == 0) {
    return(if (is.list(x)) 1L else 0L)
  }
  1L + max(vapply(x, list_depth_base, integer(1)))
}


################################################################################
##########################projection_tools######################################
################################################################################

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
#' @export
guess_crs <- function(shp, lon = NULL, lat = NULL, quiet = TRUE) {
  # 1) Fast path: sf with a known geographic CRS
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

  # 2) Require at least one finite pair
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

  # --- SpatRaster --------------------------------------------------------------
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

  # --- sf ----------------------------------------------------------------------
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


################################################################################
###############################create_area######################################
################################################################################

#' Create polygon of the area of interest
#'
#' Generate an area polygon (convex hull or bounding box) around input points, with optional scaling.
#' The scaling is an affine expansion from the centroid by a factor of (1 + buffer_frac).
#'
#' @param x An `sf` object or data frame with longitude/latitude in the first two columns.
#' @param crs Integer or string; EPSG or proj string for the metric transform. If `NULL`, UTM is chosen automatically.
#' @param area_type One of `"convex_hull"` or `"bbox"`. Default `"convex_hull"`.
#' @param buffer_frac Numeric ≥ 0; fractional expansion relative to the centroid (e.g., 0.5 = +50%).
#' @param quiet Logical; suppress informative messages. Default `TRUE`.
#' @return An `sf` POLYGON.
#' @importFrom sf st_crs st_bbox st_sfc st_polygon st_as_sf
#' @export
create_area <- function(x,
                        crs = NULL,
                        area_type = c("convex_hull", "bbox"),
                        buffer_frac = 0.5,
                        quiet = TRUE) {
  area_type <- match.arg(area_type)
  if (!is.numeric(buffer_frac) || length(buffer_frac) != 1L || buffer_frac < 0)
    stop("`buffer_frac` must be a single non-negative number.")

  # — accept sf or data.frame --------------------------------------------------
  if (inherits(x, "sf")) {
    x_sf <- x
  } else if (inherits(x, "data.frame")) {
    x_sf <- df_to_shp(x)
  } else {
    stop("`x` must be an sf or a data.frame.")
  }

  # — project to metric --------------------------------------------------------
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


################################################################################
###############################reclass_matrix###################################
################################################################################

#' Generate reclassification matrix for raster values
#'
#' Builds a matrix mapping value ranges to discrete class codes,
#' optionally excluding lowest 5 percent.
#'
#' @param raster A `terra::SpatRaster` object.
#' @param n_classes Integer 1–10 number of classes.
#' @param exclude_lowest Logical; if `TRUE`, values below 5 percent of max are mapped to `NA`.
#' @param lowest_prop Numeric; Defines the cutoff for the lowest proportion  (default `0.05`).
#' @returns A numeric matrix with columns `from`, `to`, and `class` for use with `terra::classify()`.
#' @importFrom terra global
#' @examples
#' # Example
#' r <- terra::rast(nrows=5, ncols=5, vals=c(0.01,2,3,4,5))
#' reclass_matrix(r)
#' @export
reclass_matrix <- function(raster, n_classes = 3, exclude_lowest = TRUE, lowest_prop = 0.05) {
  # Input checks
  if (!inherits(raster, "SpatRaster")) stop("`raster` must be a terra::SpatRaster.")
  if (!is.numeric(n_classes) || n_classes %% 1 != 0 || n_classes < 1 || n_classes > 10) {
    stop("`n_classes` must be an integer between 1 and 10.")
  }
  if (!is.logical(exclude_lowest) || length(exclude_lowest) != 1L) {
    stop("`exclude_lowest` must be a single logical (TRUE/FALSE).")
  }

  # Min/Max without loading all values
  mm <- terra::global(raster, fun = c("min", "max"), na.rm = TRUE)
  min_val <- as.numeric(mm[1, "min"])
  max_val <- as.numeric(mm[1, "max"])
  if (!is.finite(min_val) || !is.finite(max_val)) {
    stop("Raster contains no finite values (all NA/Inf).")
  }

  # Cutoff for the lowest proportion
  val_05 <- lowest_prop * max_val

  # Lower break for the classed bins
  lower_break <- if (exclude_lowest) val_05 else min_val
  if (lower_break > max_val) lower_break <- min_val  # guard weird cases

  # Equal-width breaks -> n_classes bins
  breaks <- seq(lower_break, max_val, length.out = n_classes + 1)

  # Build bins: make from/to the same length (n_classes)
  from <- breaks[-length(breaks)]
  to   <- breaks[-1]
  to[length(to)] <- Inf  # last bin open-ended
  cls  <- seq_len(n_classes)

  mat <- cbind(from = from, to = to, class = cls)

  # Prepend NA bin for excluded-lowest portion
  if (exclude_lowest) {
    mat <- rbind(c(from = -Inf, to = val_05, class = NA_real_), mat)
  }

  storage.mode(mat) <- "double"
  mat
}


################################################################################
############################get_class_kernel####################################
################################################################################

#' Kernel density estimation and reclassification wrapper
#'
#' Computes a KDE for point data in a metric CRS, reclassifies values into
#' `n_classes`, and returns a raster (`SpatRaster`), a vector (`sf` polygons),
#' or both.
#'
#' @param x An `sf` point layer or a data.frame/tibble with lon/lat in two columns.
#' @param area Bounding region (`sf`, `bbox`, or data.frame) or `NULL` (auto bbox with buffer).
#' @param n_classes Integer number of reclassification bins. Default 3.
#' @param output_layer_type One of `"shp"`, `"raster"`, or `"both"`. Default `"shp"`.
#' @param radius KDE bandwidth (same units as projected coordinates). If `NULL`,
#'   uses `radius_method`. Default `NULL`.
#' @param radius_method One of `"nndist"` (median NN × 1.5), `"ppl"` (profile
#'   likelihood via `spatstat.explore::bw.ppl`), or `"fixed"`. Default `"nndist"`.
#' @param group_size Optional column name in `x` to use as non-negative weights.
#' @param pixel_size Optional target pixel size (meters). If `NULL`, uses `dimyx`.
#' @param dimyx Optional `(ny, nx)` grid size for the KDE (default `c(512,512)`).
#'   Ignored if `pixel_size` is provided.
#' @param exclude_lowest Logical; passed to `reclass_matrix()` (default `TRUE`).
#' @param lowest_prop Numeric; passed to `reclass_matrix()` (default `0.05`).
#' @param return_crs One of `"metric"` (default; same CRS as used for KDE) or `"4326"`
#'   to reproject outputs to WGS84. Rasters use nearest-neighbor to preserve classes.
#' @param quiet Suppress messages. Default `TRUE`.
#' @return `sf`, `SpatRaster`, or `list(raster=..., shp=...)` depending on `output_layer_type`.
#' @importFrom sf st_as_sf st_as_sfc st_crs st_transform
#' @importFrom spatstat.geom as.owin ppp nndist
#' @importFrom spatstat.explore density.ppp bw.ppl
#' @importFrom terra rast crs classify as.polygons mask vect project
#' @export
get_class_kernel <- function(
    x,
    area = NULL,
    n_classes = 3,
    output_layer_type = c("shp", "raster", "both"),
    radius = NULL,
    radius_method = c("nndist", "ppl", "fixed"),
    group_size = NULL,
    pixel_size = NULL,
    dimyx = c(512, 512),
    exclude_lowest = TRUE,
    lowest_prop = 0.05,
    return_crs = c("metric","4326"),
    quiet = TRUE
) {
  output_layer_type <- match.arg(output_layer_type)
  radius_method     <- match.arg(radius_method)
  return_crs        <- match.arg(return_crs)

  # --- normalize input to sf ---------------------------------------------------
  if (inherits(x, "sf")) {
    x_sf <- x
  } else if (inherits(x, c("data.frame", "tbl_df", "tbl"))) {
    if (ncol(x) < 2 || !all(vapply(x[,1:2], is.numeric, logical(1)))) {
      stop("For data.frame input, the first two columns must be numeric lon/lat.")
    }
    if (!quiet) message("Converting input data.frame to sf...")
    x_sf <- df_to_shp(x)
  } else {
    stop("`x` must be an sf object or a data.frame with lon/lat.")
  }

  # --- area handling -----------------------------------------------------------
  if (is.null(area)) {
    if (!quiet) message("No area provided. Using buffered bounding box around observations...")
    area <- create_area(x_sf, area_type = "bbox", buffer_frac = 0.5, quiet = quiet)
  } else if (inherits(area, "bbox")) {
    area <- sf::st_as_sfc(area)
  } else if (inherits(area, "data.frame")) {
    area <- df_to_shp(area)
  } else if (!inherits(area, "sf")) {
    stop("`area` must be NULL, bbox, data.frame, or sf.")
  }

  # --- project both to the same metric CRS ------------------------------------
  area_m <- transform_to_metric(area, quiet = quiet)
  crs_m  <- sf::st_crs(area_m$shape)$epsg
  x_m    <- transform_to_metric(x_sf, metric_crs = crs_m, quiet = quiet)
  coords <- x_m$coordinates[, 1:2, drop = FALSE]
  if (nrow(coords) == 0L) stop("No points available after preprocessing.")
  if (nrow(unique(coords)) == 1L && is.null(radius)) {
    warning("Only one unique point; KDE will be a single peak. Consider setting `radius` explicitly.")
  }

  # --- weights -----------------------------------------------------------------
  weights <- NULL
  if (!is.null(group_size)) {
    if (!group_size %in% names(x_m$shape)) stop("`group_size` not found in `x`.")
    weights <- x_m$shape[[group_size]]
    if (!is.numeric(weights)) stop("`group_size` must be numeric.")
    if (any(weights < 0, na.rm = TRUE)) stop("`group_size` weights must be non-negative.")
  }

  # --- window & bandwidth ------------------------------------------------------
  W <- spatstat.geom::as.owin(area_m$shape)
  X <- spatstat.geom::ppp(x = coords[,1], y = coords[,2], window = W, check = FALSE)

  if (is.null(radius)) {
    if (radius_method == "nndist") {
      nn <- spatstat.geom::nndist(X)
      radius <- 1.5 * stats::median(nn[is.finite(nn)], na.rm = TRUE)
      if (!quiet) message(sprintf("Auto bandwidth (nndist): sigma = %.3f", radius))
    } else if (radius_method == "ppl") {
      bw <- spatstat.explore::bw.ppl(X)
      radius <- as.numeric(bw)
      if (!quiet) message(sprintf("Auto bandwidth (ppl): sigma = %.3f", radius))
    } else {
      stop("`radius` is NULL but `radius_method = 'fixed'`. Supply a numeric radius.")
    }
  }
  if (!is.numeric(radius) || length(radius) != 1L || !(radius > 0)) {
    stop("`radius` must be a single positive number (in projected units).")
  }

  # --- grid resolution ---------------------------------------------------------
  if (!is.null(pixel_size)) {
    if (!(is.numeric(pixel_size) && length(pixel_size) == 1L && pixel_size > 0)) {
      stop("`pixel_size` must be a single positive number (meters).")
    }
    yr <- W$yrange; xr <- W$xrange
    ny <- max(1L, round((yr[2] - yr[1]) / pixel_size))
    nx <- max(1L, round((xr[2] - xr[1]) / pixel_size))
    dimyx <- c(ny, nx)
  } else {
    if (!(is.numeric(dimyx) && length(dimyx) %in% c(1L,2L))) stop("`dimyx` must be length 1 or 2 numeric.")
    if (length(dimyx) == 1L) dimyx <- rep(as.integer(dimyx), 2L)
    dimyx <- pmax(1L, as.integer(dimyx))
  }

  # --- KDE ---------------------------------------------------------------------
  kde <- spatstat.explore::density.ppp(
    X,
    sigma   = radius,
    weights = weights,
    dimyx   = dimyx,
    at      = "pixels",
    edge    = TRUE
  )

  # --- to SpatRaster, classify, mask ------------------------------------------
  r_kde <- terra::rast(kde)
  terra::crs(r_kde) <- sf::st_crs(area_m$shape)$wkt

  re_mat <- reclass_matrix(r_kde, n_classes = n_classes,
                           exclude_lowest = exclude_lowest, lowest_prop = lowest_prop)

  r_cls <- terra::classify(r_kde, re_mat, include.lowest = TRUE)
  r_cls <- terra::mask(r_cls, terra::vect(area_m$shape))
  names(r_cls) <- "Rating"

  # --- to polygons if needed ---------------------------------------------------
  if (output_layer_type != "raster") {
    v_cls <- sf::st_as_sf(terra::as.polygons(r_cls, dissolve = TRUE))
    v_cls <- sf::st_set_crs(v_cls, sf::st_crs(area_m$shape))
  }

  # --- reproject outputs if requested -----------------------------------------
  if (return_crs == "4326") {
    if (output_layer_type %in% c("raster","both")) {
      # nearest neighbor to preserve class labels
      r_cls <- terra::project(r_cls, "EPSG:4326", method = "near")
    }
    if (output_layer_type %in% c("shp","both")) {
      v_cls <- sf::st_transform(v_cls, 4326)
    }
  }

  # --- return ------------------------------------------------------------------
  if (output_layer_type == "shp") {
    return(v_cls)
  } else if (output_layer_type == "raster") {
    return(r_cls)
  } else {
    return(list(raster = r_cls, shp = v_cls))
  }
}


################################################################################
############################get_overlap_kernel##################################
################################################################################

#' Compute overlap hotspots from two reclassified rasters
#'
#' Combines two class-coded rasters (e.g., exposure & consequence) into an overlap/risk
#' surface, then reclassifies to `out_classes` bins. Inputs are assumed to have integer
#' classes in `1…n_classes` (0 or NA treated as background).
#'
#' @param x,y Single-layer `terra::SpatRaster` with the same CRS; values are class codes.
#' @param n_classes Integer: number of classes in the inputs (default 3).
#' @param out_classes Integer: number of classes in the output (default = `n_classes`).
#' @param method Combination rule: `"product"` (default), `"sum"`, `"geom_mean"`, or `"max"`.
#' @param resample_method Resampling to align `y` to `x` when grids differ. For class rasters,
#'   the default `"near"` preserves class labels.
#' @param output_layer_type One of `"shp"`, `"raster"`, or `"both"`. Default `"shp"`.
#' @param quiet Logical; suppress informative messages. Default `TRUE`.
#' @return An `sf` polygon layer, a `SpatRaster`, or a list with both.
#' @importFrom terra same.crs compareGeom resample classify as.polygons global rast values mask crs overlay
#' @importFrom sf st_as_sf st_set_crs
#' @examples
#' species  <- data.frame(long = rnorm(80, 0, 10),  lat = rnorm(80, 0, 10))
#' stressor <- data.frame(long = rnorm(100, 0, 10), lat = rnorm(100, 0, 10))
#' kde_spe <- get_class_kernel(species,  output_layer_type = "raster")
#' kde_str <- get_class_kernel(stressor, output_layer_type = "raster")
#' # Ensure same grid (nearest to preserve class labels)
#' kde_str <- terra::project(kde_str, kde_spe, method = "near")
#' overlap <- get_overlap_kernel(kde_spe, kde_str, method = "product", output_layer_type = "raster")
#' terra::plot(overlap)
#' @export
get_overlap_kernel <- function(
    x, y,
    n_classes = 3,
    out_classes = n_classes,
    method = c("product","sum","geom_mean","max"),
    resample_method = "near",
    output_layer_type = c("shp","raster","both"),
    quiet = TRUE
) {
  method            <- match.arg(method)
  output_layer_type <- match.arg(output_layer_type)

  # --- basic checks ------------------------------------------------------------
  if (!inherits(x, "SpatRaster") || !inherits(y, "SpatRaster")) {
    stop("`x` and `y` must be terra::SpatRaster.")
  }
  if (terra::nlyr(x) != 1L || terra::nlyr(y) != 1L) {
    stop("`x` and `y` must be single-layer rasters.")
  }
  if (!terra::same.crs(x, y)) stop("Input rasters must have the same CRS.")

  # Align y to x grid if needed (use nearest to keep integer class codes)
  if (!terra::compareGeom(x, y, stopOnError = FALSE)) {
    if (!quiet) message("Aligning `y` to `x` with terra::resample(..., method = '", resample_method, "').")
    y <- terra::resample(y, x, method = resample_method)
  }

  # Treat 0 as background (set to NA); keep NA as NA
  x <- terra::ifel(x <= 0, NA, x)
  y <- terra::ifel(y <= 0, NA, y)

  # Validate class ranges (after resampling/cleaning)
  mmx <- terra::global(x, c("min","max"), na.rm = TRUE)
  mmy <- terra::global(y, c("min","max"), na.rm = TRUE)
  minx <- as.numeric(mmx[1,"min"]); maxx <- as.numeric(mmx[1,"max"])
  miny <- as.numeric(mmy[1,"min"]); maxy <- as.numeric(mmy[1,"max"])
  if (is.finite(minx) && (minx < 1 || maxx > n_classes)) {
    stop("`x` values must be within 1…n_classes (or NA).")
  }
  if (is.finite(miny) && (miny < 1 || maxy > n_classes)) {
    stop("`y` values must be within 1…n_classes (or NA).")
  }

  # --- combine -----------------------------------------------------------------
  # Define continuous "score" and theoretical min/max for scaling
  if (method == "product") {
    score    <- x * y
    s_min    <- 1
    s_max    <- n_classes^2
  } else if (method == "sum") {
    score    <- x + y
    s_min    <- 2
    s_max    <- 2 * n_classes
  } else if (method == "geom_mean") {
    score    <- sqrt(x * y)
    s_min    <- 1
    s_max    <- n_classes
  } else if (method == "max") {
    score    <- terra::overlay(x, y, fun = function(a, b) pmax(a, b))
    s_min    <- 1
    s_max    <- n_classes
  }

  # Normalize to [0,1] using theoretical bounds, then reclass to out_classes bins
  # (Using theory avoids data-dependent shrinkage.)
  norm <- (score - s_min) / (s_max - s_min)

  brks <- seq(0, 1, length.out = out_classes + 1)
  mat  <- cbind(from = brks[-length(brks)], to = brks[-1], class = seq_len(out_classes))
  result_reclass <- terra::classify(norm, rcl = mat, include.lowest = TRUE)
  names(result_reclass) <- "Rating"
  terra::crs(result_reclass) <- terra::crs(x)

  # --- vectorize if requested --------------------------------------------------
  if (output_layer_type != "raster") {
    v <- sf::st_as_sf(terra::as.polygons(result_reclass, dissolve = TRUE))
    v <- sf::st_set_crs(v, terra::crs(result_reclass))
  }

  # --- return ------------------------------------------------------------------
  if (output_layer_type == "shp") {
    return(v)
  } else if (output_layer_type == "raster") {
    return(result_reclass)
  } else {
    return(list(raster = result_reclass, shp = v))
  }
}


################################################################################
##################################risa_prep#####################################
################################################################################

#' Prepare and generate maps for risk assessment analysis
#'
#' Generates kernel density maps for species and stressors, plus overlap maps,
#' within an area of interest (auto or user-supplied).
#'
#' @param x Species input: `sf`, data.frame, or a list of `sf`. If a data.frame/sf,
#'   you can split it by `group_x` into multiple species layers.
#' @param y Stressor input: `sf`, data.frame, or a list of `sf`. If a data.frame/sf,
#'   you can split it by `group_y` into multiple stressor layers.
#' @param area Optional AOI polygon (`sf`) or `bbox`/data.frame; if `NULL`, it's computed.
#' @param n_classes Integer number of classes for kernel reclassification (default 3).
#' @param output_layer_type One of `"shp"`, `"raster"`, or `"both"` (default `"both"`).
#' @param radius KDE bandwidth (projected units). If `NULL`, uses `radius_method`.
#' @param radius_method One of `"nndist"` (default), `"ppl"`, or `"fixed"`.
#' @param group_x,group_y Optional grouping columns to split `x`/`y` when they aren't lists.
#' @param group_size_x,group_size_y Optional numeric weight columns for KDE.
#' @param pixel_size Optional pixel size (meters) for the KDE grid; else `dimyx` is used.
#' @param dimyx Grid size (ny, nx) for KDE when `pixel_size` is `NULL` (default c(512,512)).
#' @param exclude_lowest,lowest_prop Passed to `reclass_matrix()` (defaults TRUE, 0.05).
#' @param area_strategy `"stressor"` (default), `"species"`, or `"union"` for auto AOI.
#' @param area_type `"convex_hull"` (default) or `"bbox"` for auto AOI.
#' @param area_buffer_frac Fractional expansion for auto AOI (default 0.5).
#' @param return_crs `"metric"` (default) or `"4326"` to reproject outputs to WGS84.
#' @param overlap_method Combination rule for overlap: one of `"product"`, `"sum"`,
#'   `"geom_mean"`, or `"max"` (default `"product"`). Passed to `get_overlap_kernel()`.
#' @param quiet Suppress messages (default TRUE).
#' @return A list with: `species_distributions`, `stressor_distributions`,
#'   `species_kernel_maps`, `stressor_kernel_maps`, `overlap_maps`, and `area_of_interest`.
#' @importFrom sf st_as_sfc st_transform
#' @importFrom terra project
#' @examples
#' # Creating test data
#' spp_df <- rbind(data.frame(long = rnorm(80, 0, 10),
                           #' lat = rnorm(80, 0, 10), species = "sp1"),
                #' data.frame(long = rnorm(60, 0, 10),
                           #' lat = rnorm(60, 0, 10), species = "sp2"))
#' str_df <- rbind(data.frame(long = rnorm(100, 0, 10),
                           #' lat = rnorm(100, 0, 10), stressor = "trawling"),
                #' data.frame(long = rnorm(50, 0, 10),
                           #' lat = rnorm(100, 0, 5), stressor = "gillnet"))
#' # Create kernel maps of species and stressor distributions and overlap maps
#' risa_maps <- risa_prep(spp_df, str_df)
#'
#' # Species and Stressor distributions
#' dev.off()
#' par(mfrow = c(2,2))
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species 1")
#' plot(risa_maps$species_kernel_maps$sp1$shp, add = TRUE)
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species 2")
#' plot(risa_maps$species_kernel_maps$sp2$shp, add = TRUE)
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Gillnet")
#' plot(risa_maps$stressor_kernel_maps$gillnet$shp, add = TRUE, col=c("lightgreen", "green", "darkgreen"))
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Trawling")
#' plot(risa_maps$stressor_kernel_maps$trawling$shp, add = TRUE, col=c("lightgreen", "green", "darkgreen"))
#' # Overlap maps
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species1 vs Gillnet")
#' plot(risa_maps$overlap_maps$sp1$gillnet$shp, add = TRUE, col=c("yellow", "orange", "red"))
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species2 vs Gillnet")
#' plot(risa_maps$overlap_maps$sp2$gillnet$shp, add = TRUE, col=c("yellow", "orange", "red"))
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species1 vs Trawling")
#' plot(risa_maps$overlap_maps$sp1$trawling$shp, add = TRUE, col=c("yellow", "orange", "red"))
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species2 vs Trawling")
#' plot(risa_maps$overlap_maps$sp2$trawling$shp, add = TRUE, col=c("yellow", "orange", "red"))
#' @export
risa_prep <- function(
    x, y,
    area = NULL,
    n_classes = 3,
    output_layer_type = c("both","shp","raster"),
    radius = NULL,
    radius_method = c("nndist","ppl","fixed"),
    group_x = NULL,
    group_y = NULL,
    group_size_x = NULL,
    group_size_y = NULL,
    pixel_size = NULL,
    dimyx = c(512,512),
    exclude_lowest = TRUE,
    lowest_prop = 0.05,
    area_strategy = c("stressor","species","union"),
    area_type = c("convex_hull","bbox"),
    area_buffer_frac = 0.5,
    return_crs = c("metric","4326"),
    overlap_method = c("product","sum","geom_mean","max"),
    quiet = TRUE
) {
  output_layer_type <- match.arg(output_layer_type)
  radius_method     <- match.arg(radius_method)
  area_strategy     <- match.arg(area_strategy)
  area_type         <- match.arg(area_type)
  return_crs        <- match.arg(return_crs)
  overlap_method    <- match.arg(overlap_method)

  # -- helper: normalize to a named list of sf ---------------------------------
  as_sf_list <- function(obj, group = NULL, label_prefix = "layer") {
    # 1) sf: wrap as single-element list
    if (inherits(obj, "sf")) {
      nm <- if (!is.null(names(obj)) && any(nzchar(names(obj)))) names(obj)[1] else label_prefix
      return(setNames(list(obj), nm))
    }

    # 2) data.frame: split or convert
    if (is.data.frame(obj)) {
      if (!is.null(group) && group %in% names(obj)) {
        # split by explicit group col
        sp <- split(obj, obj[[group]], drop = TRUE)
        lst <- lapply(sp, df_to_shp)
        return(lst)
      } else {
        # use your df_to_list() (splits by 3rd column if present)
        return(df_to_list(obj))
      }
    }

    # 3) generic list: each element must be sf or data.frame
    if (is.list(obj)) {
      nms <- names(obj)
      lst <- lapply(obj, function(el) {
        if (inherits(el, "sf")) return(el)
        if (is.data.frame(el))  return(df_to_shp(el))
        stop("List elements must be `sf` or data.frame.")
      })
      if (is.null(nms) || any(!nzchar(nms))) {
        names(lst) <- if (!is.null(nms)) {
          # fill any empties
          nms[nchar(nms) == 0] <- paste0(label_prefix, seq_len(sum(nchar(nms) == 0)))
          nms
        } else {
          paste0(label_prefix, seq_along(lst))
        }
      }
      return(lst)
    }

    stop("Input must be an `sf`, data.frame, or a list of those.")
  }

  spp_list <- as_sf_list(x, group = group_x, label_prefix = "sp")
  str_list <- as_sf_list(y, group = group_y, label_prefix = "stressor")

  # -- choose / build AOI -------------------------------------------------------
  if (is.null(area)) {
    if (!quiet) message("No area provided; creating AOI from ", area_strategy, " layers (", area_type, ").")
    src <- switch(area_strategy,
                  "stressor" = merge_shp(str_list),
                  "species"  = merge_shp(spp_list),
                  "union"    = merge_shp(c(spp_list, str_list)))
    area <- create_area(src, area_type = area_type, buffer_frac = area_buffer_frac, quiet = quiet)
  } else if (inherits(area, "bbox")) {
    area <- sf::st_as_sfc(area)
  } else if (inherits(area, "data.frame")) {
    area <- df_to_shp(area)
  } else if (!inherits(area, "sf")) {
    stop("`area` must be NULL, bbox, data.frame, or sf.")
  }

  # -- kernels helper -----------------------------------------------------------
  build_kernels <- function(lst, group_size = NULL, ncls = n_classes) {
    lapply(lst, function(item) {
      get_class_kernel(
        item,
        area               = area,
        n_classes          = ncls,
        output_layer_type  = "both",
        radius             = radius,
        radius_method      = radius_method,
        group_size         = group_size,
        pixel_size         = pixel_size,
        dimyx              = dimyx,
        exclude_lowest     = exclude_lowest,
        lowest_prop        = lowest_prop,
        return_crs         = return_crs,
        quiet              = quiet
      )
    })
  }

  # Distributions (presence maps): 1 class; Kernels: n_classes
  spp_distribution_list <- build_kernels(spp_list, group_size = group_size_x, ncls = 1L)
  str_distribution_list <- build_kernels(str_list, group_size = group_size_y, ncls = 1L)
  spp_kernel_list       <- build_kernels(spp_list, group_size = group_size_x, ncls = n_classes)
  stressor_kernel_list  <- build_kernels(str_list, group_size = group_size_y, ncls = n_classes)

  # AOI CRS according to return_crs
  if (return_crs == "4326") {
    area <- sf::st_transform(area, 4326)
  }

  # -- overlap maps -------------------------------------------------------------
  if (!quiet) message("Generating overlap maps...")
  overlap_maps_list <- lapply(spp_kernel_list, function(sp) {
    lapply(stressor_kernel_list, function(st) {
      ol <- get_overlap_kernel(
        sp$raster, st$raster,
        n_classes          = n_classes,
        method             = overlap_method,
        output_layer_type  = output_layer_type,
        quiet              = quiet
      )
      if (return_crs == "4326") {
        if (inherits(ol, "SpatRaster")) {
          ol <- terra::project(ol, "EPSG:4326", method = "near")
        } else if (inherits(ol, "sf")) {
          ol <- sf::st_transform(ol, 4326)
        } else if (is.list(ol) && all(c("raster","shp") %in% names(ol))) {
          ol$raster <- terra::project(ol$raster, "EPSG:4326", method = "near")
          ol$shp    <- sf::st_transform(ol$shp, 4326)
        }
      }
      ol
    })
  })

  # -- assemble -----------------------------------------------------------------
  out <- list(
    species_distributions  = spp_distribution_list,
    stressor_distributions = str_distribution_list,
    species_kernel_maps    = spp_kernel_list,
    stressor_kernel_maps   = stressor_kernel_list,
    overlap_maps           = overlap_maps_list,
    area_of_interest       = area
  )
  class(out) <- c("risaMaps", class(out))
  out
}

################################################################################
####################################hra#########################################
################################################################################

#' Habitat Risk Assessment (HRA): single species or ecosystem
#'
#' Computes Exposure/Consequence and Risk for a single species (2-level `raster_list`)
#' or multiple species (3-level `raster_list`), and returns per-stressor maps plus
#' cumulative risk. Includes a `summary_stats` data.frame.
#'
#' @param raster_list For single species: named list of lists \[stressor -> attribute rasters\].
#'   For ecosystem: named list \[species -> (stressor -> attribute rasters)\].
#'   Each attribute raster is a single-layer `terra::SpatRaster` on the same grid as the species distribution.
#' @param species_distr For single species: `SpatRaster` presence mask (non-NA = presence).
#'   For ecosystem: **list** of `SpatRaster` (one per species), names must match `raster_list`.
#' @param criteria For single species: a data.frame with columns
#'   `STRESSOR, ATTRIBUTES, RATING, DQ, WEIGHT, E/C`.
#'   For ecosystem: **list** of such data.frames, names must match species.
#' @param equation c("euclidean","multiplicative")
#' @param r_max Integer (1–10), max criteria score for thresholds.
#' @param n_overlap Optional # of stressors used in cumulative-risk thresholds.
#'   Defaults to the number of unique stressors in `criteria` (single) or across species (ecosystem).
#' @param output_decimal_crs Logical; if TRUE, reprojects rasters in the output to EPSG:4326
#'   (nearest-neighbor to preserve classes). Summary stats are computed pre-reprojection.
#' @return A list of class `"risaHRA"`. In single-species mode:
#'   per-stressor sublists + `total_raw`, `total`, and `summary_stats`.
#'   In ecosystem mode: one such sublist per species plus `ecosys_risk_raw`,
#'   `ecosys_risk_classified`, and combined `summary_stats` (with `SPECIES`).
#' @importFrom terra ifel mosaic mask classify project compareGeom global freq nlyr resample
#' @export
hra <- function(
    raster_list, species_distr, criteria,
    equation = c("euclidean", "multiplicative"),
    r_max = 3, n_overlap = NULL, output_decimal_crs = FALSE
) {
  depth <- list_depth_base(raster_list)
  if (!(depth %in% c(2L, 3L))) {
    stop("'raster_list' must have depth 2 (single species) or 3 (ecosystem).")
  }

  # -------------------------- internal helpers --------------------------------
  .check_criteria <- function(df_or_list) {
    req <- c("STRESSOR","ATTRIBUTES","RATING","DQ","WEIGHT","E/C")
    if (is.data.frame(df_or_list)) {
      df <- df_or_list
      if (names(df)[1] != "STRESSOR") {
        message("'criteria' not in expected long format. Trying criteria_reshape()...")
        df <- criteria_reshape(df)
        if (length(df) > 1L) stop("Multiple species detected in criteria; pass a list for ecosystem mode.")
        df <- df[[1]]
      }
      if (!all(req %in% names(df))) stop("Criteria must contain: ", paste(req, collapse = ", "))
      return(df)
    } else if (is.list(df_or_list)) {
      out <- lapply(df_or_list, function(d) {
        if (!is.data.frame(d)) stop("Each criteria element must be a data.frame.")
        if (names(d)[1] != "STRESSOR") {
          d2 <- criteria_reshape(d)
          if (length(d2) != 1L) stop("Each criteria element must correspond to exactly one species.")
          d <- d2[[1]]
        }
        if (!all(req %in% names(d))) stop("Each criteria table must contain: ", paste(req, collapse = ", "))
        d
      })
      return(out)
    } else {
      stop("'criteria' must be a data.frame (single) or a list of data.frames (ecosystem).")
    }
  }

  .compute_summary_stats_single <- function(rr, total_risk_raw, total_class) {
    stressor_names <- setdiff(names(rr), c("total_raw","total","summary_stats"))
    rows <- lapply(stressor_names, function(st) {
      E  <- rr[[st]]$E_criteria
      C  <- rr[[st]]$C_criteria
      Rr <- rr[[st]]$Risk_map_raw
      Rc <- rr[[st]]$Risk_map

      gE <- terra::global(E,  c("min","max","mean"), na.rm = TRUE)
      gC <- terra::global(C,  c("min","max","mean"), na.rm = TRUE)
      gR <- terra::global(Rr, c("min","max","mean"), na.rm = TRUE)

      fq <- terra::freq(Rc, useNA = "no")
      cnt <- c(`0`=0,`1`=0,`2`=0,`3`=0)
      if (!is.null(fq) && nrow(fq)) {
        vv <- as.character(fq[,"value"]); cnt[vv] <- fq[,"count"]
      }
      tot <- sum(cnt)

      data.frame(
        STRESSOR = st,
        E_min   = as.numeric(gE[1,"min"]),
        E_max   = as.numeric(gE[1,"max"]),
        E_mean  = as.numeric(gE[1,"mean"]),
        C_min   = as.numeric(gC[1,"min"]),
        C_max   = as.numeric(gC[1,"max"]),
        C_mean  = as.numeric(gC[1,"mean"]),
        R_min   = as.numeric(gR[1,"min"]),
        R_max   = as.numeric(gR[1,"max"]),
        R_mean  = as.numeric(gR[1,"mean"]),
        `R%high`   = if (tot) 100 * cnt["3"]/tot else NA_real_,
        `R%medium` = if (tot) 100 * cnt["2"]/tot else NA_real_,
        `R%low`    = if (tot) 100 * cnt["1"]/tot else NA_real_,
        `R%None`   = if (tot) 100 * cnt["0"]/tot else NA_real_
      )
    })
    per_stressor <- do.call(rbind, rows)

    gRt <- terra::global(total_risk_raw, c("min","max","mean"), na.rm = TRUE)
    fqt <- terra::freq(total_class, useNA = "no")
    cntt <- c(`0`=0,`1`=0,`2`=0,`3`=0)
    if (!is.null(fqt) && nrow(fqt)) {
      vv <- as.character(fqt[,"value"]); cntt[vv] <- fqt[,"count"]
    }
    tot <- sum(cntt)

    overall <- data.frame(
      STRESSOR = "(FROM ALL STRESSORS)",
      E_min   = min(per_stressor$E_min,   na.rm = TRUE),
      E_max   = max(per_stressor$E_max,   na.rm = TRUE),
      E_mean  = mean(per_stressor$E_mean, na.rm = TRUE),
      C_min   = min(per_stressor$C_min,   na.rm = TRUE),
      C_max   = max(per_stressor$C_max,   na.rm = TRUE),
      C_mean  = mean(per_stressor$C_mean, na.rm = TRUE),
      R_min   = as.numeric(gRt[1,"min"]),
      R_max   = as.numeric(gRt[1,"max"]),
      R_mean  = as.numeric(gRt[1,"mean"]),
      `R%high`   = if (tot) 100 * cntt["3"]/tot else NA_real_,
      `R%medium` = if (tot) 100 * cntt["2"]/tot else NA_real_,
      `R%low`    = if (tot) 100 * cntt["1"]/tot else NA_real_,
      `R%None`   = if (tot) 100 * cntt["0"]/tot else NA_real_
    )

    rbind(overall, per_stressor)
  }

  .hra_single <- function(rlist, sp_distr, crit, equation, r_max, n_overlap, output_decimal_crs) {
    # validate basics
    if (list_depth_base(rlist) != 2L) stop("Single-species mode requires depth-2 raster_list.")
    if (!inherits(sp_distr, "SpatRaster")) stop("'species_distr' must be a SpatRaster.")
    crit <- .check_criteria(crit)

    equation <- match.arg(equation)
    stressors <- unique(crit$STRESSOR)
    stressors <- stressors[!(is.na(stressors) | stressors %in% c("", "NA"))]
    if (is.null(n_overlap)) n_overlap <- length(stressors)
    if (!setequal(names(rlist), stressors)) {
      stop("Top-level names in 'raster_list' must match criteria STRESSOR values.")
    }
    if (!is.numeric(r_max) || r_max < 1 || r_max > 10) stop("'r_max' must be in 1..10.")

    sp_presence    <- terra::ifel(!is.na(sp_distr), 1, NA)
    sp_distr_zeros <- terra::ifel(!is.na(sp_distr), 0, NA)
    zero_r         <- sp_distr * 0

    res <- list()
    for (stressor in stressors) {
      crit_stressor <- crit[is.na(crit$STRESSOR) | crit$STRESSOR %in% c("", "NA", stressor), ]
      C_df <- crit_stressor[crit_stressor$`E/C` == "C", , drop = FALSE]
      E_df <- crit_stressor[crit_stressor$`E/C` == "E", , drop = FALSE]
      if (nrow(E_df) == 0L || nrow(C_df) == 0L) {
        stop("Both E and C blocks must have at least one criterion for '", stressor, "'.")
      }
      E_df$DQ <- as.numeric(E_df$DQ); E_df$WEIGHT <- as.numeric(E_df$WEIGHT)
      C_df$DQ <- as.numeric(C_df$DQ); C_df$WEIGHT <- as.numeric(C_df$WEIGHT)
      if (any(E_df$DQ<=0|E_df$WEIGHT<=0) || any(C_df$DQ<=0|C_df$WEIGHT<=0)) stop("DQ/WEIGHT must be > 0.")

      E_df$RATING <- suppressWarnings(as.numeric(E_df$RATING))
      C_df$RATING <- suppressWarnings(as.numeric(C_df$RATING))

      E_const  <- E_df[!is.na(E_df$RATING), , drop = FALSE]
      C_const  <- C_df[!is.na(C_df$RATING), , drop = FALSE]
      E_mapped <- E_df[ is.na(E_df$RATING), , drop = FALSE]
      C_mapped <- C_df[ is.na(C_df$RATING), , drop = FALSE]

      E_numer_const <- if (nrow(E_const)) sum(E_const$RATING / (E_const$DQ * E_const$WEIGHT)) else 0
      C_numer_const <- if (nrow(C_const)) sum(C_const$RATING / (C_const$DQ * C_const$WEIGHT)) else 0

      sum_weighted <- function(df_map) {
        if (!nrow(df_map)) return(zero_r)
        parts <- vector("list", nrow(df_map))
        for (i in seq_len(nrow(df_map))) {
          att <- df_map$ATTRIBUTES[i]
          r   <- rlist[[stressor]][[att]]
          if (is.null(r) || !inherits(r, "SpatRaster")) {
            stop("Missing/invalid raster for stressor '", stressor, "', attribute '", att, "'.")
          }
          if (!terra::compareGeom(sp_distr, r, stopOnError = FALSE)) {
            stop("Raster grid for '", stressor, " / ", att, "' does not match 'species_distr'.")
          }
          parts[[i]] <- r / (df_map$DQ[i] * df_map$WEIGHT[i])
        }
        Reduce(`+`, parts)
      }

      E_numer_rast <- sum_weighted(E_mapped)
      C_numer_rast <- sum_weighted(C_mapped)

      E_denom <- sum(1 / (E_df$DQ * E_df$WEIGHT))
      C_denom <- sum(1 / (C_df$DQ * C_df$WEIGHT))

      E_score_raster <- (E_numer_const + E_numer_rast) / E_denom
      C_score_raster <- (C_numer_const + C_numer_rast) / C_denom

      E_map <- terra::mosaic(E_score_raster, sp_distr_zeros, fun = "first")
      C_map <- terra::mosaic(terra::mask(C_score_raster, sp_distr), sp_distr_zeros, fun = "first")

      if (equation == "multiplicative") {
        m_jkl    <- r_max^2
        risk_raw <- (C_score_raster * E_score_raster) * sp_presence
      } else {
        m_jkl    <- sqrt(2 * (r_max - 1)^2)
        risk_raw <- sqrt((E_score_raster - 1)^2 + (C_score_raster - 1)^2) * sp_presence
      }
      risk_raw <- terra::mosaic(risk_raw, sp_distr_zeros, fun = "first")

      risk_cls <- terra::ifel(
        risk_raw == 0, 0,
        terra::ifel(risk_raw < (1/3)*m_jkl, 1,
                    terra::ifel(risk_raw < (2/3)*m_jkl, 2, 3))
      )

      res[[stressor]] <- list(
        E_criteria   = E_map,
        C_criteria   = C_map,
        Risk_map_raw = risk_raw,
        Risk_map     = risk_cls
      )
    }

    total_raw <- Reduce(`+`, lapply(res, function(x) x$Risk_map_raw))
    total_cls <- terra::ifel(
      total_raw == 0, 0,
      terra::ifel(total_raw < (1/3)*m_jkl*n_overlap, 1,
                  terra::ifel(total_raw < (2/3)*m_jkl*n_overlap, 2, 3))
    )

    # summary (pre-reprojection)
    res$total_raw     <- total_raw
    res$total         <- total_cls
    res$summary_stats <- .compute_summary_stats_single(res, total_raw, total_cls)

    if (isTRUE(output_decimal_crs)) {
      for (st in names(res)) {
        if (is.list(res[[st]])) {
          for (nm in c("E_criteria","C_criteria","Risk_map_raw","Risk_map")) {
            res[[st]][[nm]] <- convert_to_decimal_degrees(res[[st]][[nm]])
          }
        }
      }
      res$total_raw <- convert_to_decimal_degrees(total_raw)
      res$total     <- convert_to_decimal_degrees(total_cls)
    }

    class(res) <- c("risaHRA", class(res))
    res
  }

  # -------------------------- single vs ecosystem -----------------------------
  equation <- match.arg(equation)

  if (depth == 2L) {
    # single species
    return(.hra_single(
      rlist = raster_list,
      sp_distr = species_distr,
      crit = .check_criteria(criteria),
      equation = equation,
      r_max = r_max,
      n_overlap = n_overlap,
      output_decimal_crs = output_decimal_crs
    ))
  }

  # ecosystem mode (depth 3)
  if (!is.list(species_distr)) stop("'species_distr' must be a list of SpatRasters in ecosystem mode.")
  if (!is.list(criteria))      stop("'criteria' must be a list of data.frames in ecosystem mode.")

  species <- names(raster_list)
  if (is.null(species) || any(!nzchar(species))) stop("Top-level 'raster_list' must be a named list of species.")
  if (!setequal(species, names(species_distr))) stop("Names of 'species_distr' must match 'raster_list' species.")
  if (!setequal(species, names(criteria)))      stop("Names of 'criteria' must match 'raster_list' species.")

  # Determine default n_overlap from the union of stressors across species
  if (is.null(n_overlap)) {
    all_stressors <- unique(unlist(lapply(criteria, function(df) {
      df2 <- .check_criteria(df); df2$STRESSOR
    })))
    all_stressors <- all_stressors[!(is.na(all_stressors) | all_stressors %in% c("", "NA"))]
    n_overlap <- length(all_stressors)
  }

  results <- list()
  for (sp in species) {
    results[[sp]] <- .hra_single(
      rlist = raster_list[[sp]],
      sp_distr = species_distr[[sp]],
      crit = .check_criteria(criteria[[sp]]),
      equation = equation,
      r_max = r_max,
      n_overlap = n_overlap,
      output_decimal_crs = FALSE  # postpone proj until after eco stats
    )
  }

  # Ecosystem mask (union of species presence)
  presences <- lapply(species_distr, function(d) terra::ifel(!is.na(d), 1, 0))
  sum_pres  <- Reduce(`+`, presences)
  eco_mask  <- terra::ifel(sum_pres > 0, 1, NA)

  # Ecosystem risk: average of species totals (NA treated as 0 inside mask)
  eco_raw <- NULL
  for (sp in species) {
    r <- results[[sp]]$total_raw
    r <- terra::ifel(is.na(r), 0, r)
    eco_raw <- if (is.null(eco_raw)) r else eco_raw + r
  }
  eco_raw <- eco_raw / length(species)
  eco_raw <- terra::mosaic(eco_raw, terra::ifel(!is.na(eco_mask), 0, NA), fun = "first")

  # Classification thresholds use same m_jkl & n_overlap
  if (equation == "multiplicative") {
    m_jkl <- r_max^2
  } else {
    m_jkl <- sqrt(2 * (r_max - 1)^2)
  }
  eco_cls <- terra::ifel(
    eco_raw == 0, 0,
    terra::ifel(eco_raw < (1/3)*m_jkl*n_overlap, 1,
                terra::ifel(eco_raw < (2/3)*m_jkl*n_overlap, 2, 3))
  )

  # Build combined summary_stats (pre-reprojection)
  per_species_stats <- do.call(rbind, lapply(names(results), function(sp) {
    cbind.data.frame(SPECIES = sp, results[[sp]]$summary_stats, row.names = NULL)
  }))

  # Ecosystem summary row
  gEco <- terra::global(eco_raw, c("min","max","mean"), na.rm = TRUE)
  fqE  <- terra::freq(eco_cls, useNA = "no")
  cntE <- c(`0`=0,`1`=0,`2`=0,`3`=0)
  if (!is.null(fqE) && nrow(fqE)) {
    vv <- as.character(fqE[,"value"]); cntE[vv] <- fqE[,"count"]
  }
  totE <- sum(cntE)
  eco_row <- data.frame(
    SPECIES = "ECOSYSTEM",
    STRESSOR = "(FROM ALL STRESSORS)",
    E_min = NA_real_, E_max = NA_real_, E_mean = NA_real_,
    C_min = NA_real_, C_max = NA_real_, C_mean = NA_real_,
    R_min = as.numeric(gEco[1,"min"]),
    R_max = as.numeric(gEco[1,"max"]),
    R_mean = as.numeric(gEco[1,"mean"]),
    `R%high`   = if (totE) 100 * cntE["3"]/totE else NA_real_,
    `R%medium` = if (totE) 100 * cntE["2"]/totE else NA_real_,
    `R%low`    = if (totE) 100 * cntE["1"]/totE else NA_real_,
    `R%None`   = if (totE) 100 * cntE["0"]/totE else NA_real_
  )

  summary_stats <- rbind(eco_row, per_species_stats)

  # Attach ecosystem layers and optionally reproject all
  results$ecosys_risk_raw        <- eco_raw
  results$ecosys_risk_classified <- eco_cls
  results$summary_stats          <- summary_stats

  if (isTRUE(output_decimal_crs)) {
    # project species sub-results
    for (sp in species) {
      for (nm in names(results[[sp]])) {
        if (is.list(results[[sp]][[nm]])) {
          for (k in c("E_criteria","C_criteria","Risk_map_raw","Risk_map")) {
            results[[sp]][[nm]][[k]] <- convert_to_decimal_degrees(results[[sp]][[nm]][[k]])
          }
        } else if (nm %in% c("total_raw","total")) {
          results[[sp]][[nm]] <- convert_to_decimal_degrees(results[[sp]][[nm]])
        }
      }
    }
    # project ecosystem layers
    results$ecosys_risk_raw        <- convert_to_decimal_degrees(results$ecosys_risk_raw)
    results$ecosys_risk_classified <- convert_to_decimal_degrees(results$ecosys_risk_classified)
  }

  class(results) <- c("risaHRA", class(results))
  results
}


################################################################################
############################criteria_reshape####################################
################################################################################

#' Reshape a criteria table into a list of data frames respective to each species
#'
#' Splits a wide Habitat Risk Assessment criteria table (first column with
#' attribute/stressor labels, triplets of columns per species: RATING/DQ/WEIGHT,
#' and a column with E/C type) into a list of long-format data frames—one per
#' species—with columns: STRESSOR, ATTRIBUTES, RATING, DQ, WEIGHT, `E/C`.
#'
#' The function is robust to internal header rows and blank-row separators.
#' It detects the criteria-type column as the rightmost column whose values
#' are in {E, C} (case-insensitive), and groups the species triplets from the
#' columns between the first and the criteria-type column.
#'
#' @param x A `data.frame`. Expected layout:
#'   - Column 1: attribute/stressor labels (character), containing blank rows
#'     that separate blocks.
#'   - Columns 2..(2+3k-1): for each species, three columns (RATING, DQ, WEIGHT).
#'   - One column (typically last) with `E`/`C` entries (criteria type).
#'     This is auto-detected.
#'   Internal header rows are allowed and will be dropped.
#' @return A named `list` of `data.frame` objects, one per species, each with
#'   columns `STRESSOR, ATTRIBUTES, RATING, DQ, WEIGHT, E/C`.
#' @examples
#' # crit_list <- criteria_reshape(df)
#' # str(crit_list$sp1)
#' @export
criteria_reshape <- function(x) {
  if (!is.data.frame(x) || ncol(x) < 6L) {
    stop("`x` must be a data.frame with at least 6 columns: 1 (ATTRIBUTES) + triplets per species + 1 (E/C).")
  }

  # --- detect criteria-type column (E/C) ---------------------------------------
  is_ec_col <- function(v) {
    u <- unique(na.omit(trimws(toupper(as.character(v)))))
    length(u) > 0L && all(u %in% c("E", "C"))
  }
  crit_col <- NULL
  for (j in rev(seq_len(ncol(x)))) {
    if (is_ec_col(x[[j]])) { crit_col <- j; break }
  }
  if (is.null(crit_col)) stop("Could not find a column with only 'E'/'C' values (criteria type).")

  # --- determine number of species triplets -----------------------------------
  mid_cols <- (crit_col - 1L) - 1L   # columns strictly between col 1 and E/C col
  if (mid_cols < 3L) stop("Not enough species columns between the first column and the E/C column.")
  if (mid_cols %% 3L != 0L) warning("Species columns between ATTRIBUTES and E/C are not a multiple of 3; extra columns will be ignored.")
  sp_n <- floor(mid_cols / 3L)
  if (sp_n < 1L) stop("No complete species triplets (RATING/DQ/WEIGHT) were found.")

  # Species names from the *first* column of each triplet; fall back to generic
  sp_firsts <- 2L + 3L * (0:(sp_n - 1L))
  sp_names  <- names(x)[sp_firsts]
  if (is.null(sp_names) || any(!nzchar(sp_names))) sp_names <- paste0("sp", seq_len(sp_n))

  # --- parse block structure in column 1 to get STRESSOR membership -----------
  lab <- as.character(x[[1L]])
  is_blank <- is.na(lab) | trimws(lab) == ""
  n <- nrow(x)

  # Identify stressor "header" rows: first non-blank after a blank run
  blank_idx <- which(is_blank)
  header_idx <- integer(0)
  header_name <- character(0)
  for (b in blank_idx) {
    j <- b + 1L
    while (j <= n && (is.na(lab[j]) || trimws(lab[j]) == "")) j <- j + 1L
    if (j <= n) {
      header_idx  <- c(header_idx, j)
      header_name <- c(header_name, lab[j])
    }
  }
  # Build a row-wise STRESSOR assignment: rows after a header up to next blank
  stressor_by_row <- rep(NA_character_, n)
  if (length(header_idx)) {
    for (k in seq_along(header_idx)) {
      start <- header_idx[k] + 1L
      nb <- blank_idx[blank_idx > header_idx[k]]
      end <- if (length(nb)) nb[1L] - 1L else n
      if (start <= end) stressor_by_row[start:end] <- header_name[k]
    }
  }
  # We'll treat any rows before the first header as "species attributes" (STRESSOR = NA)

  # Helper: does the first row look like an internal header row for columns?
  looks_like_internal_header <- function(v) {
    vv <- toupper(trimws(as.character(v)))
    any(grepl("^ATTR|^RAT|^DQ$|^WEI|^E/C$|^CRITER", vv))
  }

  # --- build per-species long tables ------------------------------------------
  out <- vector("list", sp_n)

  for (i in seq_len(sp_n)) {
    # ATTRIBUTES (col 1), this species' triplet, E/C column
    start <- 2L + 3L * (i - 1L)
    cols  <- c(1L, start:(start + 2L), crit_col)
    sub   <- x[, cols, drop = FALSE]

    # If row 1 contains an internal header line, use it for names then drop it
    r1 <- as.character(sub[1L, , drop = TRUE])
    if (looks_like_internal_header(r1)) {
      nm <- r1
      nm[nchar(nm) == 0L] <- NA_character_
      nm[is.na(nm)] <- c("ATTRIBUTES","RATING","DQ","WEIGHT","E/C")[seq_along(nm)]
      names(sub) <- nm
      data_rows <- 2L:nrow(sub)
    } else {
      names(sub) <- c("ATTRIBUTES","RATING","DQ","WEIGHT","E/C")
      data_rows <- 1L:nrow(sub)
    }

    # Drop blank rows and stressor header rows from data
    attrv <- as.character(sub$ATTRIBUTES)
    drop_rows <- which(is.na(attrv) | trimws(attrv) == "" | attrv %in% header_name)
    keep <- setdiff(data_rows, drop_rows)

    if (!length(keep)) {
      warning("No data rows found for species ", sQuote(sp_names[i]), ". Returning empty table.")
      out[[i]] <- data.frame(STRESSOR=character(0), ATTRIBUTES=character(0),
                             RATING=integer(0), DQ=integer(0), WEIGHT=integer(0),
                             `E/C`=character(0))
      next
    }

    df <- sub[keep, , drop = FALSE]
    # Attach STRESSOR label per original row index
    df$STRESSOR <- stressor_by_row[keep]

    # Coerce numeric columns
    suppressWarnings({
      df$RATING <- as.integer(df$RATING)
      df$DQ     <- as.integer(df$DQ)
      df$WEIGHT <- as.integer(df$WEIGHT)
    })

    # Reorder columns
    df <- df[, c("STRESSOR","ATTRIBUTES","RATING","DQ","WEIGHT","E/C")]
    rownames(df) <- NULL
    out[[i]] <- df
  }

  names(out) <- sp_names
  out
}


################################################################################
###############################export_maps######################################
################################################################################

#' Export nested list of maps to files
#'
#' Recursively writes `terra::SpatRaster` (.tif) and vector layers using the
#' chosen driver (ESRI Shapefile, GPKG, or GeoJSON) from a nested list to disk.
#' Optionally zips the top-level output directory.
#'
#' @param x A nested list containing `SpatRaster`, `sf`/`sfc`, or further lists.
#' @param out_dir Character path to output directory.
#' @param path Internal use: character vector of names under recursion.
#' @param zip_export Logical; if `TRUE`, zip the top-level output directory.
#' @param vector_driver One of `"ESRI Shapefile"`, `"GPKG"`, or `"GeoJSON"`.
#'   Controls the format for vector outputs. Default `"ESRI Shapefile"`.
#' @return Invisibly returns `NULL`. Files are written to `out_dir`.
#' @importFrom terra writeRaster
#' @importFrom sf st_write st_as_sf
#' @importFrom utils zip
#' @examples
#' # risa_maps <- risa_prep(spp_df, str_df)
#' # export_maps(risa_maps, "my_folder", vector_driver = "GPKG", zip_export = TRUE)
#' @export
export_maps <- function(
    x,
    out_dir,
    path = character(),
    zip_export = FALSE,
    vector_driver = c("ESRI Shapefile","GPKG","GeoJSON")
) {
  vector_driver <- match.arg(vector_driver)

  # Ensure destination exists
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Name helpers ---------------------------------------------------------------
  sanitize_name <- function(s) {
    s <- if (is.null(s) || !nzchar(s)) "layer" else as.character(s)
    s <- gsub("[^A-Za-z0-9_\\-]+", "_", s)
    s <- sub("^_+", "", s)
    s <- substr(s, 1L, 60L)
    if (!nzchar(s)) "layer" else s
  }
  base_name <- if (length(path)) sanitize_name(tail(path, 1)) else "layer"

  # Writers --------------------------------------------------------------------
  write_raster_safely <- function(r, fname_base) {
    out_file <- file.path(out_dir, paste0(fname_base, ".tif"))
    tryCatch(
      terra::writeRaster(
        r,
        filename  = out_file,
        filetype  = "GTiff",
        overwrite = TRUE
        # , gdal = c("COMPRESS=LZW")
      ),
      error = function(e) warning("Failed to write raster: ", out_file, " (", e$message, ")", call. = FALSE)
    )
  }

  # pick dsn + layer based on vector_driver
  vector_target <- function(folder, fname_base, driver) {
    if (driver == "ESRI Shapefile") {
      dsn   <- folder
      layer <- fname_base
    } else if (driver == "GPKG") {
      dsn   <- file.path(folder, paste0(fname_base, ".gpkg"))
      layer <- fname_base
    } else { # GeoJSON
      dsn   <- file.path(folder, paste0(fname_base, ".geojson"))
      layer <- fname_base
    }
    list(dsn = dsn, layer = layer)
  }

  write_vector_safely <- function(v, layer_base) {
    if (inherits(v, "sfc")) v <- sf::st_as_sf(v)
    tgt <- vector_target(out_dir, layer_base, vector_driver)
    drv <- vector_driver
    tryCatch(
      sf::st_write(
        obj          = v,
        dsn          = tgt$dsn,
        layer        = tgt$layer,
        driver       = drv,
        delete_layer = TRUE,
        quiet        = TRUE
      ),
      error = function(e) {
        dest <- if (drv == "ESRI Shapefile") file.path(out_dir, paste0(layer_base, ".shp")) else tgt$dsn
        warning("Failed to write vector: ", dest, " (", e$message, ")", call. = FALSE)
      }
    )
  }

  # Dispatch -------------------------------------------------------------------
  if (inherits(x, "SpatRaster")) {
    write_raster_safely(x, base_name)

  } else if (inherits(x, "sf") || inherits(x, "sfc")) {
    write_vector_safely(x, base_name)

  } else if (is.list(x)) {
    nms <- names(x)
    for (i in seq_along(x)) {
      nm <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else sprintf("item_%02d", i)
      nm <- sanitize_name(nm)
      sub_dir <- file.path(out_dir, nm)
      export_maps(
        x            = x[[i]],
        out_dir      = sub_dir,
        path         = c(path, nm),
        zip_export   = FALSE,
        vector_driver = vector_driver
      )
    }

  } else {
    message("Skipping unrecognized object at: ", paste(path, collapse = "$"))
  }

  # Zip only once at the top level --------------------------------------------
  if (isTRUE(zip_export) && length(path) == 0L) {
    zip_file <- paste0(normalizePath(out_dir, mustWork = FALSE), ".zip")
    message("Zipping export folder to: ", zip_file)
    utils::zip(
      zipfile = zip_file,
      files   = list.files(out_dir, full.names = TRUE, recursive = TRUE)
    )
  }

  invisible(NULL)
}


################################################################################
###############################input_tools######################################
################################################################################

#' Read spatial files into a nested list mirroring directory structure
#'
#' Scans a root directory (non-recursively at each level) for raster files
#' (`.tif`, `.tiff`) and vector files (`.shp`; also supports `.gpkg`, `.geojson`)
#' and reads them as `terra::SpatRaster` or `sf` objects. Subdirectories are
#' traversed recursively; the returned list mirrors the folder structure.
#'
#' @param dir_path character(1). Path to the root directory containing spatial data.
#' @return A nested `list`:
#'   - Intermediate elements are named by subdirectory names and contain further lists.
#'   - Leaf elements are named by file basenames (without extensions) and contain
#'     `SpatRaster` or `sf` objects.
#' @importFrom terra rast
#' @importFrom sf st_read
#' @importFrom tools file_ext file_path_sans_ext
#' @export
read_maps_nested <- function(dir_path) {
  if (!is.character(dir_path) || length(dir_path) != 1L) {
    stop("`dir_path` must be a single character path.")
  }
  if (!dir.exists(dir_path)) {
    stop("Directory not found: ", dir_path)
  }

  # internal recursive worker
  .read_dir <- function(d) {
    out <- list()

    # files in this directory (not recursive)
    files <- list.files(
      d,
      pattern = "(?i)\\.(tif|tiff|shp|gpkg|geojson)$",
      full.names = TRUE,
      recursive = FALSE
    )

    for (f in files) {
      ext  <- tolower(tools::file_ext(f))
      name <- tools::file_path_sans_ext(basename(f))

      # avoid name collisions within this folder
      if (name %in% names(out)) {
        name <- make.unique(c(names(out), name))[length(names(out)) + 1L]
      }

      obj <- switch(
        ext,
        tif     = tryCatch(terra::rast(f), error = function(e) { warning("Failed to read raster: ", f, " (", e$message, ")"); NULL }),
        tiff    = tryCatch(terra::rast(f), error = function(e) { warning("Failed to read raster: ", f, " (", e$message, ")"); NULL }),
        shp     = tryCatch(sf::st_read(f, quiet = TRUE), error = function(e) { warning("Failed to read vector: ", f, " (", e$message, ")"); NULL }),
        gpkg    = tryCatch(sf::st_read(f, quiet = TRUE), error = function(e) { warning("Failed to read vector: ", f, " (", e$message, ")"); NULL }),
        geojson = tryCatch(sf::st_read(f, quiet = TRUE), error = function(e) { warning("Failed to read vector: ", f, " (", e$message, ")"); NULL }),
        { warning("Unsupported extension: ", ext, " for file ", f); NULL }
      )

      if (!is.null(obj)) out[[name]] <- obj
    }

    # subdirectories (not recursive here; recurse manually)
    subdirs <- list.dirs(d, full.names = TRUE, recursive = FALSE)
    # guard against weird entries
    if (length(subdirs)) {
      is_dir <- dir.exists(subdirs)
      subdirs <- subdirs[is_dir]
    }

    for (sd in subdirs) {
      nm <- basename(sd)
      # avoid name collisions with files already added
      if (nm %in% names(out)) {
        nm <- make.unique(c(names(out), nm))[length(names(out)) + 1L]
      }
      out[[nm]] <- .read_dir(sd)  # <- correct recursive call
    }

    out
  }

  .read_dir(dir_path)
}







