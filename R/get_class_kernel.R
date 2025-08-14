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
#' @examples
#' # Creating test data
#' df <- data.frame(long = rnorm(120, 0, 10), lat = rnorm(120, 0, 10))
#'
#' # Generating reclassified Kernel densities estimates (3 classes)
#' kde <- get_class_kernel(df)
#'
#' # Plot KDE map
#' plot(kde)
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
