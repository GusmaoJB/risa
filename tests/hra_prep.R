#' Prepare rasters for HRA: project to AOI CRS and resample to a square metric grid
#'
#' @param criteria_rasters A single SpatRaster, a (named) list of SpatRasters,
#'   or a nested list (depth 2 or 3) as in the package's HRA structure.
#' @param distribution_rasters A single SpatRaster or a (named) list of SpatRasters,
#'   typically one per species for ecosystem runs.
#' @param area_of_interest An 'sf' object (polygon/collection) or a 'SpatRaster'.
#'   Its CRS/projection will define the target CRS. If geographic, it will be transformed
#'   to a metric CRS using \code{transform_to_metric()}.
#' @param resolution Numeric (meters). Target square cell size for all outputs.
#'
#' @return A list with two elements: \code{criteria_rasters} and \code{distribution_rasters},
#'   each mirroring the input container shape but with rasters (1) in a metric CRS,
#'   (2) aligned to a common, north-up, square grid defined by the AOI extent and \code{resolution},
#'   and (3) masked to the AOI.
#'
#' @details
#' - Uses helper functions available in the package:
#'   \code{list_depth_base()}, \code{guess_crs()}, \code{transform_to_metric()}.
#' - Resampling method defaults to "near" to avoid smoothing/categorical corruption;
#'   downstream code can reproject with a continuous method if/when needed.
#' - Ensures square pixels and meter units (prevents GDAL "pixels not square" distance warnings).
#'
#' @import terra sf
#' @export
hra_prep <- function(criteria_rasters,
                     distribution_rasters,
                     area_of_interest,
                     resolution) {

  stopifnot(is.numeric(resolution), length(resolution) == 1L, resolution > 0)

  # --- Helpers ---------------------------------------------------------------

  .as_raster <- function(x) {
    if (inherits(x, "SpatRaster")) return(x)
    if (inherits(x, "Raster"))     return(terra::rast(x))
    NULL
  }

  .is_rasterish <- function(x) inherits(x, "SpatRaster") || inherits(x, "Raster")

  .map_rasters <- function(x, f) {
    # Recursively apply f() to any SpatRaster; preserve list structure
    if (.is_rasterish(x)) return(f(.as_raster(x)))
    if (is.list(x)) {
      out <- x
      for (nm in names(x)) out[[nm]] <- .map_rasters(x[[nm]], f)
      if (is.null(names(x))) out <- lapply(x, .map_rasters, f = f)
      return(out)
    }
    x
  }

  .any_raster <- function(x) {
    if (.is_rasterish(x)) return(TRUE)
    if (is.list(x)) return(any(vapply(x, .any_raster, logical(1))))
    FALSE
  }

  .pick_any_raster <- function(x) {
    if (.is_rasterish(x)) return(.as_raster(x))
    if (is.list(x)) {
      for (el in x) {
        r <- .pick_any_raster(el)
        if (!is.null(r)) return(r)
      }
    }
    NULL
  }

  .normalize_aoi <- function(aoi) {
    if (inherits(aoi, "SpatRaster")) {
      return(aoi)
    } else if (inherits(aoi, "sf")) {
      if (is.na(sf::st_crs(aoi))) {
        # Try to guess if missing
        aoi <- sf::st_set_crs(aoi, guess_crs(aoi))
      }
      # Convert to SpatRaster-like extent carrier (we'll rasterize the grid later)
      return(aoi)
    } else {
      stop("'area_of_interest' must be an 'sf' object or a 'SpatRaster'.")
    }
  }

  .target_crs <- function(aoi) {
    if (inherits(aoi, "SpatRaster")) {
      cr <- terra::crs(aoi)
      if (is.na(cr) || terra::crs(aoi, describe = TRUE)$is_lonlat) {
        # Transform AOI grid to metric CRS
        aoi_m <- transform_to_metric(aoi)  # package helper
        return(terra::crs(aoi_m))
      } else {
        return(cr)
      }
    } else if (inherits(aoi, "sf")) {
      if (is.na(sf::st_crs(aoi))) {
        aoi <- sf::st_set_crs(aoi, guess_crs(aoi))
      }
      if (sf::st_is_longlat(aoi)) {
        aoi_m <- transform_to_metric(aoi)  # package helper
        return(sf::st_crs(aoi_m)$wkt)
      } else {
        return(sf::st_crs(aoi)$wkt)
      }
    }
    stop("Could not determine target CRS from 'area_of_interest'.")
  }

  .aoi_extent_vect <- function(aoi, crs_wkt) {
    # Return extent in target CRS and an AOI geometry for masking
    if (inherits(aoi, "SpatRaster")) {
      r <- if (!identical(terra::crs(aoi), crs_wkt)) terra::project(aoi, crs_wkt, method = "near") else aoi
      return(list(ext = terra::ext(r), aoi_geom = r))
    } else if (inherits(aoi, "sf")) {
      g <- if (!identical(sf::st_crs(aoi)$wkt, crs_wkt)) sf::st_transform(aoi, crs = crs_wkt) else aoi
      bb <- sf::st_bbox(g)
      ext <- terra::ext(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"])
      return(list(ext = ext, aoi_geom = g))
    }
    stop("Unsupported AOI type.")
  }

  .mk_template <- function(ext, crs_wkt, res_m) {
    # Square, north-up, meter-based grid
    terra::rast(extent = ext, crs = crs_wkt, resolution = c(res_m, res_m))
  }

  .project_resample_mask <- function(r, template, aoi_geom) {
    # Always use "near" here to avoid smoothing / class pollution.
    # Downstream functions can choose bilinear where appropriate.
    r2 <- if (!identical(terra::crs(r), terra::crs(template))) {
      terra::project(r, template, method = "near")
    } else r
    r3 <- terra::resample(r2, template, method = "near")
    # Mask to AOI
    if (inherits(aoi_geom, "sf")) {
      r4 <- terra::mask(r3, vect(aoi_geom))
    } else {
      # SpatRaster AOI: mask non-NA cells
      r4 <- terra::mask(r3, aoi_geom)
    }
    r4
  }

  # --- Validate presence of rasters -----------------------------------------

  if (!.any_raster(criteria_rasters)) {
    stop("'criteria_rasters' must contain at least one SpatRaster (possibly nested).")
  }
  if (!.any_raster(distribution_rasters)) {
    stop("'distribution_rasters' must contain at least one SpatRaster (possibly nested).")
  }

  dcrit <- list_depth_base(criteria_rasters)
  ddist <- list_depth_base(distribution_rasters)
  if (!(dcrit %in% c(1L, 2L, 3L))) {
    stop("'criteria_rasters' must be a raster, or a list of rasters (depth 2 or 3).")
  }
  if (!(ddist %in% c(1L, 2L))) {
    stop("'distribution_rasters' must be a raster or a list (depth 2).")
  }

  # --- Determine target CRS & template from AOI ------------------------------

  aoi <- .normalize_aoi(area_of_interest)
  tgt_crs <- .target_crs(aoi)
  aoi_info <- .aoi_extent_vect(aoi, tgt_crs)
  template <- .mk_template(aoi_info$ext, tgt_crs, res_m = resolution)

  # Safety: square & meter grid (avoid GDAL distance warnings later)
  rs <- terra::res(template)
  if (!isTRUE(all.equal(rs[1], rs[2]))) {
    stop("Internal: template resolution not square; please report this bug.")
  }
  if (terra::crs(template, describe = TRUE)$is_lonlat) {
    stop("Internal: template is lon/lat; expected metric CRS. Please report this bug.")
  }

  # --- Reproject, resample, mask everywhere ---------------------------------

  project_one <- function(r) .project_resample_mask(r, template, aoi_info$aoi_geom)

  criteria_out     <- .map_rasters(criteria_rasters, project_one)
  distribution_out <- .map_rasters(distribution_rasters, project_one)

  # Return shape-preserving containers
  list(
    criteria_rasters      = criteria_out,
    distribution_rasters  = distribution_out
  )
}
