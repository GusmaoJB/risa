#' Compute overlap hotspots from two reclassified or continuous rasters
#'
#' Combines two class-coded or continuous rasters (e.g., exposure & consequence) into an
#' overlap/risk surface. For discrete inputs (default), values are reclassified to
#' `out_classes` bins. For continuous inputs, the result is a continuous surface rescaled
#' to a specified range. Inputs are assumed to have values in `1…n_classes` for discrete
#' mode, or any positive continuous values for continuous mode (0 or negative values treated
#' as background/NA).
#'
#' @param x,y Single-layer `terra::SpatRaster`. For discrete mode, values are class codes
#'   in `1…n_classes`. For continuous mode, values can be any positive floats.
#' @param n_classes Integer: number of classes in the inputs (default 3). For continuous
#'   mode, this sets the normalization range for input values.
#' @param continuous Logical: if `TRUE`, treats inputs as continuous rasters and returns a
#'   continuous output raster. If `FALSE` (default), treats inputs as discrete class rasters.
#' @param output_min Numeric or NULL: minimum value for continuous output scaling. If NULL
#'   (default), uses 1. Only used when `continuous = TRUE`.
#' @param out_classes Integer: number of classes in the discrete output (default = `n_classes`).
#'   When `continuous = TRUE`, this sets the maximum value for output scaling.
#' @param method Combination rule: `"product"` (default), `"sum"`, `"geom_mean"`, or `"max"`.
#' @param resample_method Resampling method to align `y` to `x` when grids differ. For class
#'   rasters, the default `"near"` preserves class labels.
#' @param output_layer_type One of `"shp"`, `"raster"`, or `"both"`. Default `"shp"`.
#'   Ignored when `continuous = TRUE` (always returns raster).
#' @param quiet Logical; suppress informative messages. Default `TRUE`.
#'
#' @details
#' The function handles CRS mismatches automatically by reprojecting `y` to match `x`.
#' If both rasters lack a CRS, the function assumes they share the same projection.
#' An error occurs if only one raster has a defined CRS.
#'
#' The available combination rules are:
#' \describe{
#'   \item{`"product"`}{Multiplies the two rasters. For discrete inputs, scores range from 1 to `n_classes^2`.}
#'   \item{`"sum"`}{Adds the two rasters. For discrete inputs, scores range from 2 to `2 * n_classes`.}
#'   \item{`"geom_mean"`}{Computes the geometric mean, `sqrt(x * y)`.}
#'   \item{`"max"`}{Uses the maximum value between `x` and `y` for each cell.}
#' }
#'
#' For **discrete mode** (`continuous = FALSE`):
#' * Inputs must be integer class codes in the range `1…n_classes`
#' * Output is classified into `out_classes` discrete bins
#' * Can return a `SpatRaster` raster, an `sf` polygon, or both
#'
#' For **continuous mode** (`continuous = TRUE`):
#' * Inputs can be any positive continuous values
#' * Inputs are normalized to `[1, n_classes]` before combination
#' * Output is a continuous raster scaled to `[output_min (or 1), out_classes]`
#' * Always returns a raster (ignores `output_layer_type`)
#' * In continuous mode, if an input raster has a constant finite value (`min == max`),
#' all finite cells are normalized to 1 to avoid division by zero. This treats rasters
#' without internal variation as having the minimum normalized rating.
#'
#' @return
#' * If `continuous = TRUE`: A `SpatRaster` with continuous values.
#' * If `continuous = FALSE` and `output_layer_type = "raster"`: A `SpatRaster` with integer classes.
#' * If `continuous = FALSE` and `output_layer_type = "shp"`: An `sf` polygon layer.
#' * If `continuous = FALSE` and `output_layer_type = "both"`: A list with `raster` and `shp` elements.
#'
#' @importFrom terra nlyr crs same.crs project compareGeom resample ifel global mosaic clamp classify as.polygons
#' @importFrom sf st_as_sf st_set_crs
#'
#' @examples
#' # Discrete mode example
#' species  <- data.frame(long = rnorm(80, 0, 10),  lat = rnorm(80, 0, 10))
#' stressor <- data.frame(long = rnorm(100, 0, 10), lat = rnorm(100, 0, 10))
#' kde_spe <- get_class_kernel(species,  output_layer_type = "raster")
#' kde_str <- get_class_kernel(stressor, output_layer_type = "raster")
#' # Ensure same grid (nearest to preserve class labels)
#' kde_str <- terra::project(kde_str, kde_spe, method = "near")
#' overlap <- get_overlap_kernel(kde_spe, kde_str, method = "product",
#'                               output_layer_type = "raster")
#' terra::plot(overlap)
#'
#' \dontrun{
#' # Continuous mode example
#' cont_rast1 <- terra::rast(matrix(runif(100, 0, 100), 10, 10))
#' cont_rast2 <- terra::rast(matrix(runif(100, 0, 50), 10, 10))
#' cont_overlap <- get_overlap_kernel(cont_rast1, cont_rast2,
#'                                    continuous = TRUE,
#'                                    n_classes = 5,
#'                                    out_classes = 10,
#'                                    output_min = 2)
#' terra::plot(cont_overlap)
#' }
#' @export
get_overlap_kernel <- function(
    x, y,
    n_classes = 3,
    continuous = FALSE,
    output_min = NULL,
    out_classes = n_classes,
    method = c("product","sum","geom_mean","max"),
    resample_method = "near",
    output_layer_type = c("shp","raster","both"),
    quiet = TRUE) {

  method  <- match.arg(method)
  output_layer_type <- match.arg(output_layer_type)

  # Basic checks
  if (!inherits(x, "SpatRaster") || !inherits(y, "SpatRaster")) {
    stop("`x` and `y` must be terra::SpatRaster.")
  }
  if (terra::nlyr(x) != 1L || terra::nlyr(y) != 1L) {
    stop("`x` and `y` must be single-layer rasters.")
  }

  # CRS handling
  crs_x <- terra::crs(x)
  crs_y <- terra::crs(y)
  crs_x_known <- crs_x != ""
  crs_y_known <- crs_y != ""

  if (!crs_x_known && !crs_y_known) {
    if (!quiet) message("Both input rasters do not have a known CRS. Assuming they have the same projection and CRS...")
  } else if (!crs_x_known || !crs_y_known) {
    stop("One of the input rasters does not have a known CRS. Both must have known CRS, or both must have unknown CRS.")
  } else if (!terra::same.crs(x, y)) {
    if (!quiet) message("Reprojecting `y` to match CRS of `x`.")
    y <- terra::project(y, x)
  }

  # Helper for normalize rasters
  .normalize_to_classes <- function(r, min_val, max_val, n_classes, raster_name) {

    if (!is.finite(min_val) || !is.finite(max_val)) {
      stop("`", raster_name, "` has no finite values to normalize.")
    }

    if (max_val == min_val) {

      if (!quiet) {
        message(
          "`", raster_name, "` has a constant value. ",
          "Normalizing all finite cells to 1."
        )
      }

      return(terra::ifel(is.finite(r), 1, NA))
    }

    (r - min_val) / (max_val - min_val) * (n_classes - 1) + 1
  }

  # Align y to x grid if needed
  if (!terra::compareGeom(x, y, stopOnError = FALSE)) {
    if (!quiet) message("Aligning `y` to `x` with terra::resample(..., method = '", resample_method, "').")
    y <- terra::resample(y, x, method = resample_method)
  }

  # Treat 0 and negative values
  if (continuous) {
    # For continuous data, keep all finite values, but set negatives to NA
    x <- terra::ifel(x <= 0, NA, x)
    y <- terra::ifel(y <= 0, NA, y)

    # Use actual min/max for normalization
    mmx <- terra::global(x, c("min","max"), na.rm = TRUE)
    mmy <- terra::global(y, c("min","max"), na.rm = TRUE)
    minx <- as.numeric(mmx[1,"min"]); maxx <- as.numeric(mmx[1,"max"])
    miny <- as.numeric(mmy[1,"min"]); maxy <- as.numeric(mmy[1,"max"])

    # Normalize inputs to [1, n_classes] if they're not already in that range
    if (!quiet) message("Normalizing continuous inputs to [1, ", n_classes, "] range.")
    x <- .normalize_to_classes(
      r = x,
      min_val = minx,
      max_val = maxx,
      n_classes = n_classes,
      raster_name = "x")

    y <- .normalize_to_classes(
      r = y,
      min_val = miny,
      max_val = maxy,
      n_classes = n_classes,
      raster_name = "y")

  } else {
    # For integer classes
    x <- terra::ifel(x <= 0, NA, x)
    y <- terra::ifel(y <= 0, NA, y)

    # Validate class ranges
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
  }

  # Combine
  if (method == "product") {
    score <- x * y
    s_min <- 1
    s_max <- n_classes^2
  } else if (method == "sum") {
    score <- x + y
    s_min <- 2
    s_max <- 2 * n_classes
  } else if (method == "geom_mean") {
    score <- sqrt(x * y)
    s_min <- 1
    s_max <- n_classes
  } else if (method == "max") {
    score <- terra::mosaic(x, y, fun = function(a, b) pmax(a, b, na.rm = TRUE))
    s_min <- 1
    s_max <- n_classes
  }

  # Handle output based on continuous flag
  if (continuous) {
    if (!quiet) message("Rescaling continuous output to [1, ", out_classes, "] range.")
    if (is.null(output_min)) {
      out_min <- 1
    } else {
      out_min <- output_min
    }
    out_max <- out_classes

    result <- (score - s_min) / (s_max - s_min) * (out_max - out_min) + out_min

    names(result) <- "Rating"
    terra::crs(result) <- terra::crs(x)

    return(result)

  } else {
    norm <- (score - s_min) / (s_max - s_min)
    norm <- terra::clamp(norm, lower = 0, upper = 1)

    brks <- seq(0, 1, length.out = out_classes + 1)
    mat <- cbind(from = brks[-length(brks)], to = brks[-1], class = seq_len(out_classes))
    result <- terra::classify(norm, rcl = mat, include.lowest = TRUE)
  }

  names(result) <- "Rating"
  terra::crs(result) <- terra::crs(x)

  # Vectorize if requested (only for discrete output)
  if (output_layer_type != "raster") {
    v <- sf::st_as_sf(terra::as.polygons(result, dissolve = TRUE))
    v <- sf::st_set_crs(v, terra::crs(result))
  }

  # Return
  if (output_layer_type == "shp") {
    return(v)
  } else if (output_layer_type == "raster") {
    return(result)
  } else {
    return(list(raster = result, shp = v))
  }
}
