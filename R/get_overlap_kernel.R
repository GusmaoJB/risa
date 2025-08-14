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
#' @importFrom terra same.crs compareGeom resample classify as.polygons global rast values mask crs
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

  # Basic checks
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

  # Combine
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
    score    <- terra::mosaic(x, y, fun = function(a, b) pmax(a, b, na.rm = TRUE))
    s_min    <- 1
    s_max    <- n_classes
  }

  # Normalize to [0,1] using theoretical bounds, then reclass to out_classes bins
  norm <- (score - s_min) / (s_max - s_min)

  brks <- seq(0, 1, length.out = out_classes + 1)
  mat  <- cbind(from = brks[-length(brks)], to = brks[-1], class = seq_len(out_classes))
  result_reclass <- terra::classify(norm, rcl = mat, include.lowest = TRUE)
  names(result_reclass) <- "Rating"
  terra::crs(result_reclass) <- terra::crs(x)

  # Vectorize if requested
  if (output_layer_type != "raster") {
    v <- sf::st_as_sf(terra::as.polygons(result_reclass, dissolve = TRUE))
    v <- sf::st_set_crs(v, terra::crs(result_reclass))
  }

  # Return
  if (output_layer_type == "shp") {
    return(v)
  } else if (output_layer_type == "raster") {
    return(result_reclass)
  } else {
    return(list(raster = result_reclass, shp = v))
  }
}
