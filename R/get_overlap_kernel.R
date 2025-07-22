#' Compute overlap hotspots from two reclassified rasters
#'
#' Multiplies two classed rasters and reclassifies the result into risk ratings.
#'
#' @param x A reclassified `SpatRaster` of class codes.
#' @param y A reclassified `SpatRaster` with same CRS and class range.
#' @param n_classes Integer number of classes used in inputs.
#' @param output_layer_type Character; one of `"shp"`, `"raster"`, or `"both"`.
#' @return An `sf` object, `SpatRaster`, or list containing both.
#' @importFrom terra same.crs compareGeom resample classify
#' @importFrom sf st_as_sf
#' @examples
#' # Creating test data
#' species <- data.frame(long = rnorm(80, 0, 10), lat = rnorm(80, 0, 10))
#' stressor <- data.frame(long = rnorm(100, 0, 10), lat = rnorm(100, 0, 10))
#' kde_spe <- get_class_kernel(species,  output_layer_type = "raster")
#' kde_str <- get_class_kernel(stressor,  output_layer_type = "raster")
#' kde_spe <- terra::project(kde_spe, crs(kde_str), method = "bilinear")
#'
#' # Creating reclassified overlap raster
#' overlap <- get_overlap_kernel(kde_spe, kde_str, output_layer_type = "raster")
#'
#' # Plot KDE map
#' par(mfrow=c(1,3))
#' plot(kde_spe, main="Species distribution")
#' plot(kde_str, main="Stressor occurrence")
#' plot(overlap, main="Overlap")
#' @export
get_overlap_kernel <- function(x, y, n_classes = 3, output_layer_type = "shp") {
  # Check CRS compatibility
  if (!terra::same.crs(x, y)) {
    stop("Input rasters must have the same CRS.")
  }

  # Align y to x if needed (extent, resolution)
  if (!compareGeom(x, y, stopOnError = FALSE)) {
    message("Aligning raster y to x using terra::resample()")
    y <- resample(y, x, method = "bilinear")
  }

  # Check value ranges
  check_range <- function(r, n) {
    r_vals <- values(r)
    if (min(r_vals, na.rm = TRUE) < 0 || max(r_vals, na.rm = TRUE) > n) {
      stop(paste("Input raster values must be within the range 0 to", n,
                 "- compatible with number of classes (n_classes)."))
    }
  }

  check_range(x, n_classes)
  check_range(y, n_classes)

  # Multiplication of rasters
  result <- x * y

  # Reclassify summed result
  r_mat <- reclass_matrix(result, n_classes = n_classes, exclude_lowest=FALSE)
  result_reclass <- classify(result, r_mat, include.lowest = TRUE)
  names(result_reclass) <- "Rating"

  # Output
  if (output_layer_type == "shp") {
    vec_result_reclass <- st_as_sf(as.polygons(result_reclass, dissolve = TRUE))
    return(vec_result_reclass)
  } else if (output_layer_type == "raster") {
    return(result_reclass)
  } else if (output_layer_type == "both") {
    vec_result_reclass <- st_as_sf(as.polygons(result_reclass, dissolve = TRUE))
    return(list(raster=result_reclass, shp=vec_result_reclass))
  }
}
