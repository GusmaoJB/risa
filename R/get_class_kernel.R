#' Kernel density estimation and reclassification wrapper
#'
#' Computes KDE for point data, reclassifies into classes,
#' and returns a raster (`SpatRaster`), a vector (`sf`), or both.
#' The bandwidth of the KDE smoothing is based on the radius parameter.
#' If not informed, the radius is set as 10 percent of mean pairwise distance among points.
#'
#' @param x An `sf` object or data frame of points (lon/lat) or tibble.
#' @param area Bounding region (`bbox`, `sf`, or data frame) or `NULL` to auto-compute.
#' @param n_classes Integer number of reclassification bins.
#' @param output_layer_type Character; one of `"shp"`, `"raster"`, or `"both"`.
#' @param radius Numeric smoothing bandwidth; if `NULL`, set to 10 percent of mean pairwise distance.
#' @param group_size Optional column name for weights.
#' @returns An `sf` object, `SpatRaster`, or list containing both.
#' @importFrom sf st_as_sf st_crs st_as_sfc
#' @importFrom spatstat.geom as.owin ppp
#' @importFrom spatstat.explore density.ppp
#' @importFrom terra rast crs disagg classify as.polygons
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
get_class_kernel <- function(x, area = NULL, n_classes = 3, output_layer_type = "shp", radius = NULL, group_size = NULL) {

  # Convert input data to sf and metric CRS
  if (inherits(x, "sf")) {
    message("Input is a vector ('sf' object).")
    x_sf <- x
  } else if (inherits(x, c("data.frame", "tbl_df", "tbl"))) {
    if (ncol(x) < 2) {
      stop("Data frame must have at least two columns: lon and lat.")
    }
    if (!all(sapply(x[, 1:2], is.numeric))) {
      stop("First two columns in the data frame must be numeric.")
    }
    message("Input is a data.frame. Converting to vector (sf)...")
    x_sf <- df_to_shp(x)
  } else {
    stop("Input must be either an 'sf' object or a numeric data.frame with at least two columns.")
  }

  # Prepare area
  if (is.null(area)) {
    message("No area provided. Using a buffered bounding box around observations as working space...")
    area <- create_area(x_sf, area_type="bbox")
  } else {
    if (inherits(area, "bbox")) {
      message("Area provided is a bounding box. Converting to sf...")
      area <- st_as_sfc(area)
    } else if (inherits(area, "sf")) {
      message("Area provided is a sf object")
    } else if (inherits(area, "data.frame")) {
      message("Area provided is a dataframe of coordinates. Converting to sf...")
      area <- df_to_shp(area)
    } else {
      stop("Provided area must be bbox, sf or data.frame object.")
    }
  }

  area_m <- transform_to_metric(area)
  crs_proj <- st_crs(area$shape)$wkt
  x_m <- transform_to_metric(x_sf, metric_crs = st_crs(area_m$shape)$epsg)
  coords_x <- x_m$coordinates

  # Convert area to spatstat window
  window <- as.owin(area_m$shape)

  # Check for group_size
  if (!is.null(group_size)) {
    if (!group_size %in% colnames(x)) {
      stop("Specified 'group_size' column does not exist in the input data.")
    }
    weights <- x[[group_size]]
    if (!is.numeric(weights) || any(weights <= 0)) {
      stop("'group_size' column must be numeric and positive.")
    }
    ppp_x <- ppp(x = coords_x[, 1], y = coords_x[, 2], window = window, marks = weights)
  } else {
    ppp_x <- ppp(x = coords_x[, 1], y = coords_x[, 2], window = window)
  }

  # Estimate radius if needed
  if (is.null(radius)) {
    radius <- mean(dist(coords_x)) * 0.1
    message(sprintf("Auto-selected radius: %.2f units.", radius))
  }

  # Compute KDE with or without weights
  if (!is.null(group_size)) {
    kde <- spatstat.explore::density.ppp(ppp_x, weights = weights, sigma = radius, dimyx = c(800, 800))
  } else {
    kde <- spatstat.explore::density.ppp(ppp_x, sigma = radius, dimyx = c(800, 800))
  }

  r_kde <- terra::disagg(rast(kde), fact = 3, method = 'bilinear')
  terra::crs(r_kde) <- crs_proj

  # Reclassify raster
  re_mat <- reclass_matrix(r_kde, n_classes = n_classes)
  r_kde_reclass <- terra::classify(r_kde, re_mat, include.lowest = TRUE)
  terra::crs(r_kde_reclass) <- crs_proj
  names(r_kde_reclass) <- "Rating"

  # Convert to polygons
  vec_kde_reclass <- st_as_sf(as.polygons(r_kde_reclass, dissolve = TRUE), crs = crs_proj)

  # Output
  if (output_layer_type == "shp") {
    return(vec_kde_reclass)
  } else if (output_layer_type == "raster") {
    return(r_kde_reclass)
  } else if (output_layer_type == "both") {
    return(list(raster = r_kde_reclass, shp = vec_kde_reclass))
  }
}
