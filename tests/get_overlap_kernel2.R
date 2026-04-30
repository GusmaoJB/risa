library(risa)

get_overlap_kernel2 <- function(
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
    message("Both input rasters do not have a known CRS. I am assuming they have the same projection and CRS.")
  } else if (!crs_x_known || !crs_y_known) {
    stop("One of the input rasters does not have a known CRS. Both must have known CRS, or both must have unknown CRS.")
  } else if (!terra::same.crs(x, y)) {
    if (!quiet) message("Reprojecting `y` to match CRS of `x`.")
    y <- terra::project(y, x)
  }

  # Align y to x grid if needed
  if (!terra::compareGeom(x, y, stopOnError = FALSE)) {
    if (!quiet) message("Aligning `y` to `x` with terra::resample(..., method = '", resample_method, "').")
    y <- terra::resample(y, x, method = resample_method)
  }

  # Treat 0 and negative values
  if (continuous) {
    # For continuous data, keep all finite values, but set negatives to NA
    x <- terra::ifel(x < 0, NA, x)
    y <- terra::ifel(y < 0, NA, y)

    # Use actual min/max for normalization
    mmx <- terra::global(x, c("min","max"), na.rm = TRUE)
    mmy <- terra::global(y, c("min","max"), na.rm = TRUE)
    minx <- as.numeric(mmx[1,"min"]); maxx <- as.numeric(mmx[1,"max"])
    miny <- as.numeric(mmy[1,"min"]); maxy <- as.numeric(mmy[1,"max"])

    # Normalize inputs to [1, n_classes] if they're not already in that range
    if (!quiet) message("Normalizing continuous inputs to [1, ", n_classes, "] range.")
    x <- (x - minx) / (maxx - minx) * (n_classes - 1) + 1
    y <- (y - miny) / (maxy - miny) * (n_classes - 1) + 1

  } else {
    # Original behavior for integer classes
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
