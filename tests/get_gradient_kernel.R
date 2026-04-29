get_gradient_kernel <- function(
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
    return_crs = c("metric", "4326"),
    quiet = TRUE
) {
  output_layer_type <- match.arg(output_layer_type)
  radius_method <- match.arg(radius_method)
  return_crs <- match.arg(return_crs)

  # Normalize input to sf
  if (inherits(x, "sf")) {
    x_sf <- x
  } else if (inherits(x, c("data.frame", "tbl_df", "tbl"))) {
    if (ncol(x) < 2 || !all(vapply(x[, 1:2], is.numeric, logical(1)))) {
      stop("For data.frame input, the first two columns must be numeric lon/lat.")
    }
    if (!quiet) message("Converting input data.frame to sf...")
    x_sf <- df_to_shp(x)
  } else {
    stop("`x` must be an sf object or a data.frame with lon/lat.")
  }

  # Area handling
  crop_to_area <- !is.null(area)

  if (is.null(area)) {
    if (!quiet) message("No area provided. Using buffered bounding box around observations as KDE window...")
    area <- create_area(x_sf, area_type = "bbox", buffer_frac = 0.5, quiet = quiet)

  } else if (inherits(area, "bbox")) {
    area <- sf::st_as_sfc(area)

  } else if (inherits(area, "data.frame")) {
    area <- df_to_shp(area)

  } else if (!inherits(area, "sf")) {
    stop("`area` must be NULL, bbox, data.frame, or sf.")
  }

  # Project both to the same metric CRS
  area_m <- transform_to_metric(area, quiet = quiet)
  crs_m <- sf::st_crs(area_m$shape)$epsg
  x_m <- transform_to_metric(x_sf, metric_crs = crs_m, quiet = quiet)

  coords <- x_m$coordinates[, 1:2, drop = FALSE]

  if (nrow(coords) == 0L) {
    stop("No points available after preprocessing.")
  }

  if (nrow(unique(coords)) == 1L && is.null(radius)) {
    warning("Only one unique point; KDE will be a single peak. Consider setting `radius` explicitly.")
  }

  # Weights
  weights <- NULL

  if (!is.null(group_size)) {
    if (!group_size %in% names(x_m$shape)) {
      stop("`group_size` not found in `x`.")
    }

    weights <- x_m$shape[[group_size]]

    if (!is.numeric(weights)) {
      stop("`group_size` must be numeric.")
    }

    if (any(weights < 0, na.rm = TRUE)) {
      stop("`group_size` weights must be non-negative.")
    }
  }

  # Window and bandwidth
  W <- spatstat.geom::as.owin(area_m$shape)

  X <- spatstat.geom::ppp(
    x = coords[, 1],
    y = coords[, 2],
    window = W,
    check = FALSE
  )

  if (is.null(radius)) {
    if (radius_method == "nndist") {
      nn <- spatstat.geom::nndist(X)
      radius <- 1.5 * stats::median(nn[is.finite(nn)], na.rm = TRUE)

      if (!quiet) {
        message(sprintf("Auto bandwidth (nndist): sigma = %.3f", radius))
      }

    } else if (radius_method == "ppl") {
      bw <- spatstat.explore::bw.ppl(X)
      radius <- as.numeric(bw)

      if (!quiet) {
        message(sprintf("Auto bandwidth (ppl): sigma = %.3f", radius))
      }

    } else {
      stop("`radius` is NULL but `radius_method = 'fixed'`. Supply a numeric radius.")
    }
  }

  if (!is.numeric(radius) || length(radius) != 1L || !(radius > 0)) {
    stop("`radius` must be a single positive number, in projected units.")
  }

  # Grid resolution
  if (!is.null(pixel_size)) {
    if (!(is.numeric(pixel_size) && length(pixel_size) == 1L && pixel_size > 0)) {
      stop("`pixel_size` must be a single positive number, in meters.")
    }

    yr <- W$yrange
    xr <- W$xrange

    ny <- max(1L, round((yr[2] - yr[1]) / pixel_size))
    nx <- max(1L, round((xr[2] - xr[1]) / pixel_size))

    dimyx <- c(ny, nx)

  } else {
    if (!(is.numeric(dimyx) && length(dimyx) %in% c(1L, 2L))) {
      stop("`dimyx` must be length 1 or 2 numeric.")
    }

    if (length(dimyx) == 1L) {
      dimyx <- rep(as.integer(dimyx), 2L)
    }

    dimyx <- pmax(1L, as.integer(dimyx))
  }

  # KDE
  kde <- spatstat.explore::density.ppp(
    X,
    sigma = radius,
    weights = weights,
    dimyx = dimyx,
    at = "pixels",
    edge = TRUE
  )

  # Convert KDE to SpatRaster
  r_kde <- terra::rast(kde)
  terra::crs(r_kde) <- sf::st_crs(area_m$shape)$wkt

  # Mask KDE to area only when area was provided by the user
  if (crop_to_area) {
    r_kde <- terra::mask(r_kde, terra::vect(area_m$shape))
  }

  # Rescale KDE values to gradient
  vals <- terra::values(r_kde, mat = FALSE)
  vals_valid <- vals[is.finite(vals)]

  if (length(vals_valid) == 0L) {
    stop("KDE raster has no finite values after masking.")
  }

  kde_max <- max(vals_valid, na.rm = TRUE)

  if (exclude_lowest) {
    positive_vals <- vals_valid[vals_valid > 0]

    if (length(positive_vals) == 0L) {
      stop("No positive KDE values available to calculate the lower threshold.")
    }

    kde_min_threshold <- stats::quantile(
      positive_vals,
      probs = lowest_prop,
      na.rm = TRUE,
      names = FALSE
    )

  } else {
    kde_min_threshold <- 0
  }

  if (kde_max <= kde_min_threshold) {
    stop("Maximum KDE value must be greater than the minimum threshold.")
  }

  r_gradient <- r_kde

  r_gradient <- (r_gradient - kde_min_threshold) /
    (kde_max - kde_min_threshold) *
    n_classes

  # Remove values below the threshold if requested
  if (exclude_lowest) {
    r_gradient[r_kde < kde_min_threshold] <- NA
  } else {
    r_gradient[r_gradient < 0] <- 0
  }

  r_gradient[r_gradient > n_classes] <- n_classes

  names(r_gradient) <- "Gradient"

  # Convert to polygons if needed
  if (output_layer_type != "raster") {
    v_gradient <- sf::st_as_sf(
      terra::as.polygons(r_gradient, dissolve = FALSE, na.rm = TRUE)
    )

    v_gradient <- sf::st_set_crs(v_gradient, sf::st_crs(area_m$shape))
  }

  # Reproject outputs if requested
  if (return_crs == "4326") {
    if (output_layer_type %in% c("raster", "both")) {
      r_gradient <- terra::project(r_gradient, "EPSG:4326", method = "bilinear")
    }

    if (output_layer_type %in% c("shp", "both")) {
      v_gradient <- sf::st_transform(v_gradient, 4326)
    }
  }

  # Return KDE results with details as attributes
  details <- list(
    sigma = radius,
    method = radius_method,
    dimyx = dimyx,
    kde_min_threshold = kde_min_threshold,
    kde_max = kde_max,
    gradient_max = n_classes,
    exclude_lowest = exclude_lowest,
    lowest_prop = lowest_prop,
    message = paste0(
      "KDE estimated with sigma = ", radius,
      " using method = ", radius_method,
      " and grid = ", paste(dimyx, collapse = "x"),
      ". KDE values were rescaled from threshold = ", kde_min_threshold,
      " to max gradient value = ", n_classes, "."
    )
  )

  if (output_layer_type == "shp") {
    attr(v_gradient, "details") <- details
    return(v_gradient)

  } else if (output_layer_type == "raster") {
    attr(r_gradient, "details") <- details
    return(r_gradient)

  } else {
    attr(v_gradient, "details") <- details
    attr(r_gradient, "details") <- details

    return(list(
      raster = r_gradient,
      shp = v_gradient
    ))
  }
}
