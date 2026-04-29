get_class_kernel2 <- function(
    x,
    area = NULL,
    n_classes = 3,
    output_min = NULL,
    output_layer_type = c("shp", "raster", "both"),
    radius = NULL,
    radius_method = c("nndist", "ppl", "fixed"),
    group_size = NULL,
    pixel_size = NULL,
    dimyx = c(512, 512),
    exclude_lowest = TRUE,
    lowest_prop = 0.05,
    continuous = FALSE,
    input_crs = NULL,
    return_crs = c("metric", "4326"),
    quiet = TRUE) {

  output_layer_type <- match.arg(output_layer_type)
  radius_method <- match.arg(radius_method)
  return_crs <- match.arg(return_crs)

  # Helper: apply or check input CRS
  apply_input_crs <- function(obj, input_crs = NULL, object_name = "object") {

    if (is.null(input_crs)) {
      return(obj)
    }

    input_crs_obj <- sf::st_crs(input_crs)

    if (is.na(input_crs_obj)) {
      stop("`input_crs` could not be interpreted as a valid CRS.")
    }

    current_crs <- sf::st_crs(obj)

    if (is.na(current_crs)) {

      if (!quiet) {
        message(
          "Assigning input CRS to `", object_name, "`: ",
          input_crs_obj$input
        )
      }

      obj <- sf::st_set_crs(obj, input_crs_obj)

    } else if (current_crs != input_crs_obj) {

      warning(
        "`", object_name, "` already has a CRS different from `input_crs`. ",
        "Ignoring argument `input_crs` ...",
      )
    }

    obj
  }

  # Normalize input to sf
  if (inherits(x, "sf")) {

    x_sf <- x
    x_sf <- apply_input_crs(x_sf, input_crs, object_name = "x")

  } else if (inherits(x, c("data.frame", "tbl_df", "tbl"))) {

    if (ncol(x) < 2 || !all(vapply(x[, 1:2], is.numeric, logical(1)))) {
      stop("For data.frame input, the first two columns must be numeric coordinates.")
    }

    if (!quiet) message("Converting input data.frame to sf...")

    if (!is.null(input_crs)) {

      x_sf <- sf::st_as_sf(
        x,
        coords = names(x)[1:2],
        crs = input_crs,
        remove = FALSE
      )

    } else {

      x_sf <- df_to_shp(x)
    }

  } else {
    stop("`x` must be an sf object or a data.frame with coordinate columns.")
  }

  if (is.na(sf::st_crs(x_sf))) {
    stop(
      "`x` has no CRS. Please provide `input_crs`, for example ",
      "`input_crs = 4326` for longitude/latitude or ",
      "`input_crs = 32722` for WGS 84 / UTM zone 22S."
    )
  }

  # Helper: create observation-based rectangular window
  make_obs_window <- function(coords,
                              buffer_frac = 0.5,
                              buffer_dist = NULL) {

    xr <- range(coords[, 1], na.rm = TRUE)
    yr <- range(coords[, 2], na.rm = TRUE)

    x_width <- diff(xr)
    y_width <- diff(yr)

    # Handle cases with very few points or identical x/y coordinates
    if (x_width == 0) {
      x_width <- 1
      xr <- xr + c(-0.5, 0.5)
    }

    if (y_width == 0) {
      y_width <- 1
      yr <- yr + c(-0.5, 0.5)
    }

    # Relative buffer based on coordinate spread
    x_buffer_frac <- x_width * buffer_frac
    y_buffer_frac <- y_width * buffer_frac

    # Handle absolute buffer in projected units
    if (!is.null(buffer_dist)) {
      x_buffer <- max(x_buffer_frac, buffer_dist)
      y_buffer <- max(y_buffer_frac, buffer_dist)
    } else {
      x_buffer <- x_buffer_frac
      y_buffer <- y_buffer_frac
    }

    spatstat.geom::owin(
      xrange = c(xr[1] - x_buffer, xr[2] + x_buffer),
      yrange = c(yr[1] - y_buffer, yr[2] + y_buffer)
    )
  }

  # Area handling and projection
  if (is.null(area)) {

    if (!quiet) {
      message(
        "No area provided. KDE will use an observation-based window ",
        "and will not be masked or cropped."
      )
    }

    x_m <- transform_to_metric(x_sf, quiet = quiet)

    crs_ref <- sf::st_crs(x_m$shape)
    coords <- x_m$coordinates[, 1:2, drop = FALSE]

    W <- make_obs_window(coords, buffer_frac = 0.5)

    use_area_mask <- FALSE
    area_m <- NULL

  } else {

    if (inherits(area, "bbox")) {

      area <- sf::st_as_sf(sf::st_as_sfc(area))

      # bbox inherits CRS only if it had one originally
      area <- apply_input_crs(area, input_crs, object_name = "area")

    } else if (inherits(area, "data.frame")) {

      area <- df_to_shp(area)
      area <- apply_input_crs(area, input_crs, object_name = "area")

    } else if (inherits(area, "sf")) {

      area <- apply_input_crs(area, input_crs, object_name = "area")

    } else {
      stop("`area` must be NULL, bbox, data.frame, or sf.")
    }

    if (is.na(sf::st_crs(area))) {
      stop(
        "`area` has no CRS. Please provide `input_crs` or assign a CRS to `area`."
      )
    }

    area_m <- transform_to_metric(area, quiet = quiet)

    crs_ref <- sf::st_crs(area_m$shape)
    crs_m <- crs_ref$epsg

    x_m <- transform_to_metric(
      x_sf,
      metric_crs = crs_m,
      quiet = quiet
    )

    coords <- x_m$coordinates[, 1:2, drop = FALSE]

    W <- spatstat.geom::as.owin(area_m$shape)

    use_area_mask <- TRUE
  }

  if (nrow(coords) == 0L) {
    stop("No points available after preprocessing.")
  }

  if (nrow(unique(coords)) == 1L && is.null(radius)) {
    warning(
      "Only one unique point; KDE will be a single peak. ",
      "Consider setting `radius` explicitly."
    )
  }

  # Handling weights
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

  # Initial point pattern and bandwidth
  X <- spatstat.geom::ppp(
    x = coords[, 1],
    y = coords[, 2],
    window = W,
    check = FALSE
  )

  if (is.null(radius)) {

    if (radius_method == "nndist") {

      nn <- spatstat.geom::nndist(X)
      radius <- 1.5 * mean(nn[is.finite(nn)], na.rm = TRUE)

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
      stop("`radius` is NULL but `radius_method = 'fixed'. Supply a numeric radius.")
    }
  }

  if (!is.numeric(radius) || length(radius) != 1L || !(radius > 0)) {
    stop("`radius` must be a single positive number in projected units.")
  }

  # Expand window when area = NULL
  if (!use_area_mask) {

    kernel_buffer <- if (nrow(coords) < 20) {
      4 * radius
    } else if (nrow(coords) < 50) {
      3.5 * radius
    } else {
      3 * radius
    }

    W <- make_obs_window(
      coords = coords,
      buffer_frac = 0.5,
      buffer_dist = kernel_buffer
    )

    X <- spatstat.geom::ppp(
      x = coords[, 1],
      y = coords[, 2],
      window = W,
      check = FALSE
    )

    if (!quiet) {
      message(sprintf(
        "Observation-based KDE window expanded using buffer = %.3f projected units.",
        kernel_buffer
      ))
    }
  }

  # Handling grid resolution
  if (!is.null(pixel_size)) {

    if (!(is.numeric(pixel_size) && length(pixel_size) == 1L && pixel_size > 0)) {
      stop("`pixel_size` must be a single positive number in projected units.")
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

  # Convert to SpatRaster
  r_kde <- terra::rast(kde)
  terra::crs(r_kde) <- crs_ref$wkt

  # Create reclassification matrix first
  re_mat <- reclass_matrix(
    r_kde,
    n_classes = n_classes,
    exclude_lowest = exclude_lowest,
    lowest_prop = lowest_prop
  )

  # Reclassify or rescale KDE values
  if (continuous) {

    # Use the reclassification matrix to identify valid KDE classes
    r_valid <- terra::classify(
      r_kde,
      re_mat,
      include.lowest = TRUE
    )

    # Remove cells classified as NA by the reclassification matrix
    r_kde[is.na(r_valid)] <- NA

    # Get remaining KDE values
    kde_vals <- terra::values(r_kde, mat = FALSE)
    kde_vals <- kde_vals[is.finite(kde_vals)]

    if (length(kde_vals) == 0L) {
      stop("KDE raster has no finite values to rescale.")
    }

    # Use the reclassification matrix to define the lower threshold
    valid_rows <- !is.na(re_mat[, 3])

    lower_threshold <- min(re_mat[valid_rows, 1], na.rm = TRUE)
    upper_threshold <- max(kde_vals, na.rm = TRUE)

    if (upper_threshold <= lower_threshold) {
      stop(
        "Cannot rescale KDE values because the maximum KDE value is not ",
        "greater than the lower threshold."
      )
    }

    # Define output scale
    output_min <- if (is.null(output_min)) 1 else output_min
    output_max <- n_classes

    # Rescale KDE values using the same lower threshold used for classification
    r_cls <- ((r_kde - lower_threshold) /
                (upper_threshold - lower_threshold)) *
      (output_max - output_min) + output_min

  } else {

    # Use the same matrix directly for categorical reclassification
    r_cls <- terra::classify(
      r_kde,
      re_mat,
      include.lowest = TRUE
    )
  }

  # Crop/mask KDEs when an area is provided
  if (use_area_mask) {
    r_cls <- terra::mask(r_cls, terra::vect(area_m$shape))
  }

  names(r_cls) <- "Rating"

  # Convert to polygons if needed
  if (output_layer_type != "raster") {

    v_cls <- sf::st_as_sf(
      terra::as.polygons(r_cls, dissolve = TRUE)
    )

    v_cls <- sf::st_set_crs(v_cls, crs_ref)
  }

  # Reproject outputs if requested
  if (return_crs == "4326") {

    if (output_layer_type %in% c("raster", "both")) {
      r_cls <- terra::project(
        r_cls,
        "EPSG:4326",
        method = ifelse(continuous, "bilinear", "near")
      )
    }

    if (output_layer_type %in% c("shp", "both")) {
      v_cls <- sf::st_transform(v_cls, 4326)
    }
  }

  # Details
  details <- list(
    sigma = radius,
    method = radius_method,
    dimyx = dimyx,
    input_crs = sf::st_crs(x_sf)$input,
    output_crs = sf::st_crs(crs_ref)$input,
    area_provided = use_area_mask,
    n_observations = nrow(coords),
    message = paste0(
      "KDE estimated with sigma = ", radius,
      " using method = ", radius_method,
      " and grid = ", paste(dimyx, collapse = "x"),
      ifelse(
        use_area_mask,
        ". Output was masked to provided area.",
        ". No area was provided; output was not masked or cropped."
      )
    )
  )

  # Return
  if (output_layer_type == "shp") {

    attr(v_cls, "details") <- details
    return(v_cls)

  } else if (output_layer_type == "raster") {

    attr(r_cls, "details") <- details
    return(r_cls)

  } else {

    attr(v_cls, "details") <- details
    attr(r_cls, "details") <- details

    return(list(
      raster = r_cls,
      shp = v_cls
    ))
  }
}
