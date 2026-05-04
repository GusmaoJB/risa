#' Create reclassified or rescaled KDE hotspot maps
#'
#' @description
#' Creates kernel density estimation (KDE) maps from point data and converts the
#' resulting density surface into either categorical classes or a continuous
#' rescaled rating. The input can be an `sf` point object or a `data.frame` whose
#' first two columns contain numeric coordinates. The output can be returned as
#' polygons, as a raster, or as both.
#'
#' @details
#' This function estimates a two-dimensional kernel density surface using
#' [spatstat.explore::density.ppp()]. Input points are first converted to an
#' `sf` object when necessary and projected to a metric coordinate reference
#' system so that the KDE bandwidth, grid resolution, and optional pixel size are
#' expressed in projected map units, usually meters.
#'
#' If an `area` is supplied, it is used to define the observation window for the
#' point pattern. The resulting KDE raster is then masked to this area. If
#' `area = NULL`, the function creates an observation-based rectangular window
#' from the coordinate range of the input points. In this case, the KDE is not
#' cropped or masked, and the window is expanded using a buffer based on the
#' estimated or supplied bandwidth.
#'
#' The KDE bandwidth is controlled by the `radius` argument. If `radius` is
#' supplied, this value is used directly as the kernel bandwidth (`sigma`). If
#' `radius = NULL`, the bandwidth is estimated automatically using one of the
#' methods selected in `radius_method`:
#'
#' \describe{
#'   \item{`"nndist"`}{
#'   Uses the mean nearest-neighbour distance among points, multiplied by 1.5.
#'   This is a simple local-spacing rule that tends to produce bandwidths
#'   related to the observed point distribution.
#'   }
#'
#'   \item{`"nrd"`}{
#'   Uses a normal-reference-style rule based on the pooled standard deviation
#'   of the x and y coordinates, multiplied by \eqn{n^{-1/6}}, where \eqn{n}
#'   is the number of observations. This method can produce relatively smooth
#'   KDE surfaces when points are widely dispersed.
#'   }
#'
#'   \item{`"std_distance_scaled"`}{
#'   Uses the standard distance of the points around their mean center, scaled
#'   by \eqn{(2 / (3n))^{1/4}}, where \eqn{n} is the number of observations.
#'   This provides a bandwidth based on the spatial dispersion of the point
#'   cloud while reducing the bandwidth as sample size increases.
#'   }
#'
#'   \item{`"ppl"`}{
#'   Uses the likelihood cross-validation bandwidth selected by
#'   [spatstat.explore::bw.ppl()]. This method estimates the bandwidth from
#'   the point pattern using a point-process likelihood criterion.
#'   }
#'
#'   \item{`"fixed"`}{
#'   If `radius` is supplied, it overrides `radius_method` and is used directly.
#'   The `"fixed"` option is intended to make this behavior explicit and requires
#'   `radius` to be supplied.
#'   }
#' }
#'
#' The KDE surface can be converted into categorical classes or into a continuous
#' rating. When `continuous = FALSE`, KDE values are reclassified into
#' `n_classes` ordered classes (starting at 1 or any other value defined in
#' `output_min`). By default, the lowest KDE values can be excluded
#' using `exclude_lowest = TRUE` and `lowest_prop = 0.05`, which removes the
#' lowest proportion of KDE values before assigning classes. When
#' `continuous = TRUE`, `n_classes` is interpreted as the upper value of
#' the output rating scale rather than as a number of discrete classes. The lower
#' value is defined as 1, but can be changed by `output_min`.
#'
#' If `group_size` is provided, the corresponding numeric column in `x` is used
#' as weights in the KDE estimation. This allows observations to contribute
#' differently to the estimated density surface, for example when each point
#' represents a group of individuals.
#'
#' @param x An `sf` point object or a `data.frame`/`tibble` containing point
#'   coordinates. If a `data.frame` is supplied, the first two columns must be
#'   numeric coordinate columns.
#' @param area Optional spatial object used to define the KDE observation window
#'   and mask the final output. Ideally an `sf` object, a `bbox`, or a
#'   `data.frame` that can be converted to an `sf` object. If `NULL`, an
#'   observation-based rectangular window is created and the output is not
#'   masked.
#' @param n_classes Integer. Number of output classes or the maximum value of
#'   the continuous rating scale. Default is `3`.
#' @param output_min Numeric or `NULL`. Minimum value used when
#'   `continuous = TRUE`. If `NULL`, defaults to `1`.
#' @param output_layer_type Character. Type of output layer to return. Options
#'   are `"shp"` for polygons, `"raster"` for a `SpatRaster`, or `"both"` for
#'   a list containing both outputs.
#' @param radius Numeric or `NULL`. KDE bandwidth, passed as `sigma` to
#'   [spatstat.explore::density.ppp()]. The value must be a single positive
#'   number in projected units. If `NULL`, the bandwidth is estimated using
#'   `radius_method`.
#' @param radius_method Character. Method used to estimate `radius` when
#'   `radius = NULL`. Options are `"nndist"`, `"nrd"`,
#'   `"std_distance_scaled"`, `"ppl"`, and `"fixed"`. See Details.
#' @param group_size Character or `NULL`. Name of a numeric column in `x` to use
#'   as weights in the KDE estimation. If `NULL`, all points receive equal
#'   weight.
#' @param pixel_size Numeric or `NULL`. Desired raster pixel size in projected
#'   units. If supplied, it is used to derive the KDE grid dimensions from the
#'   observation window. If `NULL`, `dimyx` is used.
#' @param dimyx Numeric vector of length one or two. Grid dimensions passed to
#'   [spatstat.explore::density.ppp()]. If length one, the same value is used
#'   for both dimensions. Default is `c(512, 512)`.
#' @param exclude_lowest Logical. If `TRUE`, the lowest KDE values are excluded
#'   before reclassification or rescaling. Default is `TRUE`.
#' @param lowest_prop Numeric. Proportion of the lowest KDE values to exclude
#'   when `exclude_lowest = TRUE`. Default is `0.05`.
#' @param continuous Logical. If `FALSE`, the KDE raster is reclassified into
#'   ordered integer classes. If `TRUE`, KDE values are rescaled to a continuous
#'   rating scale. Default is `FALSE`.
#' @param input_crs Coordinate reference system to assign to `x` and/or `area`
#'   when they do not already have one. Can be any value accepted by
#'   [sf::st_crs()], such as an EPSG code.
#' @param return_crs Character. CRS of the returned output. `"metric"` keeps the
#'   projected metric CRS used for KDE estimation, while `"4326"` reprojects the
#'   output to longitude/latitude WGS84.
#' @param quiet Logical. If `FALSE`, progress messages are printed. Default is
#'   `TRUE`.
#'
#' @return
#' Depending on `output_layer_type`, returns:
#' \describe{
#'   \item{`"shp"`}{An `sf` polygon object containing the reclassified or
#'   rescaled KDE ratings.}
#'   \item{`"raster"`}{A `terra::SpatRaster` containing the reclassified or
#'   rescaled KDE ratings.}
#'   \item{`"both"`}{A list with two elements: `raster`, a `terra::SpatRaster`,
#'   and `shp`, an `sf` polygon object.}
#' }
#'
#' The returned object includes a `details` attribute with the estimated or
#' supplied bandwidth, selected bandwidth method, grid dimensions, input and
#' output CRS information, whether an area mask was used, the number of
#' observations, and a short processing message.
#'
#' @examples
#' \dontrun{
#' # Example using a data.frame with coordinates in the first two columns
#' kde_map <- get_class_kernel(
#'   x = points_df,
#'   area = study_area,
#'   n_classes = 3,
#'   radius_method = "std_distance_scaled",
#'   output_layer_type = "both",
#'   input_crs = 32722
#' )
#'
#' # Example using a fixed bandwidth
#' kde_fixed <- get_class_kernel(
#'   x = points_sf,
#'   area = study_area,
#'   radius = 1000,
#'   radius_method = "fixed",
#'   output_layer_type = "raster",
#'   return_crs = "4326"
#' )
#'
#' # Example returning a continuous rating instead of classes
#' kde_cont <- get_class_kernel(
#'   x = points_sf,
#'   area = study_area,
#'   continuous = TRUE,
#'   n_classes = 3,
#'   output_min = 1
#' )
#' }
#'
#' @export
get_class_kernel <- function(
    x,
    area = NULL,
    n_classes = 3,
    output_min = NULL,
    output_layer_type = c("shp", "raster", "both"),
    radius = NULL,
    radius_method = c("nndist", "nrd", "std_distance_scaled", "ppl", "fixed"),
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

    if (radius_method == "nrd") {

      sd_xy <- sqrt(
        (stats::var(coords[, 1], na.rm = TRUE) +
           stats::var(coords[, 2], na.rm = TRUE)) / 2
      )

      radius <- sd_xy * nrow(coords)^(-1 / 6)

      if (!quiet) {
        message(sprintf("Auto bandwidth (nrd): sigma = %.3f", radius))
      }

    } else if (radius_method == "nndist") {

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

    } else if (radius_method == "std_distance_scaled") {

      mean_x <- mean(coords[, 1], na.rm = TRUE)
      mean_y <- mean(coords[, 2], na.rm = TRUE)

      std_dist <- sqrt(
        mean(
          (coords[, 1] - mean_x)^2 +
            (coords[, 2] - mean_y)^2,
          na.rm = TRUE
        )
      )

      n <- nrow(coords)

      radius <- ((2 / (3 * n))^(1 / 4)) * std_dist

      if (!quiet) {
        message(sprintf(
          "Auto bandwidth (std_distance_scaled): sigma = %.3f",
          radius
        ))
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
    output_crs = if (return_crs == "4326") {
      sf::st_crs(4326)$input
      } else {
      sf::st_crs(crs_ref)$input
      },
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
