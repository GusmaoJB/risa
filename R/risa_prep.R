#' Prepare KDE and overlap maps for risk assessment analysis
#'
#' @description
#' Generates KDE-based distribution and hotspot maps for species or habitats and
#' stressors, and computes pairwise overlap maps between them. Outputs can be
#' discrete reclassified maps or continuous rescaled rasters, depending on the
#' value of `continuous`.
#'
#' @details
#' The function first converts species/habitat and stressor inputs into named
#' lists of `sf` objects. If `group_x` or `group_y` are supplied, the respective
#' inputs are split into multiple layers before KDE estimation.
#'
#' If `area = NULL`, an area of interest is generated automatically from the
#' species layers, stressor layers, or their union, depending on
#' `area_strategy`. The automatic area can be based on either a convex hull or
#' a bounding box, controlled by `area_type`, and expanded using
#' `area_buffer_frac`.
#'
#' KDE maps are generated with `get_class_kernel()`. Distribution maps are
#' generated using one class, while kernel hotspot maps are generated using
#' `n_classes`. Pairwise overlap maps between each species/habitat layer and
#' each stressor layer are then generated with `get_overlap_kernel()`.
#'
#' When `continuous = FALSE`, outputs are reclassified into discrete classes and
#' can be returned as rasters, `sf` polygons, or both. When `continuous = TRUE`,
#' outputs are returned as continuous rasters and `output_layer_type` is ignored.
#'
#' When `return_crs = "metric"`, the function ensures that raster outputs are
#' in a projected CRS with square cells. If `pixel_size = NULL`, a square pixel
#' size is automatically derived from the AOI extent and `dimyx`. Outputs are
#' regridded to a common metric raster template before being returned.
#'
#' @param x Species or habitat input. Can be an `sf` object, a `data.frame`, or
#'   a list of `sf` objects. If `x` is a single `sf` or `data.frame`, it can be
#'   split into multiple layers using `group_x`.
#' @param y Stressor input. Can be an `sf` object, a `data.frame`, or a list of
#'   `sf` objects. If `y` is a single `sf` or `data.frame`, it can be split into
#'   multiple layers using `group_y`.
#' @param area Optional area of interest. Can be an `sf` polygon object, a
#'   `bbox`, a `data.frame` convertible to `sf`, or `NULL`. If `NULL`, the area
#'   is generated automatically.
#' @param n_classes Integer. Number of classes used for KDE hotspot maps and
#'   overlap maps. Default is `3`.
#' @param output_min Numeric or `NULL`. Minimum value used when
#'   `continuous = TRUE`. If `NULL`, defaults to `1`.
#' @param continuous Logical. If `FALSE`, KDE and overlap maps are returned as
#'   discrete reclassified maps. If `TRUE`, they are returned as continuous
#'   rescaled rasters. Default is `FALSE`.
#' @param output_layer_type Character. Output type for discrete maps. Options
#'   are `"raster"` (generate only `SpatRaster` maps), or `"both"` (generates `SpatRaster` and `sf` vector ["shp"] maps). Default is `"raster"`. Ignored when
#'   `continuous = TRUE`.
#' @param radius Numeric or `NULL`. KDE bandwidth in projected units. If
#'   `NULL`, the bandwidth is estimated using `radius_method`.
#' @param radius_method Character. Method used to estimate `radius` when
#'   `radius = NULL`. Options are `"nndist"`, `"nrd"`,
#'   `"std_distance_scaled"`, `"ppl"`, and `"fixed"`.
#' @param group_x,group_y Optional column names used to split `x` and `y` into
#'   multiple layers when they are not already lists.
#' @param group_size_x,group_size_y Character or `NULL`. Names of numeric
#'   columns in `x` and `y`, respectively, used as weights in KDE estimation.
#' @param pixel_size Numeric or `NULL`. Desired raster pixel size in projected
#'   units. If `NULL`, `dimyx` is used to derive grid dimensions. When
#'   `return_crs = "metric"`, a square pixel size is automatically derived if
#'   needed.
#' @param dimyx Numeric vector of length one or two. Grid dimensions used for
#'   KDE estimation when `pixel_size = NULL`. Default is `c(512, 512)`.
#' @param exclude_lowest Logical. If `TRUE`, the lowest KDE values are excluded
#'   before reclassification or rescaling. Default is `TRUE`.
#' @param lowest_prop Numeric. Proportion of the lowest KDE values to exclude
#'   when `exclude_lowest = TRUE`. Default is `0.05`.
#' @param area_strategy Character. Which input layers are used to create the
#'   automatic AOI when `area = NULL`. Options are `"stressor"`, `"species"`,
#'   or `"union"`.
#' @param area_type Character. Geometry used to create the automatic AOI.
#'   Options are `"convex_hull"` or `"bbox"`.
#' @param area_buffer_frac Numeric. Fractional expansion applied when creating
#'   the automatic AOI. Default is `0.5`.
#' @param return_crs Character. CRS of the returned outputs. `"metric"` keeps
#'   outputs in a projected metric CRS, while `"4326"` reprojects outputs to
#'   longitude/latitude WGS84.
#' @param overlap_method Character. Combination rule used to calculate overlap
#'   maps. Options are `"product"`, `"sum"`, `"geom_mean"`, and `"max"`.
#'   Passed to `get_overlap_kernel()`.
#' @param quiet Logical. If `TRUE`, suppresses progress messages. Default is
#'   `TRUE`.
#'
#' @return
#' An object of class `risaMaps`, returned as a list with:
#' \describe{
#'   \item{`species_distributions`}{A named list of species or habitat
#'   distribution maps.}
#'   \item{`stressor_distributions`}{A named list of stressor distribution
#'   maps.}
#'   \item{`species_kernel_maps`}{A named list of KDE hotspot maps for species
#'   or habitats.}
#'   \item{`stressor_kernel_maps`}{A named list of KDE hotspot maps for
#'   stressors.}
#'   \item{`overlap_maps`}{A nested list of pairwise species/habitat-by-stressor
#'   overlap maps.}
#'   \item{`area_of_interest`}{The area of interest used in the analysis,
#'   returned in the requested CRS.}
#' }
#'
#' Map elements contain a `terra::SpatRaster` in `$raster`; when
#' `continuous = FALSE` and polygon outputs are requested, they may also contain
#' an `sf` object in `$shp`.
#'
#' @importFrom sf st_as_sfc st_transform st_crs st_is_longlat st_bbox st_set_crs
#' @importFrom terra project crs res ext rast resample
#'
#' @examples
#' \dontrun{
#' # Example to be added with package sample data
#' }
#'
#' @export
risa_prep <- function(
    x, y,
    area = NULL,
    n_classes = 3,
    output_min = NULL,
    continuous = FALSE,
    output_layer_type = c("raster", "both"),
    radius = NULL,
    radius_method = c("nndist", "nrd", "std_distance_scaled", "ppl", "fixed"),
    group_x = NULL,
    group_y = NULL,
    group_size_x = NULL,
    group_size_y = NULL,
    pixel_size = NULL,
    dimyx = c(512,512),
    exclude_lowest = TRUE,
    lowest_prop = 0.05,
    area_strategy = c("stressor","species","union"),
    area_type = c("convex_hull","bbox"),
    area_buffer_frac = 0.5,
    return_crs = c("metric","4326"),
    overlap_method = c("product","sum","geom_mean","max"),
    quiet = TRUE) {

  output_layer_type <- match.arg(output_layer_type)
  radius_method <- match.arg(radius_method)
  area_strategy <- match.arg(area_strategy)
  area_type <- match.arg(area_type)
  return_crs <- match.arg(return_crs)
  overlap_method <- match.arg(overlap_method)

  # Helpers
  # Fail-fast guards for metric, square outputs
  .check_square_metric <- function(r) {
    if (!inherits(r, "SpatRaster")) return(invisible(TRUE))

    crs_str <- terra::crs(r)

    if (is.na(crs_str) || !nzchar(crs_str)) {
      stop(
        "Output raster has no CRS defined. ",
        "Stamp a metric CRS before HRA (use return_crs = 'metric')."
      )
    }

    if (terra::is.lonlat(r)) {
      stop(
        "Output raster is lon/lat but HRA requires metric CRS. ",
        "Use return_crs = 'metric'."
      )
    }

    rs <- as.numeric(terra::res(r))

    if (length(rs) < 2L) {
      stop("Raster resolution is invalid (length < 2).")
    }

    if (!isTRUE(all.equal(rs[1], rs[2]))) {
      stop(sprintf(
        paste0(
          "Output raster has non-square cells ",
          "(xres = %.6f, yres = %.6f). ",
          "Set pixel_size or adjust dimyx."
        ),
        rs[1], rs[2]
      ))
    }

    invisible(TRUE)
  }

  # Force all kernel outputs to have a $raster element
  .as_raster_list_item <- function(obj) {
    if (inherits(obj, "SpatRaster")) {
      return(list(raster = obj))
    }

    if (is.list(obj) && "raster" %in% names(obj)) {
      return(obj)
    }

    stop("Kernel output must be a SpatRaster or a list containing a `raster` element.")
  }

  .standardize_kernel_list <- function(lst) {
    lapply(lst, .as_raster_list_item)
  }

  # Normalize inputs to sf lists
  spp_list <- as_sf_list(x, group = group_x, label_prefix = "sp")
  str_list <- as_sf_list(y, group = group_y, label_prefix = "stressor")

  # Choose or build AOI
  if (is.null(area)) {
    if (!quiet) message("No area provided; creating AOI from ", area_strategy, " layers (", area_type, ").")
    src <- switch(area_strategy,
                  "stressor" = merge_shp(str_list),
                  "species" = merge_shp(spp_list),
                  "union" = merge_shp(c(spp_list, str_list)))
    area <- create_area(src, area_type = area_type, buffer_frac = area_buffer_frac, quiet = quiet)
  } else if (inherits(area, "bbox")) {
    area <- sf::st_as_sfc(area)
  } else if (is.data.frame(area)) {
    area <- df_to_shp(area)
  } else if (!inherits(area, "sf")) {
    stop("`area` must be NULL, bbox, data.frame, or sf.")
  }

  # Ensure AOI has a CRS
  if (is.na(sf::st_crs(area))) {
    area <- sf::st_set_crs(area, guess_crs(area))
  }

  # If metric outputs are requested, ensure AOI is in a metric CRS
  is_ll <- suppressWarnings(sf::st_is_longlat(area))
  if (isTRUE(is_ll) && return_crs == "metric") {
    area_metric <- transform_to_metric(area)
  } else {
    area_metric <- area
  }

  # Auto-derive a square pixel_size (meters) if needed
  if (return_crs == "metric" && is.null(pixel_size)) {
    bb <- sf::st_bbox(area_metric)
    dx <- as.numeric(bb["xmax"] - bb["xmin"])
    dy <- as.numeric(bb["ymax"] - bb["ymin"])
    if (!is.finite(dx) || !is.finite(dy) || dx <= 0 || dy <= 0) {
      stop("AOI has zero or invalid extent; cannot derive pixel_size. Check 'area'.")
    }
    pixel_size <- min(dx / dimyx[2], dy / dimyx[1])
    if (!is.finite(pixel_size) || pixel_size <= 0) {
      stop("Derived pixel_size non-positive; adjust 'dimyx' or provide 'pixel_size'.")
    }
    if (!quiet) message(sprintf("Auto pixel_size (square): %.3f m", pixel_size))
  }

  # Kernels builder (passes pixel_size forward)
  if (!quiet) message("Buiding KDE maps...")
  # When continuous = TRUE, force output_layer_type to "raster"
  kernel_output_type <- if (continuous) "raster" else output_layer_type

  build_kernels <- function(lst, group_size = NULL, ncls = n_classes) {
    nms <- names(lst)

    if (is.null(nms)) {
      nms <- paste0("item_", seq_along(lst))
    } else {
      nms[nms == ""] <- paste0("item_", which(nms == ""))
    }

    out <- Map(function(item, nm) {

      if (!quiet) message("Processing: ", nm)

      get_class_kernel(
        item,
        area = if (return_crs == "metric") area_metric else area,
        n_classes = ncls,
        output_min = output_min,
        output_layer_type = kernel_output_type,
        radius = radius,
        radius_method = radius_method,
        group_size = group_size,
        pixel_size = pixel_size,
        dimyx = dimyx,
        exclude_lowest = exclude_lowest,
        lowest_prop = lowest_prop,
        continuous = continuous,
        return_crs = return_crs,
        quiet = quiet
      )

    }, lst, nms)

    names(out) <- nms
    out
  }

  map_pa <- function(item) {

    if (inherits(item, "SpatRaster")) {
      out <- terra::ifel(is.na(item), NA, 1)
      names(out) <- names(item)
      return(out)
    }

    if (inherits(item, "sf")) {
      geom_union <- sf::st_union(sf::st_geometry(item))

      out <- sf::st_sf(
        Rating = 1,
        geometry = sf::st_sfc(geom_union, crs = sf::st_crs(item))
      )

      return(out)
    }

    item
  }

  build_pa_list <- function(x) {
    lapply(x, function(item) {

      out <- list()

      if ("raster" %in% names(item)) {
        out$raster <- map_pa(item$raster)
      }

      if ("shp" %in% names(item)) {
        out$shp <- map_pa(item$shp)
      }

      out
    })
  }

  # Performing KDEs with n_classes
  spp_kernel_list <- build_kernels(spp_list, group_size = group_size_x, ncls = n_classes)
  stressor_kernel_list <- build_kernels(str_list, group_size = group_size_y, ncls = n_classes)

  # Standardize lists when continuous = TRUE
  spp_kernel_list <- .standardize_kernel_list(spp_kernel_list)
  stressor_kernel_list <- .standardize_kernel_list(stressor_kernel_list)

  # Create distribution lists
  spp_distribution_list <- build_pa_list(spp_kernel_list)
  str_distribution_list <- build_pa_list(stressor_kernel_list)

  # AOI CRS for output (only transform AOI geometry for return)
  if (return_crs == "4326") {
    area_out <- sf::st_transform(area, 4326)
  } else {
    area_out <- area_metric
  }

  # Overlap maps
  if (!quiet) message("Generating overlap maps...")

  # Force overlap output type to raster when continuous = TRUE
  overlap_output_type <- if (continuous) "raster" else output_layer_type

  overlap_maps_list <- lapply(spp_kernel_list, function(sp) {
    lapply(stressor_kernel_list, function(st) {

      sp_rast <- sp$raster
      st_rast <- st$raster

      ol <- get_overlap_kernel(
        sp_rast,
        st_rast,
        n_classes = n_classes,
        continuous = continuous,
        output_min = output_min,
        out_classes = n_classes,
        method = overlap_method,
        output_layer_type = overlap_output_type,
        quiet = quiet
      )

      # Standardize overlap output structure
      ol <- .as_raster_list_item(ol)

      if (return_crs == "4326") {

        ol$raster <- terra::project(
          ol$raster,
          "EPSG:4326",
          method = if (continuous) "bilinear" else "near"
        )

        if ("shp" %in% names(ol)) {
          ol$shp <- sf::st_transform(ol$shp, 4326)
        }
      }

      ol
    })
  })

  # Build a square template in the metric CRS and regrid outputs to it
  if (return_crs == "metric") {
    crs_metric_wkt <- sf::st_crs(area_metric)$wkt
    if (is.na(crs_metric_wkt) || !nzchar(crs_metric_wkt)) {
      stop("Target metric CRS is missing; cannot stamp CRS on outputs.")
    }

    .make_template <- function(aoi_sf, crs_wkt, res_m) {
      bb  <- sf::st_bbox(aoi_sf)
      ext <- terra::ext(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"])
      terra::rast(extent = ext, crs = crs_wkt, resolution = c(res_m, res_m))
    }
    .regrid_to_template <- function(r, template) {
      if (!inherits(r, "SpatRaster")) return(r)

      resample_method <- if (continuous) "bilinear" else "near"

      r2 <- if (!identical(terra::crs(r), terra::crs(template))) {
        terra::project(r, template, method = resample_method)
      } else {
        r
      }

      terra::resample(r2, template, method = resample_method)
    }

    template <- .make_template(area_metric, crs_metric_wkt, pixel_size)

    # Regrid distributions / kernels
    for (nm in names(spp_distribution_list)) {
      spp_distribution_list[[nm]]$raster <- .regrid_to_template(
        spp_distribution_list[[nm]]$raster,
        template
      )
    }

    for (nm in names(str_distribution_list)) {
      str_distribution_list[[nm]]$raster <- .regrid_to_template(
        str_distribution_list[[nm]]$raster,
        template
      )
    }

    for (nm in names(spp_kernel_list)) {
      spp_kernel_list[[nm]]$raster <- .regrid_to_template(
        spp_kernel_list[[nm]]$raster,
        template
      )
    }

    for (nm in names(stressor_kernel_list)) {
      stressor_kernel_list[[nm]]$raster <- .regrid_to_template(
        stressor_kernel_list[[nm]]$raster,
        template
      )
    }

    # Regrid overlaps
    for (i in names(overlap_maps_list)) {
      for (j in names(overlap_maps_list[[i]])) {
        overlap_maps_list[[i]][[j]]$raster <- .regrid_to_template(
          overlap_maps_list[[i]][[j]]$raster,
          template
        )
      }
    }

    # Stamp CRS on any rasters missing it
    .set_crs_if_missing <- function(r) {
      if (inherits(r, "SpatRaster")) {
        crs_str <- terra::crs(r)

        if (is.na(crs_str) || !nzchar(crs_str)) {
          terra::crs(r) <- crs_metric_wkt
        }
      }

      r
    }

    for (nm in names(spp_distribution_list)) {
      spp_distribution_list[[nm]]$raster <- .set_crs_if_missing(
        spp_distribution_list[[nm]]$raster
      )
    }

    for (nm in names(str_distribution_list)) {
      str_distribution_list[[nm]]$raster <- .set_crs_if_missing(
        str_distribution_list[[nm]]$raster
      )
    }

    for (nm in names(spp_kernel_list)) {
      spp_kernel_list[[nm]]$raster <- .set_crs_if_missing(
        spp_kernel_list[[nm]]$raster
      )
    }

    for (nm in names(stressor_kernel_list)) {
      stressor_kernel_list[[nm]]$raster <- .set_crs_if_missing(
        stressor_kernel_list[[nm]]$raster
      )
    }

    for (i in names(overlap_maps_list)) {
      for (j in names(overlap_maps_list[[i]])) {
        overlap_maps_list[[i]][[j]]$raster <- .set_crs_if_missing(
          overlap_maps_list[[i]][[j]]$raster
        )
      }
    }

    # Enforce square + metric
    for (nm in names(spp_distribution_list)) {
      .check_square_metric(spp_distribution_list[[nm]]$raster)
    }

    for (nm in names(str_distribution_list)) {
      .check_square_metric(str_distribution_list[[nm]]$raster)
    }

    for (nm in names(spp_kernel_list)) {
      .check_square_metric(spp_kernel_list[[nm]]$raster)
    }

    for (nm in names(stressor_kernel_list)) {
      .check_square_metric(stressor_kernel_list[[nm]]$raster)
    }

    for (i in names(overlap_maps_list)) {
      for (j in names(overlap_maps_list[[i]])) {
        .check_square_metric(overlap_maps_list[[i]][[j]]$raster)
      }
    }

    # Ensure overlap maps include raster outputs
    for (i in names(overlap_maps_list)) {
      for (j in names(overlap_maps_list[[i]])) {
        has_raster <- is.list(overlap_maps_list[[i]][[j]]) &&
          "raster" %in% names(overlap_maps_list[[i]][[j]]) &&
          inherits(overlap_maps_list[[i]][[j]]$raster, "SpatRaster")

        if (!has_raster) {
          stop("Overlap map missing raster output.")
        }
      }
    }
  }

  # Assemble
  out <- list(
    species_distributions = spp_distribution_list,
    stressor_distributions = str_distribution_list,
    species_kernel_maps    = spp_kernel_list,
    stressor_kernel_maps   = stressor_kernel_list,
    overlap_maps           = overlap_maps_list,
    area_of_interest       = area_out
  )
  class(out) <- c("risaMaps", class(out))
  out
}
