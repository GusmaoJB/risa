#' Prepare and generate maps for risk assessment analysis
#'
#' Generates kernel density maps for species (habitats) and stressors, plus overlap maps,
#' within a given area of interest (auto-generated from stressors' distributions or user-supplied).
#'
#' @param x Species (Habitat) input: `sf`, data.frame, or a list of `sf`. If a data.frame/sf,
#'   you can split it by `group_x` into multiple species layers.
#' @param y Stressor input: `sf`, data.frame, or a list of `sf`. If a data.frame/sf,
#'   you can split it by `group_y` into multiple stressor layers.
#' @param area Optional AOI polygon (`sf`) or `bbox`/data.frame; if `NULL`, it's computed.
#' @param n_classes Integer number of classes for kernel reclassification (default 3).
#' @param output_layer_type One of `"shp"`, `"raster"`, or `"both"` (default `"both"`).
#' @param radius KDE bandwidth (projected units). If `NULL`, uses `radius_method`.
#' @param radius_method One of `"nndist"` (default), `"ppl"`, or `"fixed"`.
#' @param group_x,group_y Optional grouping columns to split `x`/`y` when they aren't lists.
#' @param group_size_x,group_size_y Optional numeric weight columns for KDE.
#' @param pixel_size Optional pixel size (meters) for the KDE grid; else `dimyx` is used.
#' @param dimyx Grid size (ny, nx) for KDE when `pixel_size` is `NULL` (default c(512,512)).
#' @param exclude_lowest,lowest_prop Passed to `reclass_matrix()` (defaults TRUE, 0.05).
#' @param area_strategy `"stressor"` (default), `"species"`, or `"union"` for auto AOI.
#' @param area_type `"convex_hull"` (default) or `"bbox"` for auto AOI.
#' @param area_buffer_frac Fractional expansion for auto AOI (default 0.5).
#' @param return_crs `"metric"` (default) or `"4326"` to reproject outputs to WGS84.
#' @param overlap_method Combination rule for overlap: one of `"product"`, `"sum"`,
#'   `"geom_mean"`, or `"max"` (default `"product"`). Passed to `get_overlap_kernel()`.
#' @param quiet Suppress messages (default TRUE).
#' @return A list with: `species_distributions`, `stressor_distributions`,
#'   `species_kernel_maps`, `stressor_kernel_maps`, `overlap_maps`, and `area_of_interest`.
#' @importFrom sf st_as_sfc st_transform st_crs st_is_longlat st_bbox st_set_crs
#' @importFrom terra project crs res
#' @examples
#' # add example here
#' @export
risa_prep <- function(
    x, y,
    area = NULL,
    n_classes = 3,
    output_layer_type = c("both","shp","raster"),
    radius = NULL,
    radius_method = c("nndist","ppl","fixed"),
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
      stop("Output raster has no CRS defined. Stamp a metric CRS before HRA (use return_crs='metric').")
    }
    info <- terra::crs(r, describe = TRUE)
    is_ll <- tryCatch(isTRUE(info$is_lonlat), error = function(e) NA)
    if (isTRUE(is_ll)) stop("Output raster is lon/lat but HRA requires metric CRS. Use return_crs='metric'.")
    rs <- as.numeric(terra::res(r))
    if (length(rs) < 2L) stop("Raster resolution is invalid (length < 2).")
    if (!isTRUE(all.equal(rs[1], rs[2]))) {
      stop(sprintf("Output raster has non-square cells (xres=%.6f, yres=%.6f). Set pixel_size or adjust dimyx.",
                   rs[1], rs[2]))
    }
    # Optional: catch rotated/sheared transforms (can also trigger GDAL warnings)
    gt <- tryCatch(terra::geotransform(r), error = function(e) NULL)
    if (!is.null(gt) && (gt[3] != 0 || gt[5] != 0)) {
      stop("Raster has rotation/shear; reproject/resample to a north-up grid.")
    }
    invisible(TRUE)
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
  build_kernels <- function(lst, group_size = NULL, ncls = n_classes) {
    lapply(lst, function(item) {
      get_class_kernel(
        item,
        area = if (return_crs == "metric") area_metric else area,
        n_classes = ncls,
        output_layer_type = "both",
        radius = radius,
        radius_method = radius_method,
        group_size = group_size,
        pixel_size = pixel_size,
        dimyx = dimyx,
        exclude_lowest = exclude_lowest,
        lowest_prop = lowest_prop,
        return_crs = return_crs,
        quiet = quiet
      )
    })
  }

  # Distributions (presence maps): 1 class; Kernels: n_classes
  spp_distribution_list  <- build_kernels(spp_list, group_size = group_size_x, ncls = 1L)
  str_distribution_list  <- build_kernels(str_list, group_size = group_size_y, ncls = 1L)
  spp_kernel_list        <- build_kernels(spp_list, group_size = group_size_x, ncls = n_classes)
  stressor_kernel_list   <- build_kernels(str_list, group_size = group_size_y, ncls = n_classes)

  # AOI CRS for output (only transform AOI geometry for return)
  if (return_crs == "4326") {
    area_out <- sf::st_transform(area, 4326)
  } else {
    area_out <- area_metric
  }

  # Overlap maps
  if (!quiet) message("Generating overlap maps...")
  overlap_maps_list <- lapply(spp_kernel_list, function(sp) {
    lapply(stressor_kernel_list, function(st) {
      ol <- get_overlap_kernel(
        sp$raster, st$raster,
        n_classes = n_classes,
        method = overlap_method,
        output_layer_type = output_layer_type,
        quiet = quiet
      )
      if (return_crs == "4326") {
        if (inherits(ol, "SpatRaster")) {
          ol <- terra::project(ol, "EPSG:4326", method = "near")
        } else if (inherits(ol, "sf")) {
          ol <- sf::st_transform(ol, 4326)
        } else if (is.list(ol) && all(c("raster","shp") %in% names(ol))) {
          ol$raster <- terra::project(ol$raster, "EPSG:4326", method = "near")
          ol$shp    <- sf::st_transform(ol$shp, 4326)
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
      r2 <- if (!identical(terra::crs(r), terra::crs(template))) {
        terra::project(r, template, method = "near")
      } else r
      terra::resample(r2, template, method = "near")
    }

    template <- .make_template(area_metric, crs_metric_wkt, pixel_size)

    # regrid distributions / kernels
    for (nm in names(spp_distribution_list)) {
      spp_distribution_list[[nm]]$raster <- .regrid_to_template(spp_distribution_list[[nm]]$raster, template)
    }
    for (nm in names(str_distribution_list)) {
      str_distribution_list[[nm]]$raster <- .regrid_to_template(str_distribution_list[[nm]]$raster, template)
    }
    for (nm in names(spp_kernel_list)) {
      spp_kernel_list[[nm]]$raster <- .regrid_to_template(spp_kernel_list[[nm]]$raster, template)
    }
    for (nm in names(stressor_kernel_list)) {
      stressor_kernel_list[[nm]]$raster <- .regrid_to_template(stressor_kernel_list[[nm]]$raster, template)
    }
    # overlaps (handle both SpatRaster and list with $raster)
    for (i in names(overlap_maps_list)) for (j in names(overlap_maps_list[[i]])) {
      ol <- overlap_maps_list[[i]][[j]]
      if (inherits(ol, "SpatRaster")) {
        overlap_maps_list[[i]][[j]] <- .regrid_to_template(ol, template)
      } else if (is.list(ol) && "raster" %in% names(ol)) {
        overlap_maps_list[[i]][[j]]$raster <- .regrid_to_template(ol$raster, template)
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
      spp_distribution_list[[nm]]$raster <- .set_crs_if_missing(spp_distribution_list[[nm]]$raster)
    }
    for (nm in names(str_distribution_list)) {
      str_distribution_list[[nm]]$raster <- .set_crs_if_missing(str_distribution_list[[nm]]$raster)
    }
    for (nm in names(spp_kernel_list)) {
      spp_kernel_list[[nm]]$raster <- .set_crs_if_missing(spp_kernel_list[[nm]]$raster)
    }
    for (nm in names(stressor_kernel_list)) {
      stressor_kernel_list[[nm]]$raster <- .set_crs_if_missing(stressor_kernel_list[[nm]]$raster)
    }
    for (i in names(overlap_maps_list)) for (j in names(overlap_maps_list[[i]])) {
      ol <- overlap_maps_list[[i]][[j]]
      if (inherits(ol, "SpatRaster")) {
        overlap_maps_list[[i]][[j]] <- .set_crs_if_missing(ol)
      } else if (is.list(ol) && "raster" %in% names(ol)) {
        overlap_maps_list[[i]][[j]]$raster <- .set_crs_if_missing(ol$raster)
      }
    }

    # Enforce square + metric (final guard)
    for (nm in names(spp_distribution_list)) .check_square_metric(spp_distribution_list[[nm]]$raster)
    for (nm in names(str_distribution_list)) .check_square_metric(str_distribution_list[[nm]]$raster)
    for (nm in names(spp_kernel_list)) .check_square_metric(spp_kernel_list[[nm]]$raster)
    for (nm in names(stressor_kernel_list)) .check_square_metric(stressor_kernel_list[[nm]]$raster)
    for (i in names(overlap_maps_list)) for (j in names(overlap_maps_list[[i]])) {
      ol <- overlap_maps_list[[i]][[j]]
      if (inherits(ol, "SpatRaster")) .check_square_metric(ol)
      if (is.list(ol) && "raster" %in% names(ol)) .check_square_metric(ol$raster)
    }

    # If overlaps must include rasters for HRA, enforce:
    for (i in names(overlap_maps_list)) for (j in names(overlap_maps_list[[i]])) {
      ol <- overlap_maps_list[[i]][[j]]
      has_raster <- inherits(ol, "SpatRaster") || (is.list(ol) && "raster" %in% names(ol))
      if (!has_raster && output_layer_type == "raster") {
        stop("Overlap map missing raster output; set output_layer_type to include raster.")
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
