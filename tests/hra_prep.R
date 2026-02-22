library(sf)


byra_prep <- function(spp_maplist, str_maplist,
                     area = NULL,
                     pixel_size = NULL,
                     dimyx = c(512,512),
                     reclass = c(1, 3),
                     reclass_cat = FALSE,
                     overlaps = TRUE,
                     crs = NULL,
                     quiet = TRUE) {

  spp_list <- spp_maplist
  str_list <- str_maplist

  # Helpers
  # Checks if all rasters in a list have the same crs
  .same_crs <- function(obj_list) {
    # empty or single-element list: TRUE
    if (length(obj_list) < 2) return(TRUE)

    # get a comparable CRS representation (WKT) for each object
    crs_wkt <- lapply(obj_list, function(x) {
      if (inherits(x, c("SpatRaster", "SpatVector"))) {
        terra::crs(x, proj = TRUE) # WKT string
      } else if (inherits(x, c("sf", "sfc"))) {
        sf::st_crs(x)$wkt # WKT string
      } else {
        stop("Unsupported object of class: ", paste(class(x), collapse = ", "),
             ". Expected 'SpatRaster', 'SpatVector', 'sf' or 'sfc'.")
      }
    })

    ref <- crs_wkt[[1]]

    # all CRSs must be identical to the first one
    all(vapply(crs_wkt, function(z) identical(z, ref), logical(1)))
  }

  # Checks if input is a SpatRaster or sf object.
  .check_sf_spatrast <- function(x) {
    allowed <- c("SpatRaster", "sf", "sfc", "SpatVector")
    if (!inherits(x, allowed)) {
      if (is.list(x)) {
        for (item in x) {
          if (!inherits(item, allowed)) {
            stop("Inputs must be or contain 'sf', 'sfc', 'SpatRaster' or 'SpatVector' objects.")
          }
        }
      } else {
        stop("Inputs must be 'sf', 'sfc', 'SpatRaster' or 'SpatVector' objects.")
      }
    }
  }


  # Basic checks
  # Check area
  if (is.null(area)) {
    if (!quiet) message("Area is NULL.")
  } else if (inherits(area, "bbox")) {
    area <- sf::st_as_sfc(area)
  } else if (is.data.frame(area)) {
    area <- df_to_shp(area)
  } else if (!inherits(area, "sf")) {
    stop("`area` must be NULL, bbox, data.frame, or sf.")
  }

  # Check if inputs are 'sf' or 'SpatRaster'
  .check_sf_spatrast(c(spp_list, str_list))

  # Check if all rasters and sf have the same crs: if not, transform!
  template <- NULL
  if (!.same_crs(c(spp_list, str_list))){
    if (!quiet) message("Input rasters do not have the same coordinate system.")
    if (is.null(area)) {
      if (!quiet) message("No area was provided. Taking the CRS of the first map from spp_list as template")
      if (is.list(spp_list)) {
        template <- spp_list[[1]]
      } else {
        template <- spp_list
      }
      if (!quiet) message("Reprojecting rasters to template's CRS...")
      spp_list <- lapply(spp_list,
                         function(x) align_to(x, template = template,
                                              categorical = terra::is.int(x) || terra::is.factor(x)))
      str_list <- lapply(str_list,
                         function(y) align_to(y, template = template,
                                              categorical = terra::is.int(y) || terra::is.factor(y)))
    } else {
      if (! quiet) message("Using area's crs as template: reprojecting rasters...")
      longlat <- sf::st_is_longlat(area)
      if (longlat) {
        area <- transform_to_metric(area)
      }
      template <- area
      spp_list <- lapply(spp_list,
                         function(x) align_to(x, template = template,
                                              categorical = terra::is.int(x)))
      str_list <- lapply(str_list,
                         function(y) align_to(y, template = template,
                                              categorical = terra::is.int(y)))
    }
  }

  # Creating area of interest, if not provided
  if (is.null(template)) {
    if (is.list(spp_list)) {
      template <- spp_list[[1]]
    } else {
      template <- spp_list
    }
  }
  if (is.null(area)) {
    if (!quiet) message("No area provided; creating AOI (bounding box) from all layers.")
    e <- NULL
    if (all(is.list(spp_list), is.list(str_list))) {
      e <- do.call(terra::ext, c(spp_list, str_list))
    } else if (is.list(spp_list) && !is.list(str_list)) {
      e <- do.call(terra::ext, c(spp_list, list(str_list)))
    } else if (!is.list(spp_list) && is.list(str_list)) {
      e <- do.call(terra::ext, c(list(spp_list), str_list))
    } else {
      e <- do.call(terra::ext, list(spp_list, str_list))
    }

    # Estimate buffer as 5% of the diagonal of the combined extent
    dx <- terra::xmax(e) - terra::xmin(e)
    dy <- terra::ymax(e) - terra::ymin(e)
    diag_len <- sqrt(dx^2 + dy^2)
    buf <- 0.05 * diag_len

    crs_template <- if (inherits(template, "SpatRaster")) {
      sf::st_crs(terra::crs(template, proj = TRUE))
    } else {
      sf::st_crs(template)
    }

    e_buf <- terra::ext(
      terra::xmin(e) - buf,
      terra::xmax(e) + buf,
      terra::ymin(e) - buf,
      terra::ymax(e) + buf)

    area <- sf::st_bbox(c(
      xmin = terra::xmin(e_buf),
      ymin = terra::ymin(e_buf),
      xmax = terra::xmax(e_buf),
      ymax = terra::ymax(e_buf)), crs = crs_template)

    area <- sf::st_as_sfc(area)
  }

  # If reclass is a vector, then apply a transformation
  spp_output <- apply_scale(spp_list, reclass, cat = reclass_cat)
  str_output <- apply_scale(str_list, reclass, cat = reclass_cat)

  # Estimate overlaps
  overlap_maps <- list()
  sp_names <- names(spp_output)
  st_names <- names(str_output)

  if (overlaps) {
    if (reclass_cat) {
      for (sp in sp_names) {
        for (st in st_names) {
          overlap_maps[[sp]][[st]] <- get_overlap_kernel(spp_output[[sp]],
                                                     str_output[[st]],
                                                     n_classes = max(reclass),
                                                     output_layer_type = "raster")
        }
      }
    } else {
      for (sp in sp_names) {
        for (st in st_names) {
          product_map <- spp_output[[sp]] * str_output[[st]]
          overlap_maps[[sp]][[st]] <- scale_vars(product_map, reclass)
        }
      }
    }
  }

  # Preparing outputs
  spp_distribution_list <- list()
  str_distribution_list <- list()

  for (sp in sp_names) {
    spp_distribution_list[[sp]] <- terra::ifel(spp_list[[sp]] > 0, 1, NA)
  }

  for (st in st_names) {
    str_distribution_list[[st]] <- terra::ifel(str_list[[st]] > 0, 1, NA)
  }

  # Output
  output <- list(
    species_distributions = spp_distribution_list,
    stressor_distributions = str_distribution_list,
    species_kernel_maps = spp_output,
    stressor_kernel_maps = str_output,
    overlap_maps = overlap_maps,
    area_of_interest = area)

  class(output) <- c("risaMaps", class(output))

  return(output)
}
