#' Prepare Species–Stressor Spatial Input Maps for risa
#'
#' Preprocesses lists of species and stressor spatial layers (SpatRaster, sf/sfc,
#' or SpatVector objects) for use in *risa* analyses.
#' The function harmonizes coordinate reference systems (CRS), generates an
#' area of interest (AOI) when none is supplied, scales/reclassifies variables,
#' computes overlap maps, and outputs species/stressor distribution and kernel
#' maps in a common projected CRS.
#'
#' @param spp_maplist A list of species layers. Each element must be a
#'   `SpatRaster`, `sf`, `sfc`, or `SpatVector` object. Names are auto-generated
#'   if missing.
#' @param str_maplist A list of stressor layers. Same requirements as
#'   `spp_maplist`. Names are auto-generated if missing.
#' @param area Optional area of interest. May be:
#'   - `NULL` (default) → AOI is built automatically from bounding boxes of
#'     species and stressor data.
#'   - A `bbox`,
#'   - A `data.frame` convertible to a polygon (via `df_to_shp()`),
#'   - An `sf`/`sfc` object.
#' @param pixel_size Not currently used internally. Reserved for future control
#'   of raster resolution.
#' @param dimyx A length-2 numeric vector giving output raster dimensions
#'   (rows, columns). Currently passed only implicitly to downstream helpers.
#' @param reclass Numeric vector used by `apply_scale()` to reclassify or scale
#'   species and stressor rasters. Default: `c(1, 3)`.
#' @param reclass_cat Logical. If `TRUE`, variables are treated as categorical
#'   and overlap is computed using class-based kernels. If `FALSE`, overlap is
#'   based on scaled continuous products. Default: `FALSE`.
#' @param overlaps Logical. Whether to compute species–stressor overlap maps.
#'   Default: `TRUE`.
#' @param crs Optional target CRS for outputs. Must be interpretable by
#'   `sf::st_crs()` and **must be projected (metric)**. If `NULL`, a suitable
#'   projected CRS (e.g. UTM) is auto-selected when inputs are in long/lat.
#' @param quiet Logical. If `FALSE`, progress messages are printed. Default:
#'   `TRUE`.
#'
#' @details
#' The function:
#'
#' 1. Validates that input layers are spatial (`SpatRaster`, `sf`, etc.).
#' 2. Ensures all inputs share a CRS; if not, they are reprojected to the CRS of
#'    the AOI or the first species map.
#' 3. If `area` is not provided, constructs an AOI by merging and buffering the
#'    extents of all rasters.
#' 4. Applies scaling/reclassification to species and stressor layers.
#' 5. Optionally computes species–stressor overlap maps with continuous or
#'    categorical kernels.
#' 6. Derives or accepts a target **metric CRS** and reprojects all outputs.
#'
#' Helper functions used internally include:
#' - `apply_scale()`
#' - `scale_vars()`
#' - `get_overlap_kernel()`
#' - `align_to()`
#' - `transform_to_metric()`
#' - `df_to_shp()`
#'
#' These must be available in the namespace where `byra_prep()` is used.
#'
#' @return
#' A list of class `"risaMaps"` with components:
#'
#' \describe{
#'   \item{`species_distributions`}{List of binary species presence/absence rasters.}
#'   \item{`stressor_distributions`}{List of binary stressor presence/absence rasters.}
#'   \item{`species_kernel_maps`}{List of scaled/reclassified species rasters.}
#'   \item{`stressor_kernel_maps`}{List of scaled/reclassified stressor rasters.}
#'   \item{`overlap_maps`}{Nested list of species–stressor overlap rasters
#'     (empty if `overlaps = FALSE`).}
#'   \item{`area_of_interest`}{`sf` polygon defining the spatial AOI in the final
#'     projected CRS.}
#' }
#'
#' @seealso
#' `apply_scale()`, `get_overlap_kernel()`, `align_to()`,
#' `transform_to_metric()`
#'
#' @import terra
#' @import sf
#'
#' @export
byra_prep <- function(spp_maplist, str_maplist,
                      area = NULL,
                      pixel_size = NULL,
                      dimyx = c(512,512),
                      reclass = c(1, 3),
                      reclass_cat = FALSE,
                      overlaps = TRUE,
                      crs = NULL,
                      quiet = TRUE) {

  # Preparing input
  spp_list <- spp_maplist
  str_list <- str_maplist

  if (is.null(names(spp_list))) {
    names(spp_list) <- c(paste("species_", 1:length(spp_list), sep=""))
  }
  if (is.null(names(str_list))) {
    names(str_list) <- c(paste("stressor_", 1:length(str_list), sep=""))
  }

  # Helpers
  # Small helper: x %||% y
  `%||%` <- function(x, y) {
    if (is.null(x) || (is.atomic(x) && length(x) == 1L && is.na(x))) y else x
  }

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

  # Put a single SpatRaster into a list!
  .as_rast_list <- function(x) {
    if (is.null(x)) {
      list()
    } else if (inherits(x, "SpatRaster")) {
      list(x)
    } else if (is.list(x) && all(vapply(x, inherits, logical(1), "SpatRaster"))) {
      x
    } else {
      stop("Input must be a SpatRaster or a list of SpatRaster objects (or NULL).")
    }
  }

  # Project a list of SpatRasters to a target CRS (terra-style string)
  .project_raster_list <- function(lst, target_crs) {
    lapply(lst, function(r) {
      if (is.null(r)) return(NULL)
      if (!inherits(r, "SpatRaster")) return(r)
      current <- terra::crs(r, proj = TRUE)
      if (!is.na(current) && current == target_crs) return(r)
      method <- if (terra::is.int(r) || terra::is.factor(r)) "near" else "bilinear"
      terra::project(r, target_crs, method = method)
    })
  }

  # Determine a representative object to infer current CRS
  .get_reference_obj <- function(area, spp_list) {
    if (!is.null(area)) return(area)
    if (is.list(spp_list)) return(spp_list[[1]])
    spp_list
  }

  # Basic checks
  # Check area
  if (is.null(area)) {
    if (!quiet) message("Area is NULL.")
  } else if (inherits(area, "bbox")) {
    area <- sf::st_as_sfc(area)
  } else if (is.data.frame(area)) {
    area <- df_to_shp(area)
  } else if (!inherits(area, c("sf", "sfc"))) {
    stop("`area` must be NULL, bbox, data.frame, or sf.")
  }

  # Check if inputs are 'sf' or 'SpatRaster'
  .check_sf_spatrast(c(spp_list, str_list))

  # Harmonize CRS among inputs
  template <- NULL
  if (!.same_crs(c(spp_list, str_list))) {
    if (!quiet) message("Input rasters do not have the same coordinate system.")
    if (is.null(area)) {
      if (!quiet) message("No area was provided. Taking the CRS of the first map from spp_list as template")
      if (is.list(spp_list)) {
        template <- spp_list[[1]]
      } else {
        template <- spp_list
      }
      if (!quiet) message("Reprojecting rasters to template's CRS...")

      # When area is NULL and we reproject to template CRS
      spp_list <- lapply(
        spp_list,
        function(x) {
          categorical_flag <- if (inherits(x, "SpatRaster")) {
            terra::is.int(x) || terra::is.factor(x)
          } else {
            TRUE
          }
          align_to(x, template = template, categorical = categorical_flag)
        }
      )

      str_list <- lapply(
        str_list,
        function(y) {
          categorical_flag <- if (inherits(y, "SpatRaster")) {
            terra::is.int(y) || terra::is.factor(y)
          } else {
            TRUE
          }
          align_to(y, template = template, categorical = categorical_flag)
        }
      )


    } else {
      if (!quiet) message("Using area's crs as template: reprojecting rasters...")
      longlat <- sf::st_is_longlat(area)
      if (longlat) {
        area_m <- transform_to_metric(area, quiet = quiet)
        area <- area_m$shape
      }
      template <- area

      # When area is provided and used as template
      spp_list <- lapply(
        spp_list,
        function(x) {
          categorical_flag <- if (inherits(x, "SpatRaster")) {
            terra::is.int(x) || terra::is.factor(x)
          } else {
            TRUE
          }
          align_to(x, template = template, categorical = categorical_flag)
        }
      )

      str_list <- lapply(
        str_list,
        function(y) {
          categorical_flag <- if (inherits(y, "SpatRaster")) {
            terra::is.int(y) || terra::is.factor(y)
          } else {
            TRUE
          }
          align_to(y, template = template, categorical = categorical_flag)
        }
      )

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
    if (!quiet) message("No area provided. Creating AOI (bounding box) from all layers.")
    # Flatten to a single list of SpatRaster objects
    rasters_to_unite <- c(
      .as_rast_list(spp_list),
      .as_rast_list(str_list)
    )

    # Get extents for each raster
    exts <- lapply(
      rasters_to_unite,
      function(r) {
        r_trim <- terra::trim(r) # drops NA-only borders
        terra::ext(r_trim) # extent of non-NA area
      }
    )

    # Compute overall bounding box
    xmin_vals <- vapply(exts, terra::xmin, numeric(1))
    xmax_vals <- vapply(exts, terra::xmax, numeric(1))
    ymin_vals <- vapply(exts, terra::ymin, numeric(1))
    ymax_vals <- vapply(exts, terra::ymax, numeric(1))

    e <- terra::ext(
      min(xmin_vals),
      max(xmax_vals),
      min(ymin_vals),
      max(ymax_vals)
    )

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
      terra::ymax(e) + buf
    )

    area_bbox <- sf::st_bbox(
      c(
        xmin = terra::xmin(e_buf),
        ymin = terra::ymin(e_buf),
        xmax = terra::xmax(e_buf),
        ymax = terra::ymax(e_buf)
      ),
      crs = crs_template
    )

    area <- sf::st_as_sf(sf::st_as_sfc(area_bbox))
  }

  # Kernel / scaling
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
          overlap_maps[[sp]][[st]] <- get_overlap_kernel(
            spp_output[[sp]],
            str_output[[st]],
            n_classes = max(reclass),
            output_layer_type = "raster"
          )
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

  # Distribution (binary) maps
  spp_distribution_list <- list()
  str_distribution_list <- list()

  for (sp in sp_names) {
    spp_distribution_list[[sp]] <- terra::ifel(spp_list[[sp]] > 0, 1, NA)
  }

  for (st in st_names) {
    str_distribution_list[[st]] <- terra::ifel(str_list[[st]] > 0, 1, NA)
  }

  # Force all outputs to a common metric CRS
  # Determine target metric CRS as sf::crs object
  if (!is.null(crs)) {
    # User-specified CRS
    target_crs_sf <- sf::st_crs(crs)
    if (is.na(target_crs_sf)) {
      stop("Could not interpret 'crs' argument. Provide something understood by sf::st_crs(), e.g. 'EPSG:32722'.")
    }
    # Check that it is projected (not geographic)
    dummy <- sf::st_sfc(sf::st_point(c(0, 0)), crs = target_crs_sf)
    if (sf::st_is_longlat(dummy)) {
      stop("The provided 'crs' is geographic (lon/lat). Please provide a projected (metric) CRS in 'crs'.")
    }
  } else {
    # Auto-select metric CRS if current is geographic
    ref_obj <- .get_reference_obj(area, spp_list)

    ref_is_geo <- if (inherits(ref_obj, c("sf", "sfc"))) {
      sf::st_is_longlat(ref_obj)
    } else if (inherits(ref_obj, "SpatRaster")) {
      terra::is.lonlat(ref_obj)
    } else {
      FALSE
    }

    if (ref_is_geo) {
      # Use transform_to_metric() to pick a suitable metric CRS
      if (inherits(ref_obj, c("sf", "sfc"))) {
        info <- transform_to_metric(ref_obj, quiet = quiet)
        target_crs_sf <- sf::st_crs(info$crs)
        # also update area to metric if area was the reference
        if (!is.null(area) && identical(ref_obj, area)) {
          area <- info$shape
        }
      } else {
        # SpatRaster case
        rst_m <- transform_to_metric(ref_obj, quiet = quiet)
        target_crs_sf <- sf::st_crs(terra::crs(rst_m, proj = TRUE))
      }
      if (!quiet) {
        message("Auto-selected metric CRS (UTM/polar) for outputs: EPSG:",
                target_crs_sf$epsg %||% NA_integer_)
      }
    } else {
      # Already projected: just reuse that CRS
      if (inherits(ref_obj, c("sf", "sfc"))) {
        target_crs_sf <- sf::st_crs(ref_obj)
      } else if (inherits(ref_obj, "SpatRaster")) {
        target_crs_sf <- sf::st_crs(terra::crs(ref_obj, proj = TRUE))
      } else {
        stop("Could not determine a reference CRS for auto metric selection.")
      }
    }
  }

  # Build a terra-compatible CRS string
  target_crs_terra <- target_crs_sf$wkt %||%
    target_crs_sf$proj4string %||%
    if (!is.null(target_crs_sf$epsg)) paste0("EPSG:", target_crs_sf$epsg) else NA_character_

  if (is.na(target_crs_terra)) {
    stop("Could not derive a valid target CRS string from the chosen metric CRS.")
  }

  # Reproject raster outputs to target metric CRS
  spp_distribution_list <- .project_raster_list(spp_distribution_list, target_crs_terra)
  str_distribution_list <- .project_raster_list(str_distribution_list, target_crs_terra)
  spp_output <- .project_raster_list(spp_output, target_crs_terra)
  str_output <- .project_raster_list(str_output, target_crs_terra)

  if (length(overlap_maps) > 0) {
    overlap_maps <- lapply(
      overlap_maps,
      function(sp_list) {
        lapply(sp_list, function(r) {
          if (!inherits(r, "SpatRaster")) return(r)
          current <- terra::crs(r, proj = TRUE)
          if (!is.na(current) && current == target_crs_terra) return(r)
          method <- if (terra::is.int(r) || terra::is.factor(r)) "near" else "bilinear"
          terra::project(r, target_crs_terra, method = method)
        })
      }
    )
  }

  # Reproject area_of_interest to target metric CRS

  if (inherits(area, c("sf", "sfc"))) {
    area_crs <- sf::st_crs(area)
    # Compare via EPSG or WKT
    if (is.na(area_crs$epsg) || is.na(target_crs_sf$epsg) ||
        area_crs$epsg != target_crs_sf$epsg) {
      area <- sf::st_transform(area, target_crs_sf)
    }
  }

  # Output
  output <- list(
    species_distributions = spp_distribution_list,
    stressor_distributions = str_distribution_list,
    species_kernel_maps = spp_output,
    stressor_kernel_maps = str_output,
    overlap_maps = overlap_maps,
    area_of_interest = area
  )

  class(output) <- c("risaMaps", class(output))

  return(output)
}

