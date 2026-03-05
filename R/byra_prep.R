#' harmonizes input maps into a common CRS, extent, origin, and resolution
#'
#' @description
#' `byra_prep()` harmonizes (aligns) input species and stressor layers so they share a
#' common coordinate reference system (CRS), extent, origin and resolution (pixel size),
#' converting vector inputs (`sf`/`sfc`/`SpatVector`) to rasters on a single template grid.
#' It then applies a kernel/scaling step (via `apply_scale()`), optionally computes
#' pairwise overlap maps, and returns a structured object of class `risaMaps` ready to be
#' used as input by the function `quick_byra()`.
#'
#' @details
#' The function accepts lists of maps for species and stressors. Each element may be a
#' `SpatRaster`, `sf`, `sfc`, or `SpatVector`. If an AOI is provided
#' via `area`, it is used as the template domain; otherwise the template is inferred from
#' the first species layer and, later, an AOI (bounding box) is derived from the union of
#' non-`NA` extents of all layers.
#'
#' When `area` (or the inferred template) is geographic (lon/lat), the function attempts
#' to transform it to a projected (metric) CRS (via `transform_to_metric()`) before
#' building the template grid. If `pixel_size` is not supplied and the template is a
#' vector object, the function tries to infer a suitable pixel size from the resolution
#' of any raster inputs in a matching CRS; otherwise it approximates pixel size from the
#' AOI bounding box and the requested grid dimensions (`dimyx`).
#'
#' Vector layers are rasterized onto the template using `terra::rasterize()` (with
#' `touches = TRUE`), while rasters are reprojected/resampled to the template with
#' `terra::project()` using nearest-neighbour for categorical/integer data and bilinear
#' for continuous data. After kernels and overlaps are computed, outputs are optionally
#' forced to a user-supplied projected CRS (`crs`) while preserving a consistent grid.
#'
#' @param spp_maplist A list of species maps. Each element must be a
#'   `terra::SpatRaster`, `sf::sf`, `sf::sfc`, or `terra::SpatVector`.
#'   If unnamed, default names `"species_1"`, `"species_2"`, ... are assigned.
#'
#' @param str_maplist A list of stressor maps. Each element must be a
#'   `terra::SpatRaster`, `sf::sf`, `sf::sfc`, or `terra::SpatVector`.
#'   If unnamed, default names `"stressor_1"`, `"stressor_2"`, ... are assigned.
#'
#' @param area Optional area of interest (AOI). Can be `NULL`, an `sf::bbox`,
#'   a `data.frame` (converted via `df_to_shp()`), or an `sf`/`sfc` object.
#'   If `NULL`, an AOI is created as a buffered bounding box encompassing all layers
#'   after alignment.
#'
#' @param pixel_size Optional numeric pixel size (in template units). If `area`/template
#'   is vector-based and `pixel_size` is `NULL`/`NA`, the function will attempt to infer
#'   a suitable value from raster inputs or from `dimyx`.
#'
#' @param dimyx Integer vector of length 2 giving the target grid dimensions
#'   `(nrow, ncol)` used only when `pixel_size` must be inferred from a vector template.
#'
#' @param reclass Integer vector of length 2 giving the minimum and maximum class values
#'   used by `apply_scale()`/`scale_vars()` to reclassify kernel maps (e.g. `c(1, 3)`).
#'
#' @param reclass_cat Logical; if `TRUE`, treats kernel layers as categorical and
#'   computes overlaps using `get_overlap_kernel()` with `n_classes = max(reclass)`.
#'   If `FALSE` (default), overlaps are computed as the product of scaled rasters and
#'   then re-binned with `scale_vars()`.
#'
#' @param overlaps Logical; if `TRUE` (default) compute pairwise overlap maps for all
#'   species × stressor combinations.
#'
#' @param crs Optional target CRS for outputs. Anything understood by `sf::st_crs()`
#'   (e.g. `"EPSG:32722"`). Must be projected (metric), not lon/lat. If `NULL`, a metric
#'   CRS is auto-selected when inputs are geographic; otherwise the template CRS is used.
#'
#' @param quiet Logical; if `FALSE`, prints progress messages and basic geometry checks.
#'
#' @return
#' An object of class `risaMaps` (a named list) with components:
#' \describe{
#'   \item{species_distributions}{List of binary (`1`/`NA`) species distribution rasters.}
#'   \item{stressor_distributions}{List of binary (`1`/`NA`) stressor distribution rasters.}
#'   \item{species_kernel_maps}{List of scaled species kernel maps (rasters).}
#'   \item{stressor_kernel_maps}{List of scaled stressor kernel maps (rasters).}
#'   \item{overlap_maps}{Nested list of overlap maps `[species][stressor]` (possibly empty).}
#'   \item{area_of_interest}{AOI as an `sf` object in the target CRS.}
#' }
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(sf)
#'
#' # Example inputs (rasters or vectors)
#' spp <- list(sp1 = rast("sp1.tif"), sp2 = rast("sp2.tif"))
#' str <- list(fishery = rast("trawl.tif"))
#'
#' # Optional AOI
#' aoi <- st_read("aoi.shp")
#'
#' out <- byra_prep(
#'   spp_maplist = spp,
#'   str_maplist = str,
#'   area = aoi,
#'   pixel_size = 100,
#'   reclass = c(1, 3),
#'   overlaps = TRUE,
#'   quiet = FALSE
#' )
#'
#' class(out)
#' names(out)
#' }
#'
#' @export
byra_prep <- function(spp_maplist, str_maplist,
                      area = NULL,
                      pixel_size = NULL,
                      dimyx = c(512, 512),
                      reclass = c(1, 3),
                      reclass_cat = FALSE,
                      overlaps = TRUE,
                      crs = NULL,
                      quiet = TRUE) {

  # Preparing input
  spp_list <- spp_maplist
  str_list <- str_maplist

  if (is.null(names(spp_list))) names(spp_list) <- paste0("species_", seq_along(spp_list))
  if (is.null(names(str_list))) names(str_list) <- paste0("stressor_", seq_along(str_list))

  # Helpers
  # Small helper opperator
  `%||%` <- function(x, y) {
    if (is.null(x) || (is.atomic(x) && length(x) == 1L && is.na(x))) y else x
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

  # Checks if all rasters in a list have the same crs
  .same_crs <- function(obj_list) {
    if (length(obj_list) < 2) return(TRUE)

    crs_wkt <- lapply(obj_list, function(x) {
      if (inherits(x, c("SpatRaster", "SpatVector"))) {
        terra::crs(x, proj = TRUE)
      } else if (inherits(x, c("sf", "sfc"))) {
        sf::st_crs(x)$wkt
      } else {
        stop("Unsupported class: ", paste(class(x), collapse = ", "))
      }
    })

    ref <- crs_wkt[[1]]
    all(vapply(crs_wkt, function(z) identical(z, ref), logical(1)))
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

  # Determine a representative object to infer current CRS
  .get_reference_obj <- function(area, spp_list) {
    if (!is.null(area)) return(area)
    if (is.list(spp_list)) return(spp_list[[1]])
    spp_list
  }

  # Check if input has categorical (or class/code) values
  .is_categorical <- function(x) {
    if (!inherits(x, "SpatRaster")) return(TRUE)
    terra::is.int(x) || terra::is.factor(x)
  }

  # choose a default field in sf/spatvector if none provided
  .guess_field <- function(x) {
    if (inherits(x, "sf")) {
      nms <- names(x)
      nms <- nms[nms != attr(sf::st_geometry(x), "sf_column")]
      if (length(nms) == 0) return(NULL)
      return(nms[1])
    }
    if (inherits(x, "SpatVector")) {
      nms <- names(x)
      if (length(nms) == 0) return(NULL)
      return(nms[1])
    }
    NULL
  }

  # Build a single raster template from an sf/sfc bbox + pixel_size
  .make_template_raster_from_sf <- function(sf_obj, pixel_size) {
    stopifnot(inherits(sf_obj, c("sf", "sfc")))
    if (!is.finite(pixel_size) || is.na(pixel_size) || pixel_size <= 0) {
      stop("pixel_size must be a positive finite number.")
    }
    bb <- sf::st_bbox(sf_obj)
    tpl <- terra::rast(
      terra::ext(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]),
      resolution = pixel_size,
      crs = sf::st_crs(sf_obj)$wkt
    )
    tpl <- terra::extend(tpl, tpl)
    tpl
  }

  # rasterize a vector object onto a template raster grid
  .rasterize_to_template <- function(v, template_raster, field = NULL,
                                     background = NA, touches = TRUE) {
    stopifnot(inherits(template_raster, "SpatRaster"))

    if (inherits(v, "sf")) {
      v <- terra::vect(v)
    } else if (!inherits(v, "SpatVector")) {
      stop("Expected sf or SpatVector in .rasterize_to_template().")
    }

    if (is.null(field)) field <- .guess_field(v)

    # If no attribute field exists, use '1' values where feature is present
    if (is.null(field)) {
      r <- terra::rasterize(v, template_raster, field = 1, background = background, touches = touches)
    } else {
      r <- terra::rasterize(v, template_raster, field = field, background = background, touches = touches)
    }

    # lock extent exactly
    r <- terra::crop(r, template_raster, snap = "out")
    terra::ext(r) <- terra::ext(template_raster)
    r
  }

  # Align any object (sf, sfc, SpatRaster) to a SpatRaster template and return a raster
  .align_to_template_raster <- function(x, template_raster,
                                        categorical = TRUE,
                                        vector_field = NULL) {
    stopifnot(inherits(template_raster, "SpatRaster"))

    # If rasters: project-to-template
    if (inherits(x, "SpatRaster")) {
      method <- if (categorical) "near" else "bilinear"
      out <- terra::project(x, template_raster, method = method)
      out <- terra::crop(out, template_raster, snap = "out")
      terra::ext(out) <- terra::ext(template_raster)
      return(out)
    }

    # If sf/sfc/spatvector: transform + rasterize-to-template
    if (inherits(x, c("sf", "sfc"))) {
      # ensure sf has CRS and transform
      x_sf <- if (inherits(x, "sfc")) sf::st_as_sf(x) else x
      target <- sf::st_crs(terra::crs(template_raster, proj = TRUE))
      x_sf <- sf::st_transform(x_sf, target)
      return(.rasterize_to_template(x_sf, template_raster, field = vector_field))
    }

    if (inherits(x, "SpatVector")) {
      x_v <- terra::project(x, terra::crs(template_raster, proj = TRUE))
      return(.rasterize_to_template(x_v, template_raster, field = vector_field))
    }

    stop("Unsupported input class: ", paste(class(x), collapse = ", "))
  }

  # Project a list of SpatRasters to a template raster
  .project_raster_list_to_template <- function(lst, template_raster) {
    lapply(lst, function(r) {
      if (is.null(r)) return(NULL)
      if (!inherits(r, "SpatRaster")) return(r)
      method <- if (terra::is.int(r) || terra::is.factor(r)) "near" else "bilinear"
      out <- terra::project(r, template_raster, method = method)
      out <- terra::crop(out, template_raster, snap = "out")
      terra::ext(out) <- terra::ext(template_raster)
      out
    })
  }

  # Basic checks and area handling
  if (is.null(area)) {
    if (!quiet) message("Area is NULL.")
  } else if (inherits(area, "bbox")) {
    area <- sf::st_as_sfc(area)
  } else if (is.data.frame(area)) {
    area <- df_to_shp(area)
  } else if (!inherits(area, c("sf", "sfc"))) {
    stop("`area` must be NULL, bbox, data.frame, or sf.")
  }

  .check_sf_spatrast(c(spp_list, str_list))

  # Harmonize CRS + build ONE template grid
  template <- NULL
  template_raster <- NULL

  if (!.same_crs(c(spp_list, str_list))) {
    if (!quiet) message("Input maps do not have the same coordinate system.")
  }

  if (!is.null(area)) {
    if (sf::st_is_longlat(area)) {
      if (!quiet) message("Area is lon/lat: transforming to metric CRS...")
      area_m <- transform_to_metric(area, quiet = quiet)
      area <- area_m$shape
    }
    template <- area
  } else {
    template <- if (is.list(spp_list)) spp_list[[1]] else spp_list
  }

  # pixel_size inference if template is sf/sfc and pixel_size missing
  if (inherits(template, c("sf", "sfc")) &&
      (is.null(pixel_size) || (is.atomic(pixel_size) && length(pixel_size) == 1 && is.na(pixel_size)))) {

    if (!quiet) message("Template is a sf object and pixel_size is not provided. Inferring pixel size...")

    tpl_wkt <- sf::st_crs(template)$wkt
    rasters_in <- c(spp_list, str_list)
    rasters_in <- rasters_in[vapply(rasters_in, inherits, logical(1), "SpatRaster")]

    pix <- NA_real_

    if (length(rasters_in) > 0) {
      crs_match <- vapply(rasters_in, function(r) terra::crs(r, proj = TRUE) == tpl_wkt, logical(1))
      not_lonlat <- vapply(rasters_in, function(r) !terra::is.lonlat(r), logical(1))

      idx <- which(crs_match & not_lonlat)
      if (length(idx) == 0) idx <- which(not_lonlat)
      if (length(idx) == 0) idx <- 1L

      r0 <- rasters_in[[idx[1]]]

      if (terra::crs(r0, proj = TRUE) != tpl_wkt) {
        method <- if (terra::is.int(r0) || terra::is.factor(r0)) "near" else "bilinear"
        r0 <- terra::project(r0, tpl_wkt, method = method)
      }

      rr <- terra::res(r0)
      pix <- suppressWarnings(as.numeric(rr[1]))
    }

    if (!is.finite(pix) || is.na(pix)) {
      bb <- sf::st_bbox(template)
      xres <- as.numeric((bb["xmax"] - bb["xmin"]) / dimyx[2])
      yres <- as.numeric((bb["ymax"] - bb["ymin"]) / dimyx[1])
      pix <- max(xres, yres)
    }

    pixel_size <- pix
    if (!quiet) message("Using pixel_size = ", format(pixel_size, digits = 12), " (template units).")
  }

  # Build the raster template grid
  if (inherits(template, c("sf", "sfc"))) {
    if (is.null(pixel_size) || !is.finite(pixel_size) || is.na(pixel_size)) {
      stop("Could not infer pixel_size for sf template. Please provide pixel_size explicitly.")
    }
    template_raster <- .make_template_raster_from_sf(template, pixel_size)
  } else if (inherits(template, "SpatRaster")) {
    template_raster <- template
  } else if (inherits(template, "SpatVector")) {
    if (is.null(pixel_size) || !is.finite(pixel_size) || is.na(pixel_size)) {
      stop("Template is SpatVector and pixel_size is not provided. Please provide pixel_size explicitly.")
    }
    ex <- terra::ext(template)
    template_raster <- terra::rast(ex, resolution = pixel_size, crs = terra::crs(template, proj = TRUE))
    template_raster <- terra::extend(template_raster, template_raster)
  } else {
    stop("Could not build a template raster from the provided inputs.")
  }

  # Align EVERYTHING to template_raster and convert vectors into rasters
  spp_list <- lapply(spp_list, function(x) {
    .align_to_template_raster(x, template_raster, categorical = .is_categorical(x))
  })
  str_list <- lapply(str_list, function(x) {
    .align_to_template_raster(x, template_raster, categorical = .is_categorical(x))
  })

  if (!quiet) {
    rasters_now <- c(spp_list, str_list)
    bad <- vapply(rasters_now, function(r) !terra::compareGeom(r, template_raster, stopOnError = FALSE), logical(1))
    if (any(bad)) {
      message("WARNING: Some rasters still do not match template geometry: ",
              paste(names(rasters_now)[bad], collapse = ", "))
    } else {
      message("All rasters match the template geometry (extent/res/origin).")
    }
  }

  # Create AOI (bounding box) if area not provided
  if (is.null(area)) {
    if (!quiet) message("No area provided. Creating AOI (bounding box) from all raster layers.")

    rasters_to_unite <- c(.as_rast_list(spp_list), .as_rast_list(str_list))
    exts <- lapply(rasters_to_unite, function(r) terra::ext(terra::trim(r)))

    xmin_vals <- vapply(exts, terra::xmin, numeric(1))
    xmax_vals <- vapply(exts, terra::xmax, numeric(1))
    ymin_vals <- vapply(exts, terra::ymin, numeric(1))
    ymax_vals <- vapply(exts, terra::ymax, numeric(1))

    e <- terra::ext(min(xmin_vals), max(xmax_vals), min(ymin_vals), max(ymax_vals))

    dx <- terra::xmax(e) - terra::xmin(e)
    dy <- terra::ymax(e) - terra::ymin(e)
    diag_len <- sqrt(dx^2 + dy^2)
    buf <- 0.05 * diag_len

    e_buf <- terra::ext(terra::xmin(e) - buf, terra::xmax(e) + buf,
                        terra::ymin(e) - buf, terra::ymax(e) + buf)

    crs_template <- sf::st_crs(terra::crs(template_raster, proj = TRUE))

    area_bbox <- sf::st_bbox(
      c(xmin = terra::xmin(e_buf), ymin = terra::ymin(e_buf),
        xmax = terra::xmax(e_buf), ymax = terra::ymax(e_buf)),
      crs = crs_template
    )

    area <- sf::st_as_sf(sf::st_as_sfc(area_bbox))
  }

  # Kernel / scaling
  spp_output <- apply_scale(spp_list, reclass, cat = reclass_cat)
  str_output <- apply_scale(str_list, reclass, cat = reclass_cat)

  # Overlaps
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

  # Force outputs to target metric CRS while preserving grid
  if (!is.null(crs)) {
    target_crs_sf <- sf::st_crs(crs)
    if (is.na(target_crs_sf)) {
      stop("Could not interpret 'crs'. Use something sf::st_crs() understands, e.g. 'EPSG:32722'.")
    }
    dummy <- sf::st_sfc(sf::st_point(c(0, 0)), crs = target_crs_sf)
    if (sf::st_is_longlat(dummy)) stop("Provided 'crs' is lon/lat. Please provide a projected (metric) CRS.")
  } else {
    ref_obj <- .get_reference_obj(area, spp_list)

    ref_is_geo <- if (inherits(ref_obj, c("sf", "sfc"))) {
      sf::st_is_longlat(ref_obj)
    } else if (inherits(ref_obj, "SpatRaster")) {
      terra::is.lonlat(ref_obj)
    } else if (inherits(ref_obj, "SpatVector")) {
      terra::is.lonlat(ref_obj)
    } else {
      FALSE
    }

    if (ref_is_geo) {
      info <- transform_to_metric(ref_obj, quiet = quiet)
      if (inherits(info, c("sf", "sfc"))) {
        target_crs_sf <- sf::st_crs(info)
      } else if (inherits(info, "SpatRaster")) {
        target_crs_sf <- sf::st_crs(terra::crs(info, proj = TRUE))
      } else {
        target_crs_sf <- sf::st_crs(info$crs %||% info$shape)
      }
      if (!quiet) message("Auto-selected metric CRS for outputs: EPSG:", target_crs_sf$epsg %||% NA_integer_)
    } else {
      # If already projected, reuse template raster CRS
      target_crs_sf <- sf::st_crs(terra::crs(template_raster, proj = TRUE))
    }
  }

  target_crs_terra <- target_crs_sf$wkt %||%
    target_crs_sf$proj4string %||%
    if (!is.null(target_crs_sf$epsg)) paste0("EPSG:", target_crs_sf$epsg) else NA_character_

  if (is.na(target_crs_terra)) stop("Could not derive a valid target CRS string.")

  # Build a target template raster (grid) in target CRS
  target_template_raster <- terra::project(template_raster, target_crs_terra, method = "near")
  target_template_raster <- terra::extend(target_template_raster, target_template_raster)

  spp_distribution_list <- .project_raster_list_to_template(spp_distribution_list, target_template_raster)
  str_distribution_list <- .project_raster_list_to_template(str_distribution_list, target_template_raster)
  spp_output <- .project_raster_list_to_template(spp_output, target_template_raster)
  str_output <- .project_raster_list_to_template(str_output, target_template_raster)

  if (length(overlap_maps) > 0) {
    overlap_maps <- lapply(overlap_maps, function(sp_list) {
      lapply(sp_list, function(r) {
        if (!inherits(r, "SpatRaster")) return(r)
        method <- if (terra::is.int(r) || terra::is.factor(r)) "near" else "bilinear"
        out <- terra::project(r, target_template_raster, method = method)
        out <- terra::crop(out, target_template_raster, snap = "out")
        terra::ext(out) <- terra::ext(target_template_raster)
        out
      })
    })
  }

  # Reproject area_of_interest to target CRS
  if (inherits(area, c("sf", "sfc"))) {
    if (sf::st_crs(area)$wkt != target_crs_sf$wkt) {
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
  output
}


