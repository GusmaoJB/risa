#' Prepare and generate maps for risk assessment analysis
#'
#' Generates kernel density maps for species and stressors, plus overlap maps,
#' within an area of interest (auto or user-supplied).
#'
#' @param x Species input: `sf`, data.frame, or a list of `sf`. If a data.frame/sf,
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
#' @importFrom sf st_as_sfc st_transform
#' @importFrom terra project
#' @examples
#' # Creating test data
#' spp_df <- rbind(data.frame(long = rnorm(80, 0, 10),
#' lat = rnorm(80, 0, 10), species = "sp1"),
#' data.frame(long = rnorm(60, 0, 10),
#' lat = rnorm(60, 0, 10), species = "sp2"))
#' str_df <- rbind(data.frame(long = rnorm(100, 0, 10),
#' lat = rnorm(100, 0, 10), stressor = "trawling"),
#' data.frame(long = rnorm(50, 0, 10),
#' lat = rnorm(100, 0, 5), stressor = "gillnet"))
#' # Create kernel maps of species and stressor distributions and overlap maps
#' risa_maps <- risa_prep(spp_df, str_df)
#'
#' # Species and Stressor distributions
#' dev.off()
#' par(mfrow = c(2,2))
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species 1")
#' plot(risa_maps$species_kernel_maps$sp1$shp, add = TRUE)
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species 2")
#' plot(risa_maps$species_kernel_maps$sp2$shp, add = TRUE)
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Gillnet")
#' plot(risa_maps$stressor_kernel_maps$gillnet$shp, add = TRUE, col=c("lightgreen", "green", "darkgreen"))
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Trawling")
#' plot(risa_maps$stressor_kernel_maps$trawling$shp, add = TRUE, col=c("lightgreen", "green", "darkgreen"))
#' # Overlap maps
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species1 vs Gillnet")
#' plot(risa_maps$overlap_maps$sp1$gillnet$shp, add = TRUE, col=c("yellow", "orange", "red"))
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species2 vs Gillnet")
#' plot(risa_maps$overlap_maps$sp2$gillnet$shp, add = TRUE, col=c("yellow", "orange", "red"))
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species1 vs Trawling")
#' plot(risa_maps$overlap_maps$sp1$trawling$shp, add = TRUE, col=c("yellow", "orange", "red"))
#' plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species2 vs Trawling")
#' plot(risa_maps$overlap_maps$sp2$trawling$shp, add = TRUE, col=c("yellow", "orange", "red"))
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
    quiet = TRUE
) {
  output_layer_type <- match.arg(output_layer_type)
  radius_method <- match.arg(radius_method)
  area_strategy <- match.arg(area_strategy)
  area_type <- match.arg(area_type)
  return_crs <- match.arg(return_crs)
  overlap_method <- match.arg(overlap_method)

  # Helper: normalize to a named list of sf
  as_sf_list <- function(obj, group = NULL, label_prefix = "layer") {
    # sf: wrap as single-element list
    if (inherits(obj, "sf")) {
      nm <- if (!is.null(names(obj)) && any(nzchar(names(obj)))) names(obj)[1] else label_prefix
      return(setNames(list(obj), nm))
    }

    # If data.frame: split or convert
    if (is.data.frame(obj)) {
      if (!is.null(group) && group %in% names(obj)) {
        # split by explicit group col
        sp <- split(obj, obj[[group]], drop = TRUE)
        lst <- lapply(sp, df_to_shp)
        return(lst)
      } else {
        return(df_to_list(obj))
      }
    }

    # If generic list: each element must be sf or data.frame
    if (is.list(obj)) {
      nms <- names(obj)
      lst <- lapply(obj, function(el) {
        if (inherits(el, "sf")) return(el)
        if (is.data.frame(el))  return(df_to_shp(el))
        stop("List elements must be `sf` or data.frame.")
      })
      if (is.null(nms) || any(!nzchar(nms))) {
        names(lst) <- if (!is.null(nms)) {
          # fill any empties
          nms[nchar(nms) == 0] <- paste0(label_prefix, seq_len(sum(nchar(nms) == 0)))
          nms
        } else {
          paste0(label_prefix, seq_along(lst))
        }
      }
      return(lst)
    }

    stop("Input must be an `sf`, data.frame, or a list of those.")
  }

  spp_list <- as_sf_list(x, group = group_x, label_prefix = "sp")
  str_list <- as_sf_list(y, group = group_y, label_prefix = "stressor")

  # Choose or build area of interst
  if (is.null(area)) {
    if (!quiet) message("No area provided; creating AOI from ", area_strategy, " layers (", area_type, ").")
    src <- switch(area_strategy,
                  "stressor" = merge_shp(str_list),
                  "species" = merge_shp(spp_list),
                  "union" = merge_shp(c(spp_list, str_list)))
    area <- create_area(src, area_type = area_type, buffer_frac = area_buffer_frac, quiet = quiet)
  } else if (inherits(area, "bbox")) {
    area <- sf::st_as_sfc(area)
  } else if (inherits(area, "data.frame")) {
    area <- df_to_shp(area)
  } else if (!inherits(area, "sf")) {
    stop("`area` must be NULL, bbox, data.frame, or sf.")
  }

  # Kernels helper
  build_kernels <- function(lst, group_size = NULL, ncls = n_classes) {
    lapply(lst, function(item) {
      get_class_kernel(
        item,
        area = area,
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
  spp_distribution_list <- build_kernels(spp_list, group_size = group_size_x, ncls = 1L)
  str_distribution_list <- build_kernels(str_list, group_size = group_size_y, ncls = 1L)
  spp_kernel_list <- build_kernels(spp_list, group_size = group_size_x, ncls = n_classes)
  stressor_kernel_list <- build_kernels(str_list, group_size = group_size_y, ncls = n_classes)

  # AOI CRS according to return_crs
  if (return_crs == "4326") {
    area <- sf::st_transform(area, 4326)
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

  # Assemble
  out <- list(
    species_distributions = spp_distribution_list,
    stressor_distributions = str_distribution_list,
    species_kernel_maps = spp_kernel_list,
    stressor_kernel_maps = stressor_kernel_list,
    overlap_maps = overlap_maps_list,
    area_of_interest = area
  )
  class(out) <- c("risaMaps", class(out))
  out
}
