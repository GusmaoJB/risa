#' Quick Bycatch Risk Assessment (via HRA)
#'
#' Performs Habitat Risk Assessment (HRA) analysis following the Bycatch Risk
#' Assessment (ByRA) described in Hines et al. (2020) by combining
#' `risa_prep()` to prepare kernel density and overlap maps with
#' `hra()` to perform the InVEST's HRA computation. This wrapper
#' allows rapid testing of different criteria tables, decay settings, and
#' area-of-interest, returning the full HRA results and associated summary statistics.
#'
#' @param x Species or habitat input. Can be an `sf` object, a `data.frame`,
#'   a list of `sf` objects, or a `risaMaps` object produced by `risa_prep()`.
#'   If `x` is a `data.frame` or `sf` object, it can be split into multiple
#'   species/habitat layers using `group_x`. If `x` is a `risaMaps` object,
#'   `y` is ignored and the precomputed maps stored in `x` are used directly.
#' @param y Stressor input: `sf`, data.frame, or list of `sf`.
#'   If a data.frame/`sf`, can be split into multiple layers via `group_y`.
#' @param criteria A `data.frame` formatted according to InVEST's HRA criteria
#'   input, or a named list of such `data.frame`s.
#' @param area Optional AOI polygon (`sf`) or `bbox`/data.frame. If `NULL`,
#'   the AOI is auto-computed based on `area_strategy`, `area_type`,
#'   and `area_buffer_frac`.
#' @param n_classes Integer number of classes for kernel reclassification
#'   (default `3`). Passed to `risa_prep()` and used as `r_max` in `hra()`.
#' @param output_min Numeric or `NULL`. Minimum value used when
#'   `continuous = TRUE`. If `NULL`, defaults to `1`.
#' @param continuous Logical. If `FALSE`, KDE, overlap, and HRA input maps are
#'   based on discrete reclassified rasters. If `TRUE`, KDE and overlap maps are
#'   treated as continuous rescaled rasters. Default is `FALSE`.
#' @param n_overlap Optional integer: number of overlapping stressors used
#'   for classification scaling. By default, it is equal to the number of stressors.
#' @param radius Kernel density bandwidth (projected units). If `NULL`, uses
#'   `radius_method`.
#' @param radius_method Character. Method used to estimate KDE `radius` from `x`
#'   and `y` when `radius = NULL`. Options are `"nndist"`, `"nrd"`,
#'   `"std_distance_scaled"`, `"ppl"`, and `"fixed"`. See help page of `get_class_kernel()` for
#'   details.
#' @param group_x,group_y Optional grouping columns to split `x`/`y` when
#'   they aren’t lists.
#' @param group_size_x,group_size_y Optional numeric weight columns used
#'   as marks in KDE estimation.
#' @param pixel_size Optional pixel size (meters) for the KDE grid; otherwise
#'   `dimyx` is used.
#' @param dimyx Grid size `(ny, nx)` for KDE when `pixel_size` is `NULL`
#'   (default `c(512,512)`).
#' @param exclude_lowest Logical; if `TRUE` (default), the lowest class is
#'   excluded from kernel reclassification to reduce noise.
#' @param lowest_prop Proportion for lowest-value exclusion (default `0.05`).
#' @param area_strategy `"stressor"` (default), `"species"`, or `"union"`.
#'   Determines how the AOI is defined when not provided. Default is defined
#'   after the distribution of the stressors.
#' @param area_type `"convex_hull"` (default) or `"bbox"`. Shape used for
#'   the computed AOI.
#' @param area_buffer_frac Fractional expansion of the computed AOI
#'   (default `0.5`).
#' @param return_crs `"metric"` (default) or `"4326"`. If `"4326"`,
#'   outputs are reprojected to WGS84 geographic coordinates.
#' @param overlap_method Combination rule for stressor overlap: one of
#'   `"product"` (default), `"sum"`, `"geom_mean"`, or `"max"`.
#' @param equation Risk equation used in HRA: `"euclidean"` or
#'   `"multiplicative"`. Passed to `hra()`.
#' @param decay Decay function applied to stressor influence:
#'   `"none"` (default), `"linear"`, `"exponential"`,
#'   `"polynomial_2nd"`, `"polynomial_3rd"`,
#'   `"complementary_decay_2nd"`, or `"complementary_decay_3rd"`.
#'   See `hra()` help page for details.
#' @param buffer_m Named numeric vector of buffer distances (meters)
#'   for each stressor. Required when `decay`is different than "none" for a given stressor.
#' @param quiet Logical; if `TRUE` (default), suppresses console messages
#'   from map preparation.
#'
#' @return
#' An object of class `risaHRA`, as returned by `hra()`, with additional
#' elements:
#' \describe{
#'   \item{`kde_maps`}{The `risaMaps` object produced by `risa_prep()` or
#'   supplied directly through `x`.}
#'   \item{`area_of_interest`}{The AOI polygon used for map preparation and
#'   risk analysis.}
#' }
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Prepares KDE, distribution, and overlap maps with `risa_prep()`,
#'   unless `x` is already a `risaMaps` object.
#'   \item Reshapes the criteria table using `criteria_reshape()`.
#'   \item Reshapes prepared maps with `reshape_risa_maps()`.
#'   \item Runs `hra()` to compute exposure, consequence, and risk maps.
#' }
#'
#' Internally, `risa_prep()` is called with `output_layer_type = "both"` so that
#' both raster and polygon outputs are available when discrete maps are used.
#' When `continuous = TRUE`, raster outputs are used.
#'
#' @seealso `risa_prep()`, `hra()`
#'
#' @references
#' Hines, E., Ponnampalam, L. S., Junchompoo, C., Peter, C., Vu, L.,
#' Huynh, T., et al. (2020). Getting to the bottom of bycatch: a GIS-based toolbox
#' to assess the risk of marine mammal bycatch. Endang. Spec. Res. 42, 37–57.
#' doi: 10.3354/esr01037
#'
#' @examples
#' \dontrun{
#' quick_byra(x = species_sf,
#'            y = stressors_sf,
#'            criteria = my_criteria,
#'            n_classes = 4,
#'            overlap_method = "geom_mean",
#'            decay = "linear",
#'            buffer_m = c(fishing = 5000))
#' }
#' @export
quick_byra <- function(x,
                        y = NULL,
                        criteria,
                        area = NULL,
                        n_classes = 3,
                        output_min = NULL,
                        continuous = FALSE,
                        n_overlap = NULL,
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
                        equation = c("euclidean","multiplicative"),
                        decay = c("none", "linear", "exponential",
                                  "polynomial_2nd", "polynomial_3rd",
                                  "complementary_decay_2nd", "complementary_decay_3rd"),
                        buffer_m = NULL,
                        quiet = TRUE) {

  radius_method <- match.arg(radius_method)
  return_crs <- match.arg(return_crs)
  overlap_method <- match.arg(overlap_method)
  equation <- match.arg(equation)
  decay <- match.arg(decay)
  area_strategy <- match.arg(area_strategy)
  area_type <- match.arg(area_type)
  input_maps <- NULL

  # Reshape criteria input for analysis
  criteria <- criteria_reshape(criteria)
  sample_crit <- criteria[[1]]
  crit_names <- unique(sample_crit[is.na(sample_crit$RATING),"ATTRIBUTES"])

  # If x is a risaMaps object
  if (inherits(x, "risaMaps")) {
    if(!quiet) message("Input x is a risaMaps object.")
    if(!is.null(y)){
      if(!quiet) message("Input y is not null. Ignoring y object...")
    }
    input_maps <- x
  } else {
    input_maps <- risa_prep(
      x = x,
      y = y,
      area = area,
      n_classes = n_classes,
      output_min = output_min,
      output_layer_type = "both",
      radius = radius,
      radius_method = radius_method,
      group_x = group_x,
      group_y = group_y,
      group_size_x = group_size_x,
      group_size_y = group_size_y,
      pixel_size = pixel_size,
      dimyx = dimyx,
      exclude_lowest = exclude_lowest,
      lowest_prop = lowest_prop,
      area_strategy = area_strategy,
      area_type = area_type,
      area_buffer_frac = area_buffer_frac,
      continuous = continuous,
      return_crs = return_crs,
      overlap_method = overlap_method,
      quiet = quiet
    )
  }

  # Reshape lists for HRA analysis
  raster_list <- reshape_risa_maps(input_maps, crit_names)
  species_distr <- input_maps$species_distributions

  # Perform HRA
  byra_hra <- hra(
    raster_list,
    species_distr,
    criteria,
    equation,
    r_max = n_classes,
    n_overlap,
    output_decimal_crs = return_crs != "metric",
    decay,
    buffer_m)

  # Outputs
  output <- byra_hra
  output[["kde_maps"]] <- input_maps
  output$area_of_interest <- input_maps$area_of_interest
  return(output)
}
