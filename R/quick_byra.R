#' Quick Bycatch Risk Assessment (via HRA)
#'
#' Performs Habitat Risk Assessment (HRA) analysis following the Bycatch Risk
#' Assessment (ByRA) described in Hines et al. (2020) by combining
#' `risa_prep()` to prepare kernel density and overlap maps with
#' `hra()` to perform the InVEST's HRA computation. This wrapper
#' allows rapid testing of different criteria tables, decay settings, and
#' area-of-interest, returning the full HRA results and associated summary statistics.
#'
#' @param x Species (habitat) input: `sf`, data.frame of individual occurrences, or list of `sf`.
#'   If a data.frame/`sf`, can be split into multiple layers via `group_x`.
#' @param y Stressor input: `sf`, data.frame, or list of `sf`.
#'   If a data.frame/`sf`, can be split into multiple layers via `group_y`.
#' @param criteria A `data.frame` (formated after InVEST's HRA criteria input)
#'   or a named list of `data.frame`s (ecosystem) with columns `STRESSOR`,
#'   `ATTRIBUTES`, `RATING`, `DQ`, `WEIGHT`, and `E/C`. These tables define the
#'   exposure/consequence criteria for each stressor.
#' @param area Optional AOI polygon (`sf`) or `bbox`/data.frame. If `NULL`,
#'   the AOI is auto-computed based on `area_strategy`, `area_type`,
#'   and `area_buffer_frac`.
#' @param n_classes Integer number of classes for kernel reclassification
#'   (default `3`). Passed to `risa_prep()` and used as `r_max` in `hra()`.
#' @param n_overlap Optional integer: number of overlapping stressors assumed
#'   for classification scaling. By default inferred from the criteria.
#' @param radius Kernel density bandwidth (projected units). If `NULL`, uses
#'   `radius_method`.
#' @param radius_method One of `"nndist"` (default), `"ppl"`, or `"fixed"`.
#'   Bandwidth selection method for KDE. See help page of `get_class_kernel()` for
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
#' @return An object of class "risaHRA" (see `hra()` for details),
#'   with an added element: `area_of_interest` (The AOI polygon used for analysis).
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Prepares kernel density and overlap maps with `risa_prep()`.
#'   \item Reshapes the stressor–criteria structure for HRA.
#'   \item Runs `hra()` to compute exposure, consequence, and risk maps.
#' }
#'
#' It is intended as a convenience wrapper for exploratory analyses,
#' rather than full customization.
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
quick_byra <- function(x, y,
                       criteria,
                       area = NULL,
                       n_classes = 3,
                       n_overlap = NULL,
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

  # Generate Kernel Density and ditribution maps
  input_maps <- risa_prep(x, y,
                          area,
                          n_classes,
                          output_layer_type = "both",
                          radius,
                          radius_method,
                          group_x,
                          group_y,
                          group_size_x,
                          group_size_y,
                          pixel_size,
                          dimyx,
                          exclude_lowest,
                          lowest_prop,
                          area_strategy,
                          area_type,
                          area_buffer_frac,
                          return_crs,
                          overlap_method,
                          quiet)

  # Reshape lists and criteria for HRA analysis
  criteria <- criteria_reshape(criteria)
  sample_crit <- criteria[[1]]
  crit_names <- unique(sample_crit[is.na(sample_crit$RATING),"ATTRIBUTES"])

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
  output$area_of_interest <- input_maps$area_of_interest
  return(output)
}
