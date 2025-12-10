#' Reshape nested species–stressor maps into a standardized structure
#'
#' This function reorganizes a complex list of nested species and stressor maps
#' (a `risaMaps` object, output from `risa_maps`) into a simplified structure. For each
#' species–stressor pair, it extracts the stressor intensity map
#' (`stressor_kernel_maps`) and the overlap map (`overlap_maps`) and returns
#' them as a named list with user-defined criteria names.
#'
#' @param risa_maps A nested list (class `"risaMaps"`) containing species
#'   distributions, stressor distributions, kernel maps, overlap maps, and other
#'   elements as produced by the \pkg{risa} workflow.
#' @param criteria_names Character vector of length 2. Names assigned to the two
#'   elements in each species–stressor pair: the first corresponds to the
#'   stressor intensity map (from `stressor_kernel_maps`), and the second to the
#'   likelihood of interaction map (from `overlap_maps`). For example:
#'   `c("intensity", "likelihood of interaction")`.
#' @param only_species Optional character vector. Subset of species names to
#'   include. Defaults to all available species.
#' @param only_stressors Optional character vector. Subset of stressor names to
#'   include. Defaults to all available stressors.
#'
#' @return A nested list structured as
#'   \preformatted{
#'   list(
#'     species1 = list(
#'       stressor1 = list(
#'         <criteria_names[1]> = SpatRaster,
#'         <criteria_names[2]> = SpatRaster
#'       ),
#'       stressor2 = list(...)
#'     ),
#'     species2 = list(...)
#'   )
#'   }
#'   Species or stressors with no valid data are omitted.
#'
#' @examples
#' \dontrun{
#' new_maps <- reshape_risa_maps(
#'   risa_maps,
#'   criteria_names = c("intensity", "likelihood of interaction")
#' )
#'
#' # Subset to only species1 and stressor2
#' new_maps <- reshape_risa_maps(
#'   risa_maps,
#'   criteria_names = c("intensity", "likelihood of interaction"),
#'   only_species = "species1",
#'   only_stressors = "stressor2"
#' )
#' }
#'
#' @export
reshape_risa_maps <- function(risa_maps,
                              criteria_names,
                              only_species = NULL,
                              only_stressors = NULL) {
  # Checks
  if (!inherits(risa_maps, "risaMaps")) {
    stop("Input map list must be a 'risaMaps' object.")
  }
  if (missing(criteria_names)) {
    stop("`criteria_names` is required and must be a character vector of length 2.")
  }
  if (!is.character(criteria_names) || length(criteria_names) != 2L) {
    stop("`criteria_names` must be a character vector of length 2, e.g. c('intensity','likelihood of interaction').")
  }
  if (any(!nzchar(criteria_names))) {
    stop("`criteria_names` cannot contain empty strings.")
  }
  if (anyDuplicated(criteria_names)) {
    stop("`criteria_names` must be unique.")
  }

  # Helper
  get2 <- function(x, path) {
    for (p in path) {
      if (is.null(x) || is.null(x[[p]])) return(NULL)
      x <- x[[p]]
    }
    x
  }
  `%||%` <- function(a, b) if (is.null(a)) b else a

  # Discover available keys
  species_all   <- names(get2(risa_maps, c("overlap_maps"))) %||% character()
  stressors_all <- names(get2(risa_maps, c("stressor_kernel_maps"))) %||% character()
  if (length(species_all) == 0L || length(stressors_all) == 0L) return(list())

  # Optional filters
  species   <- if (is.null(only_species)) species_all else intersect(species_all, only_species)
  stressors <- if (is.null(only_stressors)) stressors_all else intersect(stressors_all, only_stressors)

  # Build output: species -> stressor -> {criteria_names[1], criteria_names[2]}
  out <- setNames(lapply(species, function(sp) {
    per_sp <- lapply(stressors, function(st) {
      exposure <- get2(risa_maps, c("stressor_kernel_maps", st, "raster"))
      overlap <- get2(risa_maps, c("overlap_maps", sp, st, "raster"))
      if (is.null(exposure) || is.null(overlap)) return(NULL)
      stats::setNames(list(exposure, overlap), criteria_names)
    })
    names(per_sp) <- stressors
    per_sp <- per_sp[!vapply(per_sp, is.null, FALSE)]
    if (length(per_sp) == 0L) NULL else per_sp
  }), species)

  out[!vapply(out, is.null, FALSE)]
}
