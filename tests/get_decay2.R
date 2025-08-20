#' Apply decay coefficients to a raster
#'
#' Combines a raster of values with a raster of decay coefficients,
#' scaling cell values according to the decay. Useful for applying
#' distance-based attenuation to spatial attributes.
#'
#' @param raster A `terra::SpatRaster` containing the base values (e.g., stressor intensity).
#' @param decay_coefficients A `terra::SpatRaster` with decay coefficients in \[0, 1\],
#' typically obtained from `decay_coeffs()`.
#' @param start_value Numeric scalar. Value to assign where `raster` has zeros but
#' is within the buffer (multiplied by decay).
#'
#' @return A `terra::SpatRaster` with decayed values. NA is returned where either
#' `raster` or `decay_coefficients` is NA.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r <- rast(nrows=100, ncols=100, res=100)
#' r[50,50] <- 1
#' dmap <- decay_coeffs(r, decay="linear", buffer_m=2000)
#' out <- get_decay_map(r, dmap, start_value=1)
#' plot(out)
#' }
#' @export
get_decay_map <- function(raster, decay_coefficients, start_value) {
  stopifnot(inherits(raster, "SpatRaster"),
            inherits(decay_coefficients, "SpatRaster"),
            length(start_value) == 1, is.finite(start_value))

  if (!terra::compareGeom(raster, decay_coefficients, stopOnError = FALSE))
    stop("`decay_coefficients` and `raster` must have identical extent/resolution/CRS.")

  # keep NA where either side is NA
  terra::ifel(
    is.na(raster) | is.na(decay_coefficients), NA,
    terra::ifel(raster == 0, start_value * decay_coefficients,
                raster      * decay_coefficients)
  )
}

