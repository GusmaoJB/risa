#' Create a stressor decay surface from a binary occurrence raster
#'
#' Computes a distance-based decay surface from a stressor occurrence
#' raster (`SpatRaster`), applying a user-defined decay function within
#' a specified buffer distance. The resulting raster contains weights
#' ranging from 0 to 1, representing the relative influence of the stressor
#' as a function of distance from stressor presence cells.
#'
#' The function requires a projected (metric) CRS. If the raster is in
#' longitude/latitude coordinates, it must be reprojected before use.
#'
#' @param stress_occ A `SpatRaster` indicating stressor occurrence.
#'   Non-`NA` cells represent stressor presence.
#' @param buffer_m Numeric. Buffer distance in meters defining the maximum
#'   spatial influence of the stressor.
#' @param decay Character. Type of decay function to apply. Options are:
#'   `"none"`, `"linear"`, `"exponential"`, `"polynomial_2nd"`,
#'   `"polynomial_3rd"`, `"complementary_decay_2nd"`,
#'   `"complementary_decay_3rd"`.
#' @param habitat_mask Optional `SpatRaster`. If provided, the resulting
#'   decay surface is masked to habitat pixels.
#'
#' @return A `SpatRaster` with values between 0 and 1 representing the
#'   stressor influence surface. Returns `NULL` if `buffer_m` is `NULL`
#'   or `NA`.
#'
#' @details
#' Distance is computed using `terra::distance()`, and pixel size is
#' derived from raster resolution. Polynomial and complementary decay
#' functions are extensions not implemented in InVEST.
#'
#' @examples
#' \dontrun{
#' w <- make_stressor_decay(
#'   stress_occ = stress_raster,
#'   buffer_m = 1000,
#'   decay = "linear",
#'   habitat_mask = habitat_raster
#' )
#' }
#'
#' @export
make_stressor_decay <- function(stress_occ, buffer_m, decay, habitat_mask = NULL) {
  stopifnot(inherits(stress_occ, "SpatRaster"))
  cr <- terra::crs(stress_occ, describe = TRUE)
  if (isTRUE(cr$is_lonlat)) stop("Reproject to a metric CRS first.")

  # distance to non-NA cells (stressor presence); units = meters (CRS units)
  d_m <- terra::distance(stress_occ)

  # pixel size in meters (assume square pixels - not a problem if risa_prep was used)
  rxy <- terra::res(stress_occ)
  pix_m <- as.numeric(rxy[1])

  if (is.null(buffer_m) || is.na(buffer_m)) {
    return(NULL)
  }

  if (buffer_m == 0) {
    w <- terra::ifel(terra::is.na(stress_occ), 0, 1)
  } else {
    d_pix <- d_m / pix_m
    B_pix <- buffer_m / pix_m

    w <- switch(match.arg(decay, c("none","linear","exponential",
                                   "polynomial_2nd","polynomial_3rd",
                                   "complementary_decay_2nd","complementary_decay_3rd")),
                "none"  = terra::ifel(d_pix <= B_pix, 1, 0),

                "linear" = {
                  w0 <- 1 - (d_pix / B_pix)
                  terra::ifel(d_pix <= B_pix, terra::ifel(w0 < 0, 0, w0), 0)
                },

                "exponential" = {
                  w0 <- exp(-d_pix)
                  w0 <- terra::ifel(w0 < 1e-6, 0, w0)
                  terra::ifel(d_pix <= B_pix, w0, 0)
                },

                # These are extensions of decay functions (not present in InVEST)
                "polynomial_2nd" = {
                  w0 <- (1 - d_pix / B_pix)^2
                  terra::ifel(d_pix <= B_pix, terra::ifel(w0 < 0, 0, w0), 0)
                },
                "polynomial_3rd" = {
                  w0 <- (1 - d_pix / B_pix)^3
                  terra::ifel(d_pix <= B_pix, terra::ifel(w0 < 0, 0, w0), 0)
                },
                "complementary_decay_2nd" = {
                  w0 <- 1 - (d_pix / B_pix)^2
                  terra::ifel(d_pix <= B_pix, terra::ifel(w0 < 0, 0, w0), 0)
                },
                "complementary_decay_3rd" = {
                  w0 <- 1 - (d_pix / B_pix)^3
                  terra::ifel(d_pix <= B_pix, terra::ifel(w0 < 0, 0, w0), 0)
                }
    )
  }

  # Limit to habitat pixels if a mask is provided
  if (!is.null(habitat_mask)) {
    w <- terra::mask(w, habitat_mask)
  }
  w
}
