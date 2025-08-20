#' Generate decay coefficient rasters
#'
#' Computes spatially explicit decay coefficients around non-NA cells of a
#' reference raster. The decay describes how values decrease with distance
#' from the source cells, up to a specified buffer distance (meters).
#'
#' Supported decay types:
#' \describe{
#'   \item{`none`}{No decay. All cells within the buffer are set to 1, outside 0.}
#'   \item{`linear`}{Linear decrease from 1 at the source to 0 at the buffer edge.}
#'   \item{`exponential`}{Exponential decrease from 1 at the source to near 0 at the buffer edge
#'   (rate chosen so that value ≈ 1e-6 at `buffer_m`).}
#'   \item{`polynomial_2nd`}{Quadratic polynomial decay: \eqn{(1 - d / buffer_m)^2}, convex decreasing.}
#'   \item{`polynomial_3rd`}{Cubic polynomial decay: \eqn{(1 - d / buffer_m)^3}, convex decreasing.}
#'   \item{`complementary_decay_2nd`}{Complementary quadratic: \eqn{1 - (d / buffer_m)^2},
#'   flat near 1 then dropping sharply near buffer edge.}
#'   \item{`complementary_decay_3rd`}{Complementary cubic: \eqn{1 - (d / buffer_m)^3},
#'   even flatter near 1, sharper drop near buffer edge.}
#' }
#'
#' @param raster A [terra::SpatRaster] with non-NA source cells from which distances are measured.
#' @param raster2 Optional [terra::SpatRaster] to mask the output decay coefficients (only keep cells where `raster2` is not NA).
#' @param decay Character string specifying the decay function (see above).
#' @param buffer_m Numeric. Buffer distance in meters. Beyond this distance decay values are set to 0 (or NA if masked).
#'
#' @return A [terra::SpatRaster] with decay coefficients in \[0, 1\].
#'
#' @details
#' The function requires a projected (metric) CRS. Geographic (lon/lat) rasters must
#' be reprojected before use.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r <- rast(nrows=100, ncols=100, res=100)
#' r[50,50] <- 1
#' dmap <- decay_coeffs(r, decay="exponential", buffer_m=2000)
#' plot(dmap)
#' }
#' @export
decay_coeffs <- function(raster, raster2 = NULL,
                         decay = c("none","linear","exponential",
                                   "polynomial_2nd","polynomial_3rd",
                                   "complementary_decay_2nd","complementary_decay_3rd"),
                         buffer_m) {
  stopifnot(inherits(raster, "SpatRaster"))
  cr <- terra::crs(raster, describe = TRUE)
  if (isTRUE(cr$is_lonlat)) {
    stop("CRS is geographic (lon/lat). Reproject rasters to a metric CRS (meters) first.")
  }

  decay <- match.arg(decay)

  # distance from non-NA cells of 'raster'
  d <- terra::distance(raster)

  # limit to buffer -> outside becomes NA
  d <- terra::ifel(d > buffer_m, NA, d)
  within_buf <- !is.na(d)

  # choose decay (all return 0 outside buffer)
  decay_map <- switch(
    decay,
    "none"              = terra::ifel(within_buf, 1, 0),
    "linear"            = terra::ifel(within_buf, 1 - (d / buffer_m), 0),
    "exponential"       = { k <- -log(1e-6) / buffer_m; terra::ifel(within_buf, exp(-k * d), 0) },
    "polynomial_2nd"    = terra::ifel(within_buf, (1 - d / buffer_m)^2, 0),
    "polynomial_3rd"    = terra::ifel(within_buf, (1 - d / buffer_m)^3, 0),
    "complementary_decay_2nd" = terra::ifel(within_buf, 1 - (d / buffer_m)^2, 0),
    "complementary_decay_3rd" = terra::ifel(within_buf, 1 - (d / buffer_m)^3, 0)
  )

  # If you prefer NA outside the buffer, replace the 0’s with NA via:
  # decay_map <- terra::ifel(within_buf, decay_map, NA)

  if (is.null(raster2)) {
    return(decay_map)  # already limited to buffer above
  } else {
    # keep only cells where raster2 is non-NA
    return(terra::mask(decay_map, raster2))
  }
}
