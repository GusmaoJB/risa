decay_coeffs <- function(species_kd, stressor_kd, buffer_m,
                        decay = c("none", "linear", "exponential")) {
  decay <- match.arg(decay)

  # Basic checks
  stopifnot(inherits(species_kd, "SpatRaster"),
            inherits(stressor_kd, "SpatRaster"))

  if (!terra::compareGeom(species_kd, stressor_kd, stopOnError = FALSE)) {
    stop("`species_kd` and `stressor_kd` must have identical geom (extent, res, CRS).")
  }

  # Check if CRS is metric
  cr <- terra::crs(species_kd, describe = TRUE)
  if (isTRUE(cr$is_lonlat)) {
    stop("CRS is geographic (lon/lat). Reproject both rasters to a metric CRS (meters) first.")
  }

  # Masks describing occurrences of species and stressors
  spp_mask  <- terra::classify(species_kd, rcl = matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = TRUE))
  str_mask  <- terra::classify(stressor_kd, rcl = matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = TRUE))

  # distance (meters) from every cell to nearest stressor cell (>0)
  d <- terra::distance(str_mask)
  d <- terra::ifel(d > buffer_m, NA, d)

  # Decay inside buffer
  within_buf <- d <= buffer_m

  Djkl <- switch(
    decay,
    "none" = terra::ifel(within_buf, 1, 0),
    "linear" = terra::ifel(within_buf, 1 - (d / buffer_m), 0),
    "exponential" = terra::ifel(within_buf, 1 - exp(log10(1e-6) / d), 0)
  )

  # set values to NA where species is absent/NA
  Djkl <- terra::mask(Djkl, spp_mask)

  names(Djkl) <- paste0("Djkl_", decay)
  Djkl
}


# Add buffer to target raster given another rasters with decay coefficients
get_decay <- function(raster, decay_coefficients, start_value) {
  stopifnot(inherits(raster, "SpatRaster"),
            inherits(decay_coefficients, "SpatRaster"))
  if (!terra::compareGeom(raster, decay_coefficients, stopOnError = FALSE))
    stop("`decay_coefficients` and `raster` must have identical extent/resolution/CRS.")

  pos_only  <- terra::ifel(raster > 0, raster, NA)

  out <- terra::ifel(raster == 0, decay_coefficients * start_value, raster * decay_coefficients)

  names(out) <- names(raster)
  return(out)
}
