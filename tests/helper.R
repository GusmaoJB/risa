# Stressor decay kernel in the style of InVEST.
# stressor_layers: list of SpatRaster layers for a given stressor (attributes); 1 where stressor occurs, NA elsewhere
# template: a SpatRaster to align to (use species distribution)
stressor_decay_kernel <- function(stressor_layers, template, buffer_m, decay = c("none","linear","exponential"),
                                  cutoff = 0.01) {
  decay <- match.arg(decay)
  if (is.null(buffer_m) || is.na(buffer_m)) return(1)  # no buffer provided
  
  # Mosaic all stressor attribute rasters -> presence mask (1 on stressor, NA elsewhere)
  collection <- terra::sprc(stressor_layers)
  r_mos <- terra::mosaic(collection, fun = max)
  r_mos <- terra::ifel(!is.na(r_mos), 1, NA)
  r_mos <- if (inherits(template, "SpatRaster")) terra::project(r_mos, template, method = "near") else r_mos
  
  # If zero buffer -> keep only the footprint (value 1 on stressor, NA elsewhere)
  if (buffer_m == 0) {
    return(r_mos)
  }
  
  # Distance (map units of CRS!). For meters like InVEST, use a projected CRS.
  # terra::distance() returns distance to non-NA cells in 'r_mos'
  d <- terra::distance(r_mos)
  
  # Build kernel
  B <- as.numeric(buffer_m)
  if (decay == "none") {
    w <- terra::ifel(d <= B, 1, NA)
  } else if (decay == "linear") {
    w <- 1 - d / B
    w <- terra::ifel(d <= B, w, NA)
  } else { # exponential
    k <- -log(cutoff) / B     # so that w(B) = cutoff
    w <- exp(-k * d)
    w <- terra::ifel(d <= B, w, NA)  # outside buffer -> NA
  }
  
  # (Optional) tiny values to 0 (cosmetic)
  w <- terra::ifel(w < 1e-6, 0, w)
  
  # Ensure exact 1 on stressor footprint
  w <- terra::cover(terra::ifel(!is.na(r_mos), 1, NA), w)
  
  w
}
