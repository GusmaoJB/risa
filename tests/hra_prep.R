library(sf)


byra_prep <- function(spp_list, str_list,
                     area = NULL,
                     pixel_size = NULL,
                     dimyx = c(512,512),
                     reclass = NULL,
                     reclass_cat = FALSE,
                     crs = NULL,
                     quiet = TRUE) {
  print("Nice!")
}

# Basic checks
# Check area
if (inherits(area, "bbox")) {
  area <- sf::st_as_sfc(area)
} else if (is.data.frame(area)) {
  area <- df_to_shp(area)
} else if (!inherits(area, "sf")) {
  stop("`area` must be NULL, bbox, data.frame, or sf.")
}

# Check if inputs are 'sf' objects and convert them to 'SpatRaster'
lapply()

# Check raster inputs
spp_list <- .as_raster(spp_list)
if (!inherits(sp_distr, "SpatRaster")) stop("'spp_list' must be or contain a SpatRaster object.")
str_list <- .as_raster(str_list)
if (!inherits(sp_distr, "SpatRaster")) stop("'str_list' must be or contain a SpatRaster object.")

# Defining important objects
all_rasters <- c(spp_list, str_list)

# Check if all rasters have the same crs
if (!.same_crs(all_rasters)){
  if (!quiet) {
    message("Input rasters do not have the same coordinate system.")
    if (is.null(area)) {
      stop("Cannot define a template crs without area.")
    } else {
      message("Using area's crs as template: reprojecting rasters...")
      longlat <- sf::st_is_longlat(area)
      if (longlat) {
        area <- transform_to_metric(area)
      }
      spp_list <- lapply(spp_list,
        function(x) align_to(x, template = template, categorical = terra::is.int(x)))
      str_list <- lapply(str_list,
        function(y) align_to(y, template = template, categorical = terra::is.int(y)))
      all_rasters <- lapply(all_rasters,
        function(z) align_to(z, template = template, categorical = terra::is.int(z)))
    }
  }
} else if (is.null(area)) {
  if (!quiet) message("No area provided; creating AOI (bounding box) from all layers.")
  e_all <- do.call(ext, all_rasters)

  # Estimate buffer as 5% of the diagonal of the combined extent
  dx <- xmax(e) - xmin(e)
  dy <- ymax(e) - ymin(e)
  diag_len <- sqrt(dx^2 + dy^2)
  buf <- 0.05 * diag_len

  e_buf <- ext(
    xmin(e) - buf,
    xmax(e) + buf,
    ymin(e) - buf,
    ymax(e) + buf)

  area <- st_bbox(c(
    xmin = xmin(e_buf),
    ymin = ymin(e_buf),
    xmax = xmax(e_buf),
    ymax = ymax(e_buf)
  ), crs = st_crs(all_rasters[[1]]))
}

# If







# Helpers
# Checks if all rasters in a list have the same crs
.same_crs <- function(rlist) {
  crs_list <- lapply(rlist, crs)
  ref <- crs_list[[1]]
  all(vapply(crs_list, function(x) x == ref, logical(1)))
}

# Find/extract a SpatRaster if user passed a list container
.as_raster <- function(x) {
  if (inherits(x, "SpatRaster")) return(x)
  if (is.list(x)) {
    for (item in x) {
      r <- .as_raster(item)
      if (inherits(r, "SpatRaster")) return(r)
    }
  }
  NULL
}

# Scale continuous values
.scale_vars <- function(x, mim_max) {
  if (!is.numeric(min_max)) {
    stop("Error: min_max must be a numeric vector.")
  }
  if (mim_max[1] >= mim_max[2]) {
    stop("Error: minimum value for rescaling is not valid.")
  }
  scaled_x <- mim_max[1] + (x-min(x) * (mim_max[2] - mim_max[1])) / (max(x) - min(x))
  return(round(scaled_x, 4))
}









# Prepare a list of rasters to hra
# Checks if all rasters have valid crs



  # if yes, checks if area has a valid crs
    # if not, guess its crs
      # if cannot guess, stop
# 2 - If there is no area, try define it from bounding box of all rasters
  # Checks if all rasters in each list have the same crs
    # If not, reproject all rasters to a common crs
      # if crs is NULL, grab the first crs of the first stressor rasters as template
        # reproject all rasters
# 3 - If reclass is a vector, then apply a transformation








# Check if spp_list and str_list are lists of rasters (SpatRaster objects)




# Accept list containers for sp_distr
sp_distr <- .as_raster(sp_distr)
if (!inherits(sp_distr, "SpatRaster")) stop("'species_distr' must be or contain a SpatRaster.")





# takes as input SpatRaster, sf, or data.frame objects.
# If area is not provided, it will be derived as a convex hull from the stressor blob
# All rasters and shp will be reprojected after the CRS of the area of interest
# All rasters and shp will be cropped after the area of interest
# The output is a list ready for hra() function
