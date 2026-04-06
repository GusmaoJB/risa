# Stressor mapping: from google earth to risa
library(risa)
library(sf)
library(terra)
library(dplyr)

#### For just one kml object
# load file
x <- read_sf("/home/jojo/Downloads/Untitled map.kml")

# Reproject to metric CRS
x_m <- transform_to_metric(x)

# Convert to SpatVector object
v_m <- vect(x_m$shape)

# Giving Rating attribute to object
v_m$Rating <- 1

# choose pixel size in meters (e.g., 100 m)
r <- rast(ext(v_m), res = 100, crs = crs(v_m))

ov_sum <- rasterize(v_m, r, field = "Rating", fun = "sum", background = NA)

plot(ov_sum)

# For multiple kml objects

dir_kml <- "/path/to/your/folder"

# 1) Find files
files <- list.files(
  dir_kml,
  pattern = "\\.(kml|kmz)$",
  full.names = TRUE,
  ignore.case = TRUE
)

# 2) Read & combine into one sf
x_all <- files |>
  lapply(\(f) {
    st_read(f, quiet = TRUE) |>
      mutate(source_file = basename(f))
  }) |>
  bind_rows()

# Optional: keep only polygons (drop points/lines if any)
x_all <- x_all |> filter(st_geometry_type(geometry) %in% c("POLYGON", "MULTIPOLYGON"))

# 3) Add rating
x_all$Rating <- 1

# Transform to metric
x_m <- transform_to_metric(x_all)

# 4) Convert to terra vector
v_m <- vect(x_m$shp)

# 6) Template raster + overlap sum
r <- rast(ext(v_m), res = 100, crs = crs(v_m))   # 100 m pixels (adjust)
ov_sum <- rasterize(v_m, r, field = "Rating", fun = "sum", background = 0)

plot(ov_sum)



