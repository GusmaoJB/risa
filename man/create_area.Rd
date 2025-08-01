% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_area.R
\name{create_area}
\alias{create_area}
\title{Create polygon of the area of interest}
\usage{
create_area(
  x,
  crs = NULL,
  area_type = c("convex_hull", "bbox"),
  buffer_frac = 0.5
)
}
\arguments{
\item{x}{An `sf` object or data frame with longitude/latitude in the first two columns.}

\item{crs}{Integer or string; EPSG code or proj4 string for metric transformation. If `NULL`, select UTM zone automatically.}

\item{area_type}{Character; one of \code{"convex_hull"} or \code{"bbox"}. Default is \code{"convex_hull"}.}

\item{buffer_frac}{Numeric; fraction by which to buffer the hull or bounding box. (e.g. 0.5 for 50\%.)}
}
\value{
An `sf` object representing the area polygon.
}
\description{
Generate an area polygon (convex hull or bounding box) around input points, with optional buffering.
}
\examples{
# Create test data
coords <- data.frame(long = c(1,2,2,4), lat = c(4,4,2,2))
coords_vec <- df_to_shp(coords)
coords_vec_metric <- transform_to_metric(coords_vec)$shape

# Create area of interest (default: convex hull)
aoi <- create_area(coords)

# Plot results
plot(st_geometry(aoi), border = "blue", axes=TRUE)
plot(st_geometry(coords_vec_metric), add = TRUE, col = "red")
}
