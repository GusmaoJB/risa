% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/risa_prep.R
\name{risa_prep}
\alias{risa_prep}
\title{Prepare RISA kernel and overlap maps}
\usage{
risa_prep(
  x,
  y,
  area = NULL,
  n_classes = 3,
  output_layer_type = "both",
  radius = NULL,
  output_decimal_crs = FALSE
)
}
\arguments{
\item{x}{Species input as `sf`, data frame, or list.}

\item{y}{Stressor input as `sf`, data frame, or list.}

\item{area}{Optional area polygon (`sf`); auto-computed if `NULL`.}

\item{n_classes}{Integer number of classes.}

\item{output_layer_type}{Character; one of `"shp"`, `"raster"`, or `"both"`.}

\item{radius}{Numeric smoothing bandwidth.}

\item{output_decimal_crs}{Logical; if `TRUE`, reproject outputs to EPSG:4326.}
}
\value{
A list with lists for `species_kernel_maps`, `stressor_kernel_maps`, `overlap_maps`, and `area_of_interest`.
}
\description{
Generates kernel density maps for species and stressors, merges, and computes overlap within an area of interest.
}
\examples{
# Creating test data
spp_df <- rbind(data.frame(long = rnorm(80, 0, 10),
                           lat = rnorm(80, 0, 10), species = "sp1"),
                data.frame(long = rnorm(60, 0, 10),
                           lat = rnorm(60, 0, 10), species = "sp2"))
str_df <- rbind(data.frame(long = rnorm(100, 0, 5),
                           lat = rnorm(100, 0, 10), stressor = "trawling"),
                data.frame(long = rnorm(50, 0, 10),
                           lat = rnorm(100, 0, 5), stressor = "gillnet"))

# Create kernel maps of species and stressor distributions and overlap maps
risa_maps <- risa_prep(spp_df, str_df)

# Species and Stressor distributions
dev.off()
par(mfrow = c(2,2))

plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species 1")
plot(risa_maps$species_kernel_maps$sp1$shp, add = TRUE)

plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species 2")
plot(risa_maps$species_kernel_maps$sp2$shp, add = TRUE)

plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Gillnet")
plot(risa_maps$stressor_kernel_maps$gillnet$shp, add = TRUE, col=c("lightgreen", "green", "darkgreen"))

plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Trawling")
plot(risa_maps$stressor_kernel_maps$trawling$shp, add = TRUE, col=c("lightgreen", "green", "darkgreen"))

# Overlap maps
plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species1 vs Gillnet")
plot(risa_maps$overlap_maps$sp1$gillnet, add = TRUE, col=c("yellow", "orange", "red"))

plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species2 vs Gillnet")
plot(risa_maps$overlap_maps$sp2$gillnet, add = TRUE, col=c("yellow", "orange", "red"))

plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species1 vs Trawling")
plot(risa_maps$overlap_maps$sp1$gillnet, add = TRUE, col=c("yellow", "orange", "red"))

plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species2 vs Trawling")
plot(risa_maps$overlap_maps$sp2$gillnet, add = TRUE, col=c("yellow", "orange", "red"))
}
