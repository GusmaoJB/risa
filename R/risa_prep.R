#' Prepare RISA kernel and overlap maps
#'
#' Generates kernel density maps for species and stressors, merges, and computes overlap within an area of interest.
#'
#' @param x Species input as `sf`, data frame, or list.
#' @param y Stressor input as `sf`, data frame, or list.
#' @param area Optional area polygon (`sf`); auto-computed if `NULL`.
#' @param n_classes Integer number of classes.
#' @param output_layer_type Character; one of `"shp"`, `"raster"`, or `"both"`.
#' @param radius Numeric smoothing bandwidth.
#' @param output_decimal_crs Logical; if `TRUE`, reproject outputs to EPSG:4326.
#' @return A list with lists for `species_kernel_maps`, `stressor_kernel_maps`, `overlap_maps`, and `area_of_interest`.
#' @importFrom sf st_crs
#' @examples
#' # Creating test data
#' spp_df <- rbind(data.frame(long = rnorm(80, 0, 10),
#'                            lat = rnorm(80, 0, 10), species = "sp1"),
#'                 data.frame(long = rnorm(60, 0, 10),
#'                            lat = rnorm(60, 0, 10), species = "sp2"))
#' str_df <- rbind(data.frame(long = rnorm(100, 0, 5),
#'                            lat = rnorm(100, 0, 10), stressor = "trawling"),
#'                 data.frame(long = rnorm(50, 0, 10),
#'                            lat = rnorm(100, 0, 5), stressor = "gillnet"))
#'
#' # Create kernel maps of species and stressor distributions and overlap maps
#' risa_maps <- risa_prep(spp_df, str_df)
#'
#' # Species and Stressor distributions
#' dev.off()
#' par(mfrow = c(2,2))
#'
#' plot(st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species 1")
#' plot(risa_maps$species_kernel_maps$sp1$shp, add = TRUE)
#'
#' plot(st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species 2")
#' plot(risa_maps$species_kernel_maps$sp2$shp, add = TRUE)
#'
#' plot(st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Gillnet")
#' plot(risa_maps$stressor_kernel_maps$gillnet$shp, add = TRUE, col=c("lightgreen", "green", "darkgreen"))
#'
#' plot(st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Trawling")
#' plot(risa_maps$stressor_kernel_maps$trawling$shp, add = TRUE, col=c("lightgreen", "green", "darkgreen"))
#'
#' # Overlap maps
#' plot(st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species1 vs Gillnet")
#' plot(risa_maps$overlap_maps$sp1$gillnet, add = TRUE, col=c("yellow", "orange", "red"))
#'
#' plot(st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species2 vs Gillnet")
#' plot(risa_maps$overlap_maps$sp2$gillnet, add = TRUE, col=c("yellow", "orange", "red"))
#'
#' plot(st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species1 vs Trawling")
#' plot(risa_maps$overlap_maps$sp1$gillnet, add = TRUE, col=c("yellow", "orange", "red"))
#'
#' plot(st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species2 vs Trawling")
#' plot(risa_maps$overlap_maps$sp2$gillnet, add = TRUE, col=c("yellow", "orange", "red"))
#' @export
risa_prep <- function(x, y, area = NULL, n_classes = 3,
                      output_layer_type = "shp", radius = NULL,
                      output_decimal_crs = FALSE) {

  # Checking x and y
  spp_list_shp <- spp_df
  if (inherits(spp_list_shp, "data.frame")) {
    message("Input x (species) is a data.frame or shapefile of species distributions")
    spp_list_shp <- df_to_list(spp_list_shp)
  } else if (is.list(spp_list_shp)) {
    message("Input x (species) is a list")
  } else {
    stop("Error: Species (x) input must be one of the following object types: 'sf', 'data.frame', or a list of 'sf'")
  }

  stressor_list_shp <- str_df
  if (inherits(stressor_list_shp, "data.frame")) {
    message("Input y (stressor) is a data.frame or shapefile of stressor distributions")
    stressor_list_shp <- df_to_list(stressor_list_shp)
  } else if (is.list(stressor_list_shp)) {
    message("Input y (stressor) is a list")
  } else {
    stop("Error: Stressor (y) input must be one of the following object types: 'sf', 'data.frame', or a list of 'sf'")
  }

  # Check and define area
  if (is.null(area)) {
    message("No area provided: Area calculated as a bounding box around stressor distributions.")
    merged_stressor_df <- merge_shp(stressor_list_shp)
    area <- create_area(merged_stressor_df, area_type = "convex_hull")
  } else {
    if (inherits(area, "bbox")) {
      message("Area provided is a bounding box. Converting to sf...")
      area <- st_as_sfc(area)
    } else if (inherits(area, "sf")) {
      message("Area provided is a sf object")
    } else if (inherits(area, "data.frame")) {
      message("Area provided is a dataframe of coordinates. Converting to sf...")
      area <- df_to_shp(area)
    } else {
      stop("Provided area must be bbox, sf or data.frame object.")
    }
  }

  # Helper: calculate kernel from list of 'sf' objects
  generate_kernel_list <- function(shp_list,
                                   area,
                                   n_classes,
                                   radius,
                                   output_decimal_crs) {
    lapply(shp_list, function(item) {
      ker <- get_class_kernel(item,
                              area = area,
                              n_classes = n_classes,
                              output_layer_type = "both",
                              radius = radius)
      if (output_decimal_crs) {
        ker$raster <- convert_to_decimal_degrees(ker$raster)
        ker$shp    <- convert_to_decimal_degrees(ker$shp)
      }
      return(ker)
    })
  }

  # Create Kernels
  spp_kernel_list <- generate_kernel_list(spp_list_shp, area, n_classes, radius, output_decimal_crs)
  stressor_kernel_list <- generate_kernel_list(stressor_list_shp, area, n_classes, radius, output_decimal_crs)

  # Fix area projection, if necessary
  if (output_decimal_crs) {
    area <- convert_to_decimal_degrees(area)
  }

  # Create overlap maps
  message("Generating overlap maps.")
  overlap_maps_list <- lapply(spp_kernel_list, function(sp) {
    lapply(stressor_kernel_list, function(st) {
      get_overlap_kernel(sp$raster, st$raster, n_classes=n_classes, output_layer_type=output_layer_type)
    })
  })

  # Output
  all_maps <- list(spp_kernel_list, stressor_kernel_list, overlap_maps_list, area)
  names(all_maps) <- c("species_kernel_maps", "stressor_kernel_maps", "overlap_maps", "area_of_interest")
  return(all_maps)
}
