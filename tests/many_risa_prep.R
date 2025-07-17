# This function will be modified to incorporate interest areas with multiple subareas
many_risa_prep <- function(x, y, area = NULL, n_classes = 3,
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
