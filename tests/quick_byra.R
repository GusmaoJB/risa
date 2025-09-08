quick_byra <- function(x, y, criteria, area = NULL,
                       n_classes = 3,
                       n_overlap = NULL,
                       output_layer_type = c("both","shp","raster"),
                       radius = NULL,
                       radius_method = c("nndist","ppl","fixed"),
                       group_x = NULL,
                       group_y = NULL,
                       group_size_x = NULL,
                       group_size_y = NULL,
                       pixel_size = NULL,
                       dimyx = c(512,512),
                       exclude_lowest = TRUE,
                       lowest_prop = 0.05,
                       area_strategy = c("stressor","species","union"),
                       area_type = c("convex_hull","bbox"),
                       area_buffer_frac = 0.5,
                       return_crs = c("metric","4326"),
                       overlap_method = c("product","sum","geom_mean","max"),
                       equation = c("euclidean","multiplicative"),
                       decay = c("none", "linear", "exponential",
                                 "polynomial_2nd", "polynomial_3rd",
                                 "complementary_decay_2nd", "complementary_decay_3rd"),
                       buffer_m = NULL,
                       quiet = TRUE) {

  # Generate Kernel Density and ditribution maps
  input_maps <- risa_prep(x, y,
                          area,
                          n_classes,
                          output_layer_type,
                          radius,
                          radius_method,
                          group_x,
                          group_y,
                          group_size_x,
                          group_size_y,
                          pixel_size,
                          dimyx,
                          exclude_lowest,
                          lowest_prop,
                          area_strategy,
                          area_type,
                          area_buffer_frac,
                          return_crs,
                          overlap_method,
                          quiet)

  # Reshape lists for HRA analysis


  # Perform HRA


  # Prepare outputs


  # Outputs

}
