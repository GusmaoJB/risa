############### Test Field ###############
#library(risa)

#' # Creating test data
#' set.seed(12)
#' spp_df <- rbind(data.frame(long = rnorm(80, 0, 10),
                           #' lat = rnorm(80, 0, 10), species = "species1"),
                #' data.frame(long = rnorm(60, 0, 10),
                           #' lat = rnorm(60, 0, 10), species = "species2"))
#' str_df <- rbind(data.frame(long = rnorm(100, 0, 5),
                           #' lat = rnorm(100, 0, 10), stressor = "stressor1"),
                #' data.frame(long = rnorm(50, 0, 10),
                           #' lat = rnorm(100, 0, 5), stressor = "stressor2"))
#' # Create kernel maps of species and stressor distributions and overlap maps
#' risa_maps <- risa_prep(spp_df, str_df)
#' #Load example data
#' path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
#' df <- read.csv(path)
#' #Reshape criteria table
#' crit_list <- criteria_reshape(df)
#' # Selecting spatially explicit criteria ratings
#' # Note that the rasters in the stressors's list are named after the respective attribute in the criteria table.
#' rast_list <- list(
  #' species1 = list(
    #' stressor1 = list(
      #' intensity = risa_maps$stressor_kernel_maps$stressor1$raster,
      #' `likelihood of interaction`=risa_maps$overlap_maps$species1$stressor1$raster),
    #' stressor2 = list(
      #' intensity = risa_maps$stressor_kernel_maps$stressor2$raster,
      #' `likelihood of interaction`=risa_maps$overlap_maps$species1$stressor2$raster)
  #' ),
  #' species2 = list(
    #' stressor1 = list(
      #' intensity = risa_maps$stressor_kernel_maps$stressor1$raster,
      #' `likelihood of interaction`=risa_maps$overlap_maps$species2$stressor1$raster),
    #' stressor2 = list(
      #' intensity = risa_maps$stressor_kernel_maps$stressor2$raster,
      #' `likelihood of interaction`=risa_maps$overlap_maps$species2$stressor2$raster)
  #' )
#' )
#'
#' # Species' distributions rasters
#' spp_dist <- list(species1 = risa_maps$species_distributions$species1$raster,
                 #' species2 = risa_maps$species_distributions$species2$raster)
#'
#' # Simple example with one species and one stressor
#' test1 <- hra(rast_list[[1]], spp_dist[[1]], crit_list[[1]], equation = "euclidean")
#' terra::plot(test1$total_raw)
#' test1$summary_stats
#'
#' # Now with two species and two stressors
#' many_test <- hra(rast_list, spp_dist, crit_list, equation = "euclidean")
#' many_test$summary_stats
#' terra::plot(many_test$species1$total_raw)
#' terra::plot(many_test$species2$total_raw)





terra::plot(many_test$ecosys_risk_raw)

test_stat <- get_stats(test1)

test_stat_many <- get_stats(many_test)
test_stat_many

total_raw_index <- which(names(test1) == "total_raw")
total_raw_index
stressor_names <- names(test1)[1:(total_raw_index - 1)]
stressor_names

terra::plot(many_test$ecosys_risk_raw/2)

m_jkl <- 2.83

ecosys_risk_classified <- terra::ifel(many_test$ecosys_risk_raw == 0, 0,
                                      terra::ifel(many_test$ecosys_risk_raw < (1/3)*m_jkl*length(crit_list), 1,
                                                  terra::ifel(many_test$ecosys_risk_raw < (2/3)*m_jkl*length(crit_list), 2, 3)))

terra::plot(ecosys_risk_classified)

terra::plot(test1$total_raw)
terra::plot(test2$total_raw)

max(test1$total_raw)
max(test2$total_raw)


