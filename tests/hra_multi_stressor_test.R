############### Test Field ###############

# Creating test data
set.seed(12)
spp_df <- rbind(data.frame(long = rnorm(80, 0, 10),
                           lat = rnorm(80, 0, 10), species = "species1"),
                data.frame(long = rnorm(60, 0, 10),
                           lat = rnorm(60, 0, 10), species = "species2"))
str_df <- rbind(data.frame(long = rnorm(100, 0, 5),
                           lat = rnorm(100, 0, 10), stressor = "stressor1"),
                data.frame(long = rnorm(50, 0, 10),
                           lat = rnorm(100, 0, 5), stressor = "stressor2"))
# Create kernel maps of species and stressor distributions and overlap maps
risa_maps <- risa_prep(spp_df, str_df)
# export to compare with invest
# export_maps(risa_maps, "C:/Users/gusma/Documents/research/test_hra/maps")
getwd()
#Load example data
#path <- "C:/Users/gusma/Documents/research/test_hra/1sp_2stressors.csv"
path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
df <- read.csv(path)
#Reshape criteria table
crit_list <- criteria_reshape(df)
# Selecting spatially explicit criteria ratings
# Note that the rasters in the stressors's list are named after the respective attribute in the criteria table.
rast_list <- list(
   species1 = list(
     stressor1 = list(
       intensity = risa_maps$stressor_kernel_maps$stressor1$raster,
       `likelihood of interaction`=risa_maps$overlap_maps$species1$stressor1$raster),
     stressor2 = list(
       intensity = risa_maps$stressor_kernel_maps$stressor2$raster,
       `likelihood of interaction`=risa_maps$overlap_maps$species1$stressor2$raster)
   ),
   species2 = list(
     stressor1 = list(
       intensity = risa_maps$stressor_kernel_maps$stressor1$raster,
       `likelihood of interaction`=risa_maps$overlap_maps$species2$stressor1$raster),
     stressor2 = list(
       intensity = risa_maps$stressor_kernel_maps$stressor2$raster,
       `likelihood of interaction`=risa_maps$overlap_maps$species2$stressor2$raster)
   )
 )

# Species' distributions rasters
spp_dist <- list(species1 = risa_maps$species_distributions$species1$raster,
                 species2 = risa_maps$species_distributions$species2$raster)



# single value for all stressors, linear decay within 10 km
res <- hra(rast_list[[1]], spp_dist[[1]], crit_list[[1]],
          equation = "multiplicative",
          r_max = 3, n_overlap = 2)

terra::plot(res$stressor1$Risk_map)
terra::plot(res$stressor2$Risk_map)
terra::plot(res$total_raw)
terra::plot(res$total)


# per-stressor buffers (named)
res4 <- hra7(rast_list, spp_dist, crit_list,
             equation = "multiplicative",
             decay = "linear",
             r_max = 3, n_overlap = 2,
             buffer_m = c(stressor1 = 500000, stressor2 = 1000000))

terra::plot(res4$species1$stressor1$Risk_map)
terra::plot(res4$species1$stressor2$Risk_map)
terra::plot(res4$species1$total)
terra::plot(res4$species2$total)

terra::plot(res4$ecosys_risk_classified)

res3$summary_stats


terra::plot(rast_list$species1$stressor1$intensity)
terra::plot(res3$species1$stressor1$Risk_map_raw)

par(mfrow=c(1,2))
terra::plot(res3$species1$stressor1$Risk_map)
terra::plot(res3$species1$stressor2$Risk_map)

test <- terra::mosaic(res3$species1$stressor1$Risk_map,
                      res3$species1$stressor2$Risk_map, fun=max)

terra::plot(test)
terra::plot(res3$species2$stressor1$Risk_map_raw)
terra::plot(res3$species2$stressor2$Risk_map)

terra::plot(res3$species2$total_raw)


terra::plot(res3$species1$total)
terra::plot(res3$species2$total)

terra::plot(res3$ecosys_risk_raw)
terra::plot(res3$ecosys_risk_classified)

res3$summary_stats[-1,]



export_maps(res3$ecosys_risk_raw, "C:/Users/gusma/Documents/research/test_hra/curius_output13")


test <- decay_coeffs(rast_list[[1]][[1]][[1]],
                     raster2 = spp_dist[[1]],
                         decay = "linear",
                         500000)
terra::plot(test)

test2 <- get_decay_map(rast_list[[1]][[1]][[1]], test, start_value = 1)
terra::plot(test2)

general_decay <- c(1:10)
general_decay <- 1

general_decay != 1
is.integer(general_decay)



buffer_m <- c(stressor1 = 500000, stressor2 = 1000000)

# What attributes do your criteria expect?
sort(unique(trimws(crit_list[[1]]$ATTRIBUTES)))

# What does your container actually provide for a given stressor?
names(rast_list[[1]]$stressor1)             # if it's a list of rasters, names()
names(rast_list[[1]][["stressor1"]])        # idem
names(rast_list[[1]][["stressor1"]])        # if SpatRaster stack, layer names

# Do buffer names match stressor names exactly?
names(buffer_m); names(rast_list[[1]])


