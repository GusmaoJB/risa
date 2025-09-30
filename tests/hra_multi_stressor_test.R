############### Test Field ###############
setwd("/home/jojo/Documents/pontal_projects/risa")
library(risa)

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

#Load example data
#path <- "C:/Users/gusma/Documents/research/test_hra/1sp_2stressors.csv"
path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
df <- read.csv(path)

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
res$total_hotspots_reclassified
terra::plot(res$stressor1$Risk_map)
terra::plot(res$stressor2$Risk_map)
terra::plot(res$total_raw)
terra::plot(res$total)


# per-stressor buffers (named)
res4 <- hra2(rast_list, spp_dist, crit_list,
             equation = "euclidean",
             decay = "linear",
             r_max = 3, n_overlap = 2,
             buffer_m = c(stressor1 = 500000, stressor2 = 1000000))


criteria <- criteria_reshape(df)
sample_crit <- criteria[[1]]
crit_names <- unique(sample_crit[is.na(sample_crit$RATING),"ATTRIBUTES"])

input_mapr <- risa_prep(spp_df, str_df, output_layer_type = "raster")
input_mapb <- risa_prep(spp_df, str_df, output_layer_type = "both")

list_depth_base(input_mapr)
list_depth_base(input_mapb)

input_mapr$species_kernel_maps$species1$raster
input_mapb$species_kernel_maps$species1$raster

input_mapr$overlap_maps$species1$stressor1
input_mapb$overlap_maps$species1$stressor1



raster_list <- reshape_risa_maps(input_maps, crit_names)
species_distr <- input_maps$species_distributions



byra_test <- quick_byra(spp_df, str_df, df)
byra_test$summary_stats
terra::plot(byra_test$ecosys_risk_raw)
terra::plot(byra_test$ecosys_risk_classified)

byra_test$species1$stressor1

# 1) What attribute names did HRA want?
crit_names <- unique(sample_crit[is.na(sample_crit$RATING), "ATTRIBUTES"])
crit_names <- trimws(as.character(crit_names))

# 2) What layer names do we actually have from risa_prep?
#    This tries to collect names from any SpatRaster inside input_maps (and sublists).
get_layer_names <- function(x) {
  out <- character()
  if (inherits(x, "SpatRaster")) out <- names(x)
  if (is.list(x)) out <- c(out, unlist(lapply(x, get_layer_names), use.names = FALSE))
  unique(out)
}
available <- get_layer_names(input_maps)

# 3) See what doesn’t match
missing <- setdiff(crit_names, available)
dup_in_criteria <- crit_names[duplicated(crit_names)]

list(crit_names = crit_names, example_available = head(available), missing = missing, dup = dup_in_criteria)


