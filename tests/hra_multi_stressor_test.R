############### Test Field ###############
setwd("/home/jojo/Documents/pontal_projects/risa")
library(devtools)

library(risa)
library(sf)
library(terra)

# Loading real data
avist_data <- read.csv("/home/jojo/Documents/pontal_projects/megacost/avistagens_atualizadas.csv")
sg_data <- avist_data[avist_data$sp == "sg", c("Easting", "Northing", "tam_max", "id_grupo", "data", "area_code", "best_transect")]
sg_data$tam_max <- as.numeric(sg_data$tam_max)
sg_data_clean <- aggregate(. ~ data + area_code + best_transect, data=sg_data[,-4], FUN=mean)
sg_data_clean$tam_max <- round(sg_data_clean$tam_max)
sg_data_cleaner <- sg_data_clean[,c("Easting", "Northing", "tam_max")]

# Creating test data
set.seed(12)
spp_df <- rbind(data.frame(long = rnorm(80, 0, 1),
                           lat = rnorm(80, 0, 1), species = "species1"),
                data.frame(long = rnorm(60, 0, 1),
                           lat = rnorm(60, 0, 1), species = "species2"))
str_df <- rbind(data.frame(long = rnorm(100, 0, 0.8),
                           lat = rnorm(100, 0, 1), stressor = "stressor1"),
                data.frame(long = rnorm(50, 0, 1),
                           lat = rnorm(100, 0, 0.7), stressor = "stressor2"))

test_data1 <- spp_df[spp_df$species == "species1",]
test_data2 <- test_data1[1:30,]
test_data3 <- test_data1[1:10,]
test_data4 <- test_data1[1:5,]

# Function to plot KDE and points
plot_kernel_points <- function(data,
                               x_col = "long",
                               y_col = "lat",
                               input_crs = "EPSG:4326",
                               point_col = "red",
                               point_pch = 16,
                               point_cex = 0.8,
                               ...) {

  kde_raster <- get_class_kernel(
    data,
    input_crs = input_crs,
    output_layer_type = "raster",
    ...
  )

  pts <- terra::vect(
    data,
    geom = c(x_col, y_col),
    crs = input_crs
  )

  pts_proj <- terra::project(
    pts,
    terra::crs(kde_raster)
  )

  terra::plot(kde_raster)

  terra::plot(
    pts_proj,
    add = TRUE,
    pch = point_pch,
    col = point_col,
    cex = point_cex
  )

  invisible(list(
    raster = kde_raster,
    points = pts_proj
  ))
}

sp_dat <- spp_df[spp_df$species == "species1",]
st_dat <- str_df[str_df$stressor == "stressor1",]

r1 <- get_class_kernel(sp_dat)
r2 <- get_class_kernel(st_dat, output_layer_type = "raster")

terra::plot(r1)
terra::plot(r2)

plot_kernel_points(sp_dat)
plot_kernel_points(st_dat)

overlap <- get_overlap_kernel(r1, r2, continuous = TRUE)
terra::plot(overlap)

#Load example data
#path <- "C:/Users/gusma/Documents/research/test_hra/1sp_2stressors.csv"
path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
df <- read.csv(path)
crit_list <- criteria_reshape(df)

crit_list <- criteria_reshape(df)
crit_dat <- check_criteria(crit_list[[1]])
stressors <- unique(crit_dat$STRESSOR)
n_overlap <- length(stressors)
n_overlap

# Risa prep
input_maps <- risa_prep(spp_df, str_df)
input_maps_single <- risa_prep(spp_df[,-3], str_df[,-3])

input_maps$species_kernel_maps$species1
input_maps$species_distributions$species1

#export
#risa::export_maps(input_maps, out_dir = "/home/jojo/Documents/pontal_projects/risa_example_maps")

# Create kernel maps of species and stressor distributions and overlap maps
risa_maps <- risa_prep(spp_df, str_df, quiet=FALSE)

risa_maps$species_kernel_maps$species1
risa_maps$species_distributions$species1

# export to compare with invest
# export_maps(risa_maps, "C:/Users/gusma/Documents/research/test_hra/maps")

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
          r_max = 3, n_overlap = 2, quiet = FALSE)

res2 <- hra(rast_list[[1]], spp_dist[[1]], crit_list[[1]])

# per-stressor buffers (named)
res_2_2 <- hra(rast_list, spp_dist, crit_list,
             equation = "euclidean",
             decay = "linear",
             r_max = 3, n_overlap = 2,
             buffer_m = c(stressor1 = 500000, stressor2 = 1000000),
             quiet = FALSE)

head(spp_df)

output <- quick_byra(spp_df, str_df, df, quiet = FALSE)

crit_list <- criteria_reshape(df)
stressors <- unique(crit_list[[1]]$STRESSOR)
n_overlap <- length(stressors)
n_overlap

output$kde_maps$species_distributions$species1$

risaplot(output$kde_maps)
risaplot(output)

input_maps <- risa_prep(spp_df, str_df)
input_maps_single <- risa_prep(spp_df[,-3], str_df[,-3])

test2 <- quick_byra(risa_maps, criteria = df, quiet = FALSE)
test2$kde_maps$species_kernel_maps$species1$raster$Rating
risaplot(test2)

r_stack <- terra::rast(input_maps$species_distributions)

input_maps$species_distributions$species1

raster_list <- reshape_risa_maps(input_maps, crit_names)
species_distr <- input_maps$species_distributions

"area_of_interest" %in% names(byra_test)

terra::plot(byra_test$species1$stressor1$E_criteria)
terra::plot(byra_test$species1$stressor1$C_criteria)

byra_test$summary_stats
terra::plot(byra_test$ecosys_risk_raw)
terra::plot(byra_test$ecosys_risk_classified)



