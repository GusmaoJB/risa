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
risa_maps$species_distributions$species1$raster

#Load example criteria data
path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
df <- read.csv(path)

#Reshape criteria table
crit_list <- criteria_reshape(df)

