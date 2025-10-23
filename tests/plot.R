setwd("/home/jojo/Documents/pontal_projects/risa")
library(sf)
library(risa)
library(purrr)
library(tidyr)
library(dplyr)
library(terra)
library(ggplot2)
library(patchwork)

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

#Load example data
#path <- "C:/Users/gusma/Documents/research/test_hra/1sp_2stressors.csv"
path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
df <- read.csv(path)
crit_list <- criteria_reshape(df)

# input tests
input_maps <- risa_prep(spp_df, str_df)
input_maps_single <- risa_prep(spp_df[,-3], str_df[,-3])
input_maps_shp <- risa_prep(spp_df, str_df, output_layer_type = "shp")
input_maps_raster <- risa_prep(spp_df, str_df, output_layer_type = "raster")

rast_list <- list(
  species1 = list(
    stressor1 = list(
      intensity = input_maps$stressor_kernel_maps$stressor1$raster,
      `likelihood of interaction`=input_maps$overlap_maps$species1$stressor1$raster),
    stressor2 = list(
      intensity = input_maps$stressor_kernel_maps$stressor2$raster,
      `likelihood of interaction`=input_maps$overlap_maps$species1$stressor2$raster)
  ),
  species2 = list(
    stressor1 = list(
      intensity = input_maps$stressor_kernel_maps$stressor1$raster,
      `likelihood of interaction`=input_maps$overlap_maps$species2$stressor1$raster),
    stressor2 = list(
      intensity = input_maps$stressor_kernel_maps$stressor2$raster,
      `likelihood of interaction`=input_maps$overlap_maps$species2$stressor2$raster)
  )
)

rast_list_1st <- list(species1 = list(stressor1 = rast_list[[1]][[1]]),
                      species2 = list(stressor1 = rast_list[[2]][[1]]))

spp_dist <- list(species1 = input_maps$species_distributions$species1$raster,
                 species2 = input_maps$species_distributions$species2$raster)

stressor1 <- list(stressor1 = rast_list[[1]][[1]])
crit_sp_1 <- crit_list[[1]]
crit_sp_1 <- crit_sp_1[crit_sp_1$STRESSOR %in% c(NA, "stressor1"),]

crit_list_1_str <- list()

for (i in 1:length(crit_list)){
  targ_sp <- crit_list[[i]]
  crit_list_1_str[[i]] <- targ_sp[targ_sp$STRESSOR %in% c(NA, "stressor1"),]
}

names(crit_list_1_str) <- names(crit_list)

hra_res0 <- hra(stressor1, spp_dist[[1]], crit_sp_1,
                equation = "multiplicative",
                r_max = 3, n_overlap = 2)

hra_res1 <- hra(rast_list[[1]], spp_dist[[1]], crit_list[[1]],
           equation = "multiplicative",
           r_max = 3, n_overlap = 2)

hra_res <- hra(rast_list, spp_dist, crit_list,
                equation = "multiplicative",
                r_max = 3, n_overlap = 2)

hra_res_1str <- hra(rast_list_1st, spp_dist, crit_list_1_str,
                    equation = "multiplicative",
                    r_max = 3, n_overlap = 2)

names(hra_res0)
names(hra_res1)
names(hra_res_1str)
names(hra_res)

hra_input <- hra_res0

risaplot2(input_maps)


# Helpers
# Convert a raster list into a data.frame
rast_list_to_df <- function(r_list, object_name = "raster", group_names = NULL, value_name = "value") {

  if (is.null(group_names)) {
    group_names <- names(r_list)
  }

  # Build a multi-layer SpatRaster with layer names = species
  r_stack <- terra::rast(lapply(group_names, function(r) {
    r <- r_list[[r]][[object_name]]
    r
  }))

  names(r_stack) <- group_names

  # Convert to data frame and reshape
  df <- as.data.frame(r_stack, xy = TRUE) %>%
    pivot_longer(
      cols = -c(x, y),
      names_to = "group",
      values_to = value_name
    )
}

# Convert a list of sf objects into a grouped sf object
# Add a 'group' column to each sf and row-bind them
join_shps <- function(r_list, group_name = "group") {
  sf_list <- imap(r_list, ~ {
    s <- .x$shp
    s[[group_name]] <- .y
    s
  })

  # rbind preserves sf and CRS
  all_sf <- do.call(rbind, sf_list)
  return(all_sf)
}

# Reclassify labels
reclass_labels <- function(x) {
  output <- ifelse(x == 0, "None",
                   ifelse(x == 1, "Low",
                          ifelse(x == 2, "Medium", "High")))
  output <- factor(output, levels = c("None", "Low", "Medium", "High"))
}

# if risaMaps
# Split input into useful lists
sp_list <- input_maps$species_kernel_maps
st_list <- input_maps$stressor_kernel_maps
overlap_list <- input_maps$overlap_maps
species_names <- names(sp_list)
stressor_names <- names(st_list)
aio <- input_maps$area_of_interest

# Check if iput is raster or shp
raster <- "raster" %in% objects(overlap_list[[1]][[1]])

# Prepare data for plots
# if raster
spp_df <- na.omit(rast_list_to_df(sp_list))
names(spp_df)[1:2] <- c("Longitude", "Latitude")
stressor_df <- na.omit(rast_list_to_df(st_list))
names(stressor_df)[1:2] <- c("Longitude", "Latitude")
overlap_df_list <- list()

for (species in species_names) {
  for (stressor in stressor_names) {
    overlap_df_list[[species]] <- rast_list_to_df(overlap_list[[species]])
  }
}

overlap_df <- na.omit(bind_rows(overlap_df_list, .id = "species"))
names(overlap_df)[2:3] <- c("Longitude", "Latitude")

# if raster false
spp_sf <- sp_list[[1]]
if (length(sp_list) > 1) {
  spp_sf <- join_shps(sp_list)
}

stressor_sf <- st_list[[1]]
if (length(st_list) > 1) {
  stressor_sf <- join_shps(st_list)
}

overlap_sfs <- list()
for (species in species_names) {
 overlap_sfs[[species]] <- join_shps(overlap_list[[species]])
 overlap_sfs[[species]][["species"]] <- species
}

all_overlap <- do.call(rbind, overlap_sfs)

# plots if rasters true
gg_spp <- ggplot() + geom_tile(data=spp_df, aes(x=Longitude, y=Latitude, fill = value)) +
  geom_sf(data=aio, fill="transparent", linewidth=0.5) +
  facet_wrap(~ group, scales = "fixed")

gg_stress <- ggplot() + geom_tile(data=stressor_df, aes(x=x, y=y, fill = value)) +
  geom_sf(data=aio, fill="transparent", linewidth=0.5) +
  facet_wrap(~ group, scales = "fixed")

gg_overlap <- ggplot() + geom_tile(data=overlap_df, aes(x=x, y=y, fill = value)) +
  geom_sf(data=aio, fill="transparent", linewidth=0.5) +
  facet_grid(group ~ species, scales = "fixed")

# plots if raster false
gg_spp <- ggplot() + geom_sf(data=spp_sf, aes(fill=Rating), col="transparent") +
  geom_sf(data=aio, fill="transparent", linewidth=0.5) +
  facet_wrap(~ group, scales = "fixed")

gg_stress <- ggplot() + geom_sf(data=stressor_sf, aes(fill=Rating), col="transparent") +
  geom_sf(data=aio, fill="transparent", linewidth=0.5) +
  facet_wrap(~ group, scales = "fixed")

gg_overlap <- ggplot() + geom_sf(data=all_overlap, aes(fill=Rating), col="transparent") +
  geom_sf(data=aio, fill="transparent", linewidth=0.5) +
  facet_grid(group ~ species, scales = "fixed")

# output for risaMaps
gg_kde_list <- list(gg_spp + ggtitle("Species KDE"),
                gg_stress + ggtitle("Stressor KDE"),
                gg_overlap + ggtitle("Species-Stressor Overlaps"))

gg_kde_list

# If risaHRA
depth <- list_depth_base(hra_input)
str_vecs <- unique(hra_input$summary_stats$STRESSOR)
stressor_names <- str_vecs[!str_vecs %in% "(FROM ALL STRESSORS)"]
hra_input$summary_stats
spp_vecs <- unique(hra_input$summary_stats$SPECIES)
species_names <- spp_vecs[!spp_vecs %in% "ECOSYSTEM"]

# if (depth == 2)
risk_raw <- as.data.frame(hra_input$total_raw, xy=TRUE)
names(risk_raw) <- c("Longitude", "Latitude", "Risk")
risk_reclassified <- as.data.frame(hra_input$total_hotspots_reclassified, xy = TRUE)
names(risk_reclassified) <- c("Longitude", "Latitude", "Risk (reclass.)")
risk_reclassified$`Risk (reclass.)` <- reclass_labels(risk_reclassified$`Risk (reclass.)`)

# if length(stressor_names == 1)
gg_risk_raw <- ggplot() + geom_tile(data=risk_raw, aes(x=Longitude, y=Latitude, fill = Risk))
gg_risk_reclass <- ggplot() + geom_tile(data=risk_reclassified, aes(x=Longitude, y=Latitude, fill = `Risk (reclass.)`))

gg_risk_1sp_1str <- list(gg_risk_raw + ggtitle("Total risk"),
                         gg_risk_reclass + ggtitle("Total risk (reclassified)"))

gg_risk_1sp_1str # output if 1 species and 1 stressor

# if length(stressor_names > 1)
stressor_risk_raw <- rast_list_to_df(hra_input, "Risk_map_raw", group_names = stressor_names, value_name = "Risk")
names(stressor_risk_raw)[1:2] <- c("Longitude", "Latitude")
stressor_risk_reclassified <- rast_list_to_df(hra_input, "Risk_map", group_names = stressor_names, value_name = "Risk (reclass.)")
names(stressor_risk_reclassified)[1:2] <- c("Longitude", "Latitude")
stressor_risk_reclassified$`Risk (reclass.)` <- reclass_labels(stressor_risk_reclassified$`Risk (reclass.)`)
risk_max_ratings <- as.data.frame(hra_input$total_reclassified, xy=TRUE)
names(risk_max_ratings) <- c("Longitude", "Latitude", "Highest str. risk (reclass.)")
risk_max_ratings$`Highest str. risk (reclass.)` <- reclass_labels(risk_max_ratings$`Highest str. risk (reclass.)`)

gg_str_risk_raw <- ggplot() + geom_tile(data=stressor_risk_raw, aes(x=Longitude, y=Latitude, fill = Risk)) +
  facet_wrap(~ group, scales = "fixed")

gg_str_risk_reclass <- ggplot() + geom_tile(data=stressor_risk_reclassified, aes(x=Longitude, y=Latitude, fill = `Risk (reclass.)`)) +
  facet_wrap(~ group, scales = "fixed")

gg_str_risk_max_risk <- ggplot() + geom_tile(data=risk_max_ratings, aes(x=Longitude, y=Latitude, fill=`Highest str. risk (reclass.)`))

gg_risk_raw <- ggplot() + geom_tile(data=risk_raw, aes(x=Longitude, y=Latitude, fill = Risk))

gg_risk_reclass <- ggplot() + geom_tile(data=risk_reclassified, aes(x=Longitude, y=Latitude, fill = `Risk (reclass.)`))

gg_risk_list <- list(gg_str_risk_raw + ggtitle("Stressor risk"),
                     gg_str_risk_reclass + ggtitle("Reclassified stressor risk"),
                     gg_str_risk_max_risk + ggtitle("Highest stressor risk estimates"),
                     gg_risk_raw + ggtitle("Total (combined) risk"),
                     gg_risk_reclass + ggtitle("Total (combined) risk (reclassified)"))

gg_risk_list # output if 1 species and multiple stressors

# if (depth == 3)
hra_input <- hra_res_1str

eco_risk_raw <- as.data.frame(hra_input$ecosys_risk_raw, xy=TRUE)
names(eco_risk_raw) <- c("Longitude", "Latitude", "Risk")
eco_risk_reclassified <- as.data.frame(hra_input$ecosys_risk_classified, xy = TRUE)
names(eco_risk_reclassified) <- c("Longitude", "Latitude", "Risk (reclass.)")
eco_risk_reclassified$`Risk (reclass.)` <- reclass_labels(eco_risk_reclassified$`Risk (reclass.)`)

gg_ecosys_raw <- ggplot() + geom_tile(data=eco_risk_raw, aes(x=Longitude, y=Latitude, fill = Risk))
gg_ecosys_reclass <- ggplot() + geom_tile(data=eco_risk_reclassified, aes(x=Longitude, y=Latitude, fill = `Risk (reclass.)`))

species_names <- c("species1", "species2")
raw_risk_df <- data.frame()
for (sp in species_names) {
  sp_df <- as.data.frame(hra_input[[sp]]$total_raw, xy=TRUE)
  sp_df$species <- sp
  raw_risk_df <- rbind.data.frame(raw_risk_df, sp_df)
}
names(raw_risk_df) <- c("Longitude", "Latitude", "Risk", "species")

recl_risk_df <- data.frame()
for (sp in species_names) {
  sp_df <- as.data.frame(hra_input[[sp]]$total_reclassified, xy=TRUE)
  sp_df$species <- sp
  recl_risk_df <- rbind.data.frame(recl_risk_df, sp_df)
}
names(recl_risk_df) <- c("Longitude", "Latitude", "Risk (reclass.)", "species")
recl_risk_df$`Risk (reclass.)` <- reclass_labels(recl_risk_df$`Risk (reclass.)`)

gg_raw_risk <- ggplot() + geom_tile(data=raw_risk_df, aes(x=Longitude, y=Latitude, fill=Risk)) + facet_grid(.~species)
gg_recl_risk <- ggplot() + geom_tile(data=recl_risk_df, aes(x=Longitude, y=Latitude, fill=`Risk (reclass.)`)) + facet_grid(.~species)

# if length(stressor_names == 1)
gg_spp_1str <- list(gg_raw_risk + ggtitle("Stressor risk"),
                    gg_recl_risk + ggtitle("Reclassified stressor risk"),
                    gg_ecosys_raw + ggtitle("Ecosystem risk"),
                    gg_ecosys_reclass + ggtitle("Reclassified ecosystem risk"))

# if length(stressor_names > 1)
raw_str_risk_df <- data.frame()
for (sp in species_names) {
  for (str in stressor_names) {
    spec <- as.data.frame(hra_input[[sp]][[str]]$Risk_map_raw, xy=TRUE)
    spec$species <- sp
    spec$stressor <- str
    raw_str_risk_df <- rbind.data.frame(raw_str_risk_df, spec)
  }
}

names(raw_str_risk_df) <- c("Longitude", "Latitude", "Risk", "species", "stressor")

reclass_str_risk_df <- data.frame()
for (sp in species_names) {
  for (str in stressor_names) {
    spec <- as.data.frame(hra_input[[sp]][[str]]$Risk_map, xy=TRUE)
    spec$species <- sp
    spec$stressor <- str
    reclass_str_risk_df <- rbind.data.frame(reclass_str_risk_df, spec)
  }
}

names(reclass_str_risk_df) <- c("Longitude", "Latitude", "Risk (reclass.)", "species", "stressor")
reclass_str_risk_df$`Risk (reclass.)` <- reclass_labels(reclass_str_risk_df$`Risk (reclass.)`)

gg_str_raw_risk <- ggplot() + geom_tile(data=raw_str_risk_df, aes(x=Longitude, y=Latitude, fill=Risk)) +
  facet_grid(stressor ~ species)

gg_str_recl_risk <- ggplot() + geom_tile(data=reclass_str_risk_df, aes(x=Longitude, y=Latitude, fill=`Risk (reclass.)`)) +
  facet_grid(stressor ~ species)

gg_spp_strss <- list(gg_str_raw_risk + ggtitle("Stressor risk"),
                     gg_str_recl_risk + ggtitle("Reclassified stressor risk"),
                     gg_raw_risk + ggtitle("Total (combined stress.) risk"),
                     gg_recl_risk + ggtitle("Reclassified total (combined stress.) risk"),
                     gg_ecosys_raw + ggtitle("Ecosystem risk"),
                     gg_ecosys_reclass + ggtitle("Reclassified ecosystem risk"))

# if ("area_of_interest" %in% names(hra_input))
hra_input <- byra_test

aio <- geom_sf(data=hra_input$area_of_interest, fill="transparent", linewidth=0.5)

gg_spp_strss <- list(gg_str_raw_risk + aio + ggtitle("Stressor risk"),
                     gg_str_recl_risk + aio + ggtitle("Reclassified stressor risk"),
                     gg_raw_risk + aio + ggtitle("Total (combined stress.) risk"),
                     gg_recl_risk + aio + ggtitle("Reclassified total (combined stress.) risk"),
                     gg_ecosys_raw + aio + ggtitle("Ecosystem risk"),
                     gg_ecosys_reclass + aio + ggtitle("Reclassified ecosystem risk"))


