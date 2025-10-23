library(sf)
library(risa)
library(purrr)
library(tidyr)
library(dplyr)
library(terra)
library(ggplot2)
library(patchwork)

risaplot <- function(x) {

  # Checks
  stopifnot(inherits(x, c("risaMaps", "risaHRA")))

  # Helpers
  # Convert a raster list into a data.frame
  .rast_list_to_df <- function(r_list, object_name = "raster", group_names = NULL, value_name = "value") {

    if (is.null(group_names)) {
      group_names <- names(r_list)
    }

    r_stack <- terra::rast(lapply(group_names, function(r) {
      r <- r_list[[r]][[object_name]]
      r
    }))

    names(r_stack) <- group_names

    df <- as.data.frame(r_stack, xy = TRUE) %>%
      pivot_longer(
        cols = -c(x, y),
        names_to = "group",
        values_to = value_name
      )
  }

  # Convert a list of sf objects into a grouped sf object
  .join_shps <- function(r_list, group_name = "group") {
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
  .reclass_labels <- function(x) {
    output <- ifelse(x == 0, "None",
                     ifelse(x == 1, "Low",
                            ifelse(x == 2, "Medium", "High")))
    output <- factor(output, levels = c("None", "Low", "Medium", "High"))
  }

  if (inherits(x, "risaMaps")) {
    # Split input into useful objects
    sp_list <- x$species_kernel_maps
    st_list <- x$stressor_kernel_maps
    overlap_list <- x$overlap_maps
    species_names <- names(sp_list)
    stressor_names <- names(st_list)
    aio <- x$area_of_interest

    # Check if iput is raster or shp
    raster <- "raster" %in% objects(overlap_list[[1]][[1]])

    if (raster) {
      spp_df <- na.omit(.rast_list_to_df(sp_list))
      names(spp_df)[1:2] <- c("Longitude", "Latitude")
      stressor_df <- na.omit(.rast_list_to_df(st_list))
      names(stressor_df)[1:2] <- c("Longitude", "Latitude")
      overlap_df_list <- list()

      for (species in species_names) {
        for (stressor in stressor_names) {
          overlap_df_list[[species]] <- .rast_list_to_df(overlap_list[[species]])
        }
      }

      overlap_df <- na.omit(bind_rows(overlap_df_list, .id = "species"))
      names(overlap_df)[2:3] <- c("Longitude", "Latitude")

      gg_spp <- ggplot() + geom_tile(data=spp_df, aes(x=Longitude, y=Latitude, fill = value)) +
        geom_sf(data=aio, fill="transparent", linewidth=0.5) +
        facet_wrap(~ group, scales = "fixed")

      gg_stress <- ggplot() + geom_tile(data=stressor_df, aes(x=Longitude, y=Latitude, fill = value)) +
        geom_sf(data=aio, fill="transparent", linewidth=0.5) +
        facet_wrap(~ group, scales = "fixed")

      gg_overlap <- ggplot() + geom_tile(data=overlap_df, aes(x=Longitude, y=Latitude, fill = value)) +
        geom_sf(data=aio, fill="transparent", linewidth=0.5) +
        facet_grid(group ~ species, scales = "fixed")

      gg_kde_list <- list(gg_spp + ggtitle("Species KDE"),
                          gg_stress + ggtitle("Stressor KDE"),
                          gg_overlap + ggtitle("Species-Stressor Overlaps"))

      return(gg_kde_list)

    } else {
      spp_sf <- sp_list[[1]]
      if (length(sp_list) > 1) {
        spp_sf <- .join_shps(sp_list)
      }

      stressor_sf <- st_list[[1]]
      if (length(st_list) > 1) {
        stressor_sf <- .join_shps(st_list)
      }

      overlap_sfs <- list()
      for (species in species_names) {
        overlap_sfs[[species]] <- .join_shps(overlap_list[[species]])
        overlap_sfs[[species]][["species"]] <- species
      }

      all_overlap <- do.call(rbind, overlap_sfs)

      gg_spp <- ggplot() + geom_sf(data=spp_sf, aes(fill=Rating), col="transparent") +
        geom_sf(data=aio, fill="transparent", linewidth=0.5) +
        facet_wrap(~ group, scales = "fixed")

      gg_stress <- ggplot() + geom_sf(data=stressor_sf, aes(fill=Rating), col="transparent") +
        geom_sf(data=aio, fill="transparent", linewidth=0.5) +
        facet_wrap(~ group, scales = "fixed")

      gg_overlap <- ggplot() + geom_sf(data=all_overlap, aes(fill=Rating), col="transparent") +
        geom_sf(data=aio, fill="transparent", linewidth=0.5) +
        facet_grid(group ~ species, scales = "fixed")

      gg_kde_list <- list(gg_spp + ggtitle("Species KDE"),
                          gg_stress + ggtitle("Stressor KDE"),
                          gg_overlap + ggtitle("Species-Stressor Overlaps"))

      return(gg_kde_list)
    }

  } else {
    # Split input into useful objects
    depth <- list_depth_base(x)
    str_vecs <- unique(x$summary_stats$STRESSOR)
    stressor_names <- str_vecs[!str_vecs %in% "(FROM ALL STRESSORS)"]
    x$summary_stats
    spp_vecs <- unique(x$summary_stats$SPECIES)
    species_names <- spp_vecs[!spp_vecs %in% "ECOSYSTEM"]
    area <- x$area_of_interest
    aio <- NULL
    if ("area_of_interest" %in% names(x)) {
      aio <- geom_sf(data=x$area_of_interest, fill="transparent", linewidth=0.5)
    }

    # Check if there are multiple species
    if (depth == 2) {
      risk_raw <- as.data.frame(x$total_raw, xy=TRUE)
      names(risk_raw) <- c("Longitude", "Latitude", "Risk")
      risk_reclassified <- as.data.frame(x$total_hotspots_reclassified, xy = TRUE)
      names(risk_reclassified) <- c("Longitude", "Latitude", "Risk (reclass.)")
      risk_reclassified$`Risk (reclass.)` <- .reclass_labels(risk_reclassified$`Risk (reclass.)`)

      if (is.null(aio)) {
        crs_r <- st_crs(x$total_raw)
        bb <- st_bbox(c(xmin = min(risk_raw$Longitude, na.rm = TRUE),
                        ymin = min(risk_raw$Latitude, na.rm = TRUE),
                        xmax = max(risk_raw$Longitude, na.rm = TRUE),
                        ymax = max(risk_raw$Latitude, na.rm = TRUE)),
                      crs = crs_r)
        bb_sfc <- st_as_sfc(bb)
        bb_sf <- st_sf(bb_sfc)
        aio <- geom_sf(data=bb_sf, fill="transparent", col="transparent", linewidth=0.5)
      }

      # Check if multiple stressors
      if (length(stressor_names) == 1) {
        gg_risk_raw <- ggplot() + geom_tile(data=risk_raw, aes(x=Longitude, y=Latitude, fill = Risk))
        gg_risk_reclass <- ggplot() + geom_tile(data=risk_reclassified, aes(x=Longitude, y=Latitude, fill = `Risk (reclass.)`))

        gg_risk_1sp_1str <- list(gg_risk_raw + ggtitle("Total risk") + aio,
                                 gg_risk_reclass + ggtitle("Total risk (reclassified)") + aio)

        return(gg_risk_1sp_1str)
      } else {
        stressor_risk_raw <- .rast_list_to_df(x, "Risk_map_raw", group_names = stressor_names, value_name = "Risk")
        names(stressor_risk_raw)[1:2] <- c("Longitude", "Latitude")
        stressor_risk_reclassified <- .rast_list_to_df(x, "Risk_map", group_names = stressor_names, value_name = "Risk (reclass.)")
        names(stressor_risk_reclassified)[1:2] <- c("Longitude", "Latitude")
        stressor_risk_reclassified$`Risk (reclass.)` <- .reclass_labels(stressor_risk_reclassified$`Risk (reclass.)`)
        risk_max_ratings <- as.data.frame(x$total_reclassified, xy=TRUE)
        names(risk_max_ratings) <- c("Longitude", "Latitude", "Highest str. risk (reclass.)")
        risk_max_ratings$`Highest str. risk (reclass.)` <- .reclass_labels(risk_max_ratings$`Highest str. risk (reclass.)`)

        gg_str_risk_raw <- ggplot() + geom_tile(data=stressor_risk_raw, aes(x=Longitude, y=Latitude, fill = Risk)) +
          facet_wrap(~ group, scales = "fixed")

        gg_str_risk_reclass <- ggplot() + geom_tile(data=stressor_risk_reclassified, aes(x=Longitude, y=Latitude, fill = `Risk (reclass.)`)) +
          facet_wrap(~ group, scales = "fixed")

        gg_str_risk_max_risk <- ggplot() + geom_tile(data=risk_max_ratings, aes(x=Longitude, y=Latitude, fill=`Highest str. risk (reclass.)`))

        gg_risk_raw <- ggplot() + geom_tile(data=risk_raw, aes(x=Longitude, y=Latitude, fill = Risk))

        gg_risk_reclass <- ggplot() + geom_tile(data=risk_reclassified, aes(x=Longitude, y=Latitude, fill = `Risk (reclass.)`))

        gg_risk_list <- list(gg_str_risk_raw + ggtitle("Stressor risk") + aio,
                             gg_str_risk_reclass + ggtitle("Reclassified stressor risk") + aio,
                             gg_str_risk_max_risk + ggtitle("Highest stressor risk estimates") + aio,
                             gg_risk_raw + ggtitle("Total (combined) risk") + aio,
                             gg_risk_reclass + ggtitle("Total (combined) risk (reclassified)") + aio)
        return(gg_risk_list)
      }
    } else {
      eco_risk_raw <- as.data.frame(x$ecosys_risk_raw, xy=TRUE)
      names(eco_risk_raw) <- c("Longitude", "Latitude", "Risk")
      eco_risk_reclassified <- as.data.frame(x$ecosys_risk_classified, xy = TRUE)
      names(eco_risk_reclassified) <- c("Longitude", "Latitude", "Risk (reclass.)")
      eco_risk_reclassified$`Risk (reclass.)` <- .reclass_labels(eco_risk_reclassified$`Risk (reclass.)`)

      if (is.null(aio)) {
        crs_r <- st_crs(x$ecosys_risk_raw)
        bb <- st_bbox(c(xmin = min(eco_risk_raw$Longitude, na.rm = TRUE),
                        ymin = min(eco_risk_raw$Latitude, na.rm = TRUE),
                        xmax = max(eco_risk_raw$Longitude, na.rm = TRUE),
                        ymax = max(eco_risk_raw$Latitude, na.rm = TRUE)),
                      crs = crs_r)
        bb_sfc <- st_as_sfc(bb)
        bb_sf <- st_sf(bb_sfc)
        aio <- geom_sf(data=bb_sf, fill="transparent", col="transparent", linewidth=0.5)
      }

      gg_ecosys_raw <- ggplot() + geom_tile(data=eco_risk_raw, aes(x=Longitude, y=Latitude, fill = Risk))
      gg_ecosys_reclass <- ggplot() + geom_tile(data=eco_risk_reclassified, aes(x=Longitude, y=Latitude, fill = `Risk (reclass.)`))

      species_names <- c("species1", "species2")
      raw_risk_df <- data.frame()
      for (sp in species_names) {
        sp_df <- as.data.frame(x[[sp]]$total_raw, xy=TRUE)
        sp_df$species <- sp
        raw_risk_df <- rbind.data.frame(raw_risk_df, sp_df)
      }
      names(raw_risk_df) <- c("Longitude", "Latitude", "Risk", "species")

      recl_risk_df <- data.frame()
      for (sp in species_names) {
        sp_df <- as.data.frame(x[[sp]]$total_reclassified, xy=TRUE)
        sp_df$species <- sp
        recl_risk_df <- rbind.data.frame(recl_risk_df, sp_df)
      }
      names(recl_risk_df) <- c("Longitude", "Latitude", "Risk (reclass.)", "species")
      recl_risk_df$`Risk (reclass.)` <- .reclass_labels(recl_risk_df$`Risk (reclass.)`)

      gg_raw_risk <- ggplot() + geom_tile(data=raw_risk_df, aes(x=Longitude, y=Latitude, fill=Risk)) + facet_grid(.~species)
      gg_recl_risk <- ggplot() + geom_tile(data=recl_risk_df, aes(x=Longitude, y=Latitude, fill=`Risk (reclass.)`)) + facet_grid(.~species)

      # Check if multiple stressors
      if (length(stressor_names) == 1) {
        gg_spp_1str <- list(gg_raw_risk + ggtitle("Stressor risk") + aio,
                            gg_recl_risk + ggtitle("Reclassified stressor risk") + aio,
                            gg_ecosys_raw + ggtitle("Ecosystem risk") + aio,
                            gg_ecosys_reclass + ggtitle("Reclassified ecosystem risk") + aio)
        return(gg_spp_1str)
      } else {
        raw_str_risk_df <- data.frame()
        for (sp in species_names) {
          for (str in stressor_names) {
            spec <- as.data.frame(x[[sp]][[str]]$Risk_map_raw, xy=TRUE)
            spec$species <- sp
            spec$stressor <- str
            raw_str_risk_df <- rbind.data.frame(raw_str_risk_df, spec)
          }
        }

        names(raw_str_risk_df) <- c("Longitude", "Latitude", "Risk", "species", "stressor")

        reclass_str_risk_df <- data.frame()
        for (sp in species_names) {
          for (str in stressor_names) {
            spec <- as.data.frame(x[[sp]][[str]]$Risk_map, xy=TRUE)
            spec$species <- sp
            spec$stressor <- str
            reclass_str_risk_df <- rbind.data.frame(reclass_str_risk_df, spec)
          }
        }

        names(reclass_str_risk_df) <- c("Longitude", "Latitude", "Risk (reclass.)", "species", "stressor")
        reclass_str_risk_df$`Risk (reclass.)` <- .reclass_labels(reclass_str_risk_df$`Risk (reclass.)`)

        gg_str_raw_risk <- ggplot() + geom_tile(data=raw_str_risk_df, aes(x=Longitude, y=Latitude, fill=Risk)) +
          facet_grid(stressor ~ species)

        gg_str_recl_risk <- ggplot() + geom_tile(data=reclass_str_risk_df, aes(x=Longitude, y=Latitude, fill=`Risk (reclass.)`)) +
          facet_grid(stressor ~ species)

        gg_spp_strss <- list(gg_str_raw_risk + ggtitle("Stressor risk") + aio,
                             gg_str_recl_risk + ggtitle("Reclassified stressor risk") + aio,
                             gg_raw_risk + ggtitle("Total (combined stress.) risk"),
                             gg_recl_risk + ggtitle("Reclassified total (combined stress.) risk") + aio,
                             gg_ecosys_raw + ggtitle("Ecosystem risk") + aio,
                             gg_ecosys_reclass + ggtitle("Reclassified ecosystem risk") + aio)

        return(gg_spp_strss)
      }
    }
  }

}


