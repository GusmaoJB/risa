risaplot2 <- function(x, area_of_interest = TRUE, raw_ggplot=FALSE) {

  stopifnot(inherits(x, c("risaMaps","risaHRA")))

  # Helpers
  # Reclassify risk caterogy labels
  .reclass_labels <- function(v) {
    labs <- c("0"="None","1"="Low","2"="Medium","3"="High")
    fct <- factor(unname(labs[as.character(v)]),
                  levels = c("None","Low","Medium","High"))
    fct
  }

  # Define aoi and its ggplot2 layer
  .ensure_aoi_layer <- function(aoi, xy_df = NULL, crs_like = NULL) {
    if (inherits(aoi, c("sf","sfc"))) {
      return(ggplot2::geom_sf(data = aoi, fill = "transparent", linewidth = 0.5))
    }

    stopifnot(!is.null(xy_df))
    stopifnot(!is.null(crs_like))
    bb <- st_bbox(
      c(xmin = min(xy_df$Longitude, na.rm = TRUE),
        ymin = min(xy_df$Latitude, na.rm = TRUE),
        xmax = max(xy_df$Longitude, na.rm = TRUE),
        ymax = max(xy_df$Latitude, na.rm = TRUE)),
      crs = st_crs(crs_like)
    )
    ggplot2::geom_sf(data = st_as_sf(st_as_sfc(bb)),
            fill = "transparent", col = "transparent")
  }

  # Convert a raster list into a data.frame
  # Extract raster from either: a direct SpatRaster or a list containing a $raster element
  .get_raster <- function(x, object_name = "raster") {
    if (inherits(x, "SpatRaster")) {
      return(x)
    }

    if (is.list(x) && !is.null(x[[object_name]]) &&
        inherits(x[[object_name]], "SpatRaster")) {
      return(x[[object_name]])
    }

    stop("Could not find a SpatRaster in this object.")
  }

  # Convert a raster list into a data.frame
  .rast_list_to_df <- function(r_list, object_name = "raster",
                               group_names = NULL, value_name = "value") {

    if (is.null(group_names)) group_names <- names(r_list)

    r_stack <- terra::rast(
      purrr::map(group_names, ~ .get_raster(r_list[[.x]], object_name))
    )

    names(r_stack) <- group_names

    df <- terra::as.data.frame(r_stack, xy = TRUE)

    names(df)[1:2] <- c("Longitude", "Latitude")

    df |>
      tidyr::pivot_longer(
        -c(Longitude, Latitude),
        names_to = "group",
        values_to = value_name
      ) |>
      dplyr::filter(!is.na(.data[[value_name]]))
  }

  #Merge multiple sf objects
  .join_shps <- function(r_list, group_col = "group") {
    # bind_rows is fast and preserves sf; add group from names
    imap(r_list, ~ dplyr::mutate(.x$shp, !!group_col := .y)) |>
      vctrs::vec_rbind()
  }

  # Define fill scales for each plot type
  .fill_scale <- function(type, fill_col = NULL) {

    is_reclass <- !is.null(fill_col) &&
      grepl("reclass", fill_col, ignore.case = TRUE)

    risk_cols <- c(
      "None"   = hcl.colors(4, "YlOrRd")[4],
      "Low"    = hcl.colors(4, "YlOrRd")[3],
      "Medium" = hcl.colors(4, "YlOrRd")[2],
      "High"   = hcl.colors(4, "YlOrRd")[1]
    )

    switch(
      type,

      species = scale_fill_gradientn(
        colors = hcl.colors(120, "Greens", rev=TRUE)[31:120],
        na.value = "transparent"
      ),

      stressor = scale_fill_gradientn(
        colors = hcl.colors(120, "Blues", rev=TRUE)[31:120],
        na.value = "transparent"
      ),

      overlap = scale_fill_gradientn(
        colors = hcl.colors(100, "Viridis", rev=TRUE),
        na.value = "transparent"
      ),

      risk = {
        if (is_reclass) {
          scale_fill_manual(
            values = risk_cols,
            limits = c("None", "Low", "Medium", "High"),
            breaks = c("None", "Low", "Medium", "High"),
            drop = FALSE,
            na.value = "transparent",
            guide = guide_legend(
              override.aes = list(
                fill = unname(risk_cols),
                alpha = 1
              )
            )
          )
        } else {
          scale_fill_gradientn(
            colors = hcl.colors(120, "YlOrRd", rev = TRUE)[21:120],
            na.value = "transparent"
          )
        }
      },

      stop("Unknown palette type.")
    )
  }

  # Discrete version
  .fill_scale_discrete <- function(type) {

    risk_levels <- c("None", "Low", "Medium", "High")

    risk_cols <- c(
      "None"   = hcl.colors(4, "YlOrRd")[4],
      "Low"    = hcl.colors(4, "YlOrRd")[3],
      "Medium" = hcl.colors(4, "YlOrRd")[2],
      "High"   = hcl.colors(4, "YlOrRd")[1]
    )

    switch(
      type,

      species = scale_fill_brewer(
        palette = "Greens",
        drop = FALSE,
        na.value = "transparent"
      ),

      stressor = scale_fill_brewer(
        palette = "Blues",
        drop = FALSE,
        na.value = "transparent"
      ),

      overlap = scale_fill_viridis_d(
        option = "viridis",
        drop = FALSE,
        na.value = "transparent"
      ),

      risk = scale_fill_manual(
        values = risk_cols,
        limits = risk_levels,
        breaks = risk_levels,
        drop = FALSE,
        na.value = "transparent"
      ),

      stop("Unknown palette type.")
    )
  }

  # Create raster layer in ggplot
  .gg_raster <- function(df, fill_col, aoi_layer = NULL, palette_type = "risk") {
    ggplot() +
      geom_tile(
        data = df,
        aes(x = Longitude, y = Latitude, fill = .data[[fill_col]])
      ) +
      .fill_scale(palette_type, fill_col) +
      aoi_layer
  }

  # Check whether overlap maps are rasters
  .is_raster_overlap <- function(overlap_one) {
    any(purrr::map_lgl(overlap_one, ~ {
      inherits(.x, "SpatRaster") ||
        (is.list(.x) && !is.null(.x$raster) && inherits(.x$raster, "SpatRaster"))
    }))
  }

  # Defining important objects
  risk_levels <- c("None", "Low", "Medium", "High")
  style <- NULL

  # risaMaps
  if (inherits(x, "risaMaps")) {
    sp_list <- x$species_kernel_maps
    st_list <- x$stressor_kernel_maps
    overlap_list <- x$overlap_maps
    species_names <- names(sp_list)
    stressor_names <- names(st_list)
    aoi <- x$area_of_interest

    # raster mode
    raster_mode <- .is_raster_overlap(overlap_list[[1]])

    if (raster_mode) {
      spp_df <- .rast_list_to_df(sp_list,  "raster", value_name = "value")
      stressor_df <- .rast_list_to_df(st_list,  "raster", value_name = "value")

      # overlap: map over species, each returns long df with its stressors
      overlap_df <- purrr::map_dfr(species_names, function(sp) {
        .rast_list_to_df(overlap_list[[sp]], "raster", value_name = "value") |>
          dplyr::mutate(species = sp)
      })


      aoi_layer <- NULL
      if (area_of_interest) {
        aoi_layer <- .ensure_aoi_layer(aoi, xy_df = spp_df, crs_like = sp_list[[1]]$raster)
      }

      kde_titles <- NULL
      if (!raw_ggplot) {
        kde_titles <- list(
          spp = ggtitle("Species KDE"),
          stress = ggtitle("Stressor KDE"),
          overlaps = ggtitle("Species-Stressor Overlaps")
        )
        style <- theme_bw()
      }

      gg_spp <- .gg_raster(spp_df, "value", aoi_layer, palette_type = "species") +
        facet_wrap(~ group, scales = "fixed") +
        style +
        kde_titles[["spp"]]

      gg_stress <- .gg_raster(stressor_df, "value", aoi_layer, palette_type = "stressor") +
        facet_wrap(~ group, scales = "fixed") +
        style +
        kde_titles[["stress"]]

      gg_overlap <- .gg_raster(overlap_df, "value", aoi_layer, palette_type = "overlap") +
        facet_grid(group ~ species, scales = "fixed") +
        style +
        kde_titles[["overlaps"]]

      return(list(gg_spp, gg_stress, gg_overlap))

    } else {
      # sf mode
      spp_sf <- if (length(sp_list) == 1) sp_list[[1]]$shp else .join_shps(sp_list)
      stressor_sf <- if (length(st_list) == 1) st_list[[1]]$shp else .join_shps(st_list)
      all_overlap <- purrr::map_dfr(species_names, ~ dplyr::mutate(.join_shps(overlap_list[[.x]]), species = .x))

      spp_sf <- spp_sf |>
        dplyr::mutate(Rating = .reclass_labels(Rating))

      stressor_sf <- stressor_sf |>
        dplyr::mutate(Rating = .reclass_labels(Rating))

      all_overlap <- all_overlap |>
        dplyr::mutate(Rating = .reclass_labels(Rating))

      aoi_layer <- NULL
      if (area_of_interest) {
        aoi_layer <- ggplot2::geom_sf(data = aoi, fill = "transparent", linewidth = 0.5)
      }

      if (!raw_ggplot) {
        kde_titles <- list(
          spp = ggtitle("Species KDE"),
          stress = ggtitle("Stressor KDE"),
          overlaps = ggtitle("Species-Stressor Overlaps")
        )
        style <- theme_bw()
      }

      gg_spp <- ggplot() +
        ggplot2::geom_sf(data = spp_sf, aes(fill = Rating), col = "transparent") +
        .fill_scale_discrete("species") +
        aoi_layer +
        facet_wrap(~ group, scales = "fixed") +
        style +
        kde_titles[["spp"]]

      gg_stress <- ggplot() +
        ggplot2::geom_sf(data = stressor_sf, aes(fill = Rating), col = "transparent") +
        .fill_scale_discrete("stressor") +
        aoi_layer +
        facet_wrap(~ group, scales = "fixed") +
        style +
        kde_titles[["stress"]]

      gg_overlap <- ggplot() +
        ggplot2::geom_sf(data = all_overlap, aes(fill = Rating), col = "transparent") +
        .fill_scale_discrete("overlap") +
        aoi_layer +
        facet_grid(group ~ species, scales = "fixed") +
        style +
        kde_titles[["overlaps"]]

      return(list(gg_spp, gg_stress, gg_overlap))
    }

  } else {
    # risaHRA
    depth <- list_depth_base(x)
    str_vecs <- unique(x$summary_stats$STRESSOR)
    stressor_names <- str_vecs[!str_vecs %in% "(FROM ALL STRESSORS)"]
    spp_vecs <- unique(x$summary_stats$SPECIES)
    species_names <- spp_vecs[!spp_vecs %in% "ECOSYSTEM"]

    # Prepare AOI layer when needed
    aoi_layer <- NULL

    if (depth == 2) {
      risk_raw <- terra::as.data.frame(x$total_raw, xy = TRUE) |>
        dplyr::rename(Longitude = x, Latitude = y) |>
        dplyr::rename(Risk = 3)
      risk_reclassified <- terra::as.data.frame(x$total_hotspots_reclassified, xy = TRUE) |>
        dplyr::rename(Longitude = x, Latitude = y) |>
        dplyr::rename(`Risk (reclass.)` = 3) |>
        dplyr::mutate(`Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`)) |>
        dplyr::mutate(`Risk (reclass.)` = factor(`Risk (reclass.)`, levels = risk_levels))

      if (is.null(aoi_layer)) {
        if (area_of_interest) {
          aoi_layer <- .ensure_aoi_layer(x$area_of_interest %||% NULL,
                                         xy_df = risk_raw, crs_like = x$total_raw)
        }
      }

      risk_titles <- NULL

      if (length(stressor_names) == 1) {

        if (!raw_ggplot) {
          risk_titles <- list(
            total = ggtitle("Total risk"),
            total_rec = ggtitle("Total risk (reclassified)"),
          )
          style <- theme_bw()
        }

        gg1 <- .gg_raster(risk_raw, "Risk", aoi_layer) +
          risk_titles[["total"]] + style

        gg2 <- .gg_raster(risk_reclassified, "Risk (reclass.)", aoi_layer) +
          risk_titles[["total_rec"]] + style

        return(list(gg1, gg2))
      } else {
        stressor_risk_raw <- .rast_list_to_df(x, "Risk_map_raw",
                                              group_names = stressor_names,
                                              value_name = "Risk")
        stressor_risk_reclass <- .rast_list_to_df(x, "Risk_map",
                                                  group_names = stressor_names,
                                                  value_name = "Risk (reclass.)") |>
          dplyr::mutate(`Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`))

        risk_max_ratings <- terra::as.data.frame(x$total_reclassified, xy = TRUE) |>
          dplyr::rename(Longitude = x, Latitude = y) |>
          dplyr::rename(`Highest str. risk (reclass.)` = 3) |>
          dplyr::mutate(`Highest str. risk (reclass.)` =
                          .reclass_labels(`Highest str. risk (reclass.)`)) |>
          dplyr::mutate(`Highest str. risk (reclass.)` =
                        factor(`Highest str. risk (reclass.)`, levels = risk_levels))

        if (!raw_ggplot) {
          risk_titles <- list(
            str_risk = ggtitle("Stressor risk"),
            str_risk_rec = ggtitle("Reclassified stressor risk"),
            h_str_risk = ggtitle("Highest stressor risk estimates"),
            comb_risk = ggtitle("Total (combined) risk"),
            comb_risk_rec = ggtitle("Total (combined) risk (reclassified)")
          )
          style <- theme_bw()
        }

        gg_str_risk_raw <- .gg_raster(stressor_risk_raw, "Risk", aoi_layer, palette_type = "risk") +
          facet_wrap(~ group, scales = "fixed") +
          risk_titles[["str_risk"]] + style

        gg_str_risk_reclass <- .gg_raster(stressor_risk_reclass, "Risk (reclass.)", aoi_layer, palette_type = "risk") +
          facet_wrap(~ group, scales = "fixed") +
          risk_titles[["str_risk_rec"]] + style

        gg_str_risk_max <- .gg_raster(risk_max_ratings, "Highest str. risk (reclass.)", aoi_layer, palette_type = "risk") +
          risk_titles[["h_str_risk"]] + style

        gg_risk_raw <- .gg_raster(risk_raw, "Risk", aoi_layer, palette_type = "risk") +
          risk_titles[["comb_risk"]] + style

        gg_risk_reclass <- .gg_raster(risk_reclassified, "Risk (reclass.)", aoi_layer, palette_type = "risk") +
          risk_titles[["comb_risk_rec"]] + style

        return(list(gg_str_risk_raw,
                    gg_str_risk_reclass,
                    gg_str_risk_max,
                    gg_risk_raw,
                    gg_risk_reclass))
      }

    } else {
      eco_risk_raw <- terra::as.data.frame(x$ecosys_risk_raw, xy = TRUE) |>
        dplyr::rename(Longitude = x, Latitude = y) |>
        dplyr::rename(Risk = 3)
      eco_risk_reclass <- terra::as.data.frame(x$ecosys_risk_classified, xy = TRUE) |>
        dplyr::rename(Longitude = x, Latitude = y) |>
        dplyr::rename(`Risk (reclass.)` = 3) |>
        dplyr::mutate(`Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`)) |>
        dplyr::mutate(`Risk (reclass.)` = factor(`Risk (reclass.)`, levels = risk_levels))

      aoi_layer <- NULL
      if (area_of_interest) {
        aoi_layer <- .ensure_aoi_layer(x$area_of_interest %||% NULL,
                                       xy_df = eco_risk_raw, crs_like = x$ecosys_risk_raw)
      }

      eco_titles <- NULL
      if (!raw_ggplot) {
        eco_titles <- list(
          eco_risk = ggtitle("Ecosystem risk"),
          eco_risk_rec = ggtitle("Reclassified ecosystem risk")
        )
        style <- theme_bw()
      }

      gg_ecosys_raw <- .gg_raster(eco_risk_raw, "Risk", aoi_layer, palette_type = "risk") +
        eco_titles[["eco_risk"]] + style
      gg_ecosys_reclass <- .gg_raster(eco_risk_reclass, "Risk (reclass.)", aoi_layer, palette_type = "risk") +
        eco_titles[["eco_risk_rec"]] + style

      # Per-species totals
      raw_risk_df <- purrr::map_dfr(species_names, ~ {
        terra::as.data.frame(x[[.x]]$total_raw, xy = TRUE) |>
          dplyr::rename(Longitude = x, Latitude = y) |>
          dplyr::rename(Risk = 3) |>
          dplyr::mutate(species = .x)
      })

      recl_risk_df <- purrr::map_dfr(species_names, ~ {
        terra::as.data.frame(x[[.x]]$total_reclassified, xy = TRUE) |>
          dplyr::rename(Longitude = x, Latitude = y) |>
          dplyr::rename(`Risk (reclass.)` = 3) |>
          dplyr::mutate(species = .x,
                        `Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`)) |>
          dplyr::mutate(`Risk (reclass.)` = factor(`Risk (reclass.)`, levels = risk_levels))
      })

      comb_titles <- NULL
      if (!raw_ggplot) {
        comb_titles <- list(
          comb_risk = ggtitle("Total (combined stress.) risk"),
          comb_risk_rec = ggtitle("Reclassified total (combined stress.) risk")
        )
        style <- theme_bw()
      }

      gg_raw_risk <- .gg_raster(raw_risk_df,  "Risk", aoi_layer, palette_type = "risk") +
        facet_grid(. ~ species) +
        comb_titles[["comb_risk"]] + style

      gg_recl_risk <- .gg_raster(recl_risk_df, "Risk (reclass.)", aoi_layer, palette_type = "risk") +
        facet_grid(. ~ species) +
        comb_titles[["comb_risk_rec"]] + style

      if (length(stressor_names) == 1) {
        return(list(gg_raw_risk, gg_recl_risk, gg_ecosys_raw, gg_ecosys_reclass))
      }

      # Species × stressor maps
      raw_str_risk_df <- purrr::map_dfr(species_names, function(sp) {
        purrr::map_dfr(stressor_names, function(st) {
          terra::as.data.frame(x[[sp]][[st]]$Risk_map_raw, xy = TRUE) |>
            dplyr::rename(Longitude = x, Latitude = y) |>
            dplyr::rename(Risk = 3) |>
            dplyr::mutate(species = sp, stressor = st)
        })
      })

      reclass_str_risk_df <- purrr::map_dfr(species_names, function(sp) {
        purrr::map_dfr(stressor_names, function(st) {
          terra::as.data.frame(x[[sp]][[st]]$Risk_map, xy = TRUE) |>
            dplyr::rename(Longitude = x, Latitude = y) |>
            dplyr::rename(`Risk (reclass.)` = 3) |>
            dplyr::mutate(species = sp, stressor = st,
                          `Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`)) |>
            dplyr::mutate(`Risk (reclass.)` = factor(`Risk (reclass.)`, levels = risk_levels))
        })
      })

      risk_titles <- NULL
      if (!raw_ggplot) {
        risk_titles <- list(
          str_risk = ggtitle("Stressor risk"),
          str_risk_rec = ggtitle("Reclassified stressor risk")
        )
        style <- theme_bw()
      }

      gg_str_raw_risk <- .gg_raster(raw_str_risk_df, "Risk", aoi_layer, palette_type = "risk") +
        facet_grid(stressor ~ species) +
        risk_titles[["str_risk"]] + style

      gg_str_recl_risk <- .gg_raster(reclass_str_risk_df, "Risk (reclass.)", aoi_layer, palette_type = "risk") +
        facet_grid(stressor ~ species) +
        risk_titles[["str_risk_rec"]] + style

      return(list(gg_str_raw_risk,
                  gg_str_recl_risk,
                  gg_raw_risk,
                  gg_recl_risk,
                  gg_ecosys_raw,
                  gg_ecosys_reclass))
    }
  }
}
