risaplot <- function(x, area_of_interest = TRUE, raw_ggplot = FALSE, target_crs = 4326) {
  
  stopifnot(inherits(x, c("risaMaps", "risaHRA")))
  
  # Constants
  risk_levels <- c("None", "Low", "Medium", "High")
  
  risk_cols <- c(
    "None"   = grDevices::hcl.colors(12, "YlOrRd")[10],
    "Low"    = grDevices::hcl.colors(12, "YlOrRd")[7],
    "Medium" = grDevices::hcl.colors(12, "YlOrRd")[4],
    "High"   = grDevices::hcl.colors(12, "YlOrRd")[1]
  )
  
  `%||%` <- function(a, b) {
    if (is.null(a)) b else a
  }
  
  base_theme <- if (raw_ggplot) NULL else ggplot2::theme_bw()
  
  
  # Helper functions
  # Helper function to transform rasters to target CRS
  .transform_raster_crs <- function(r, target_crs) {
    if (is.null(target_crs)) return(r)
    
    # Get current CRS
    current_crs <- terra::crs(r)
    
    # If already in target CRS, return as-is
    if (!is.na(current_crs) && 
        sf::st_crs(current_crs) == sf::st_crs(target_crs)) {
      return(r)
    }
    
    # Check if target CRS is geographic (decimal degrees)
    target_is_geo <- sf::st_crs(target_crs)$epsg %in% c(4326, 4978, 4979) ||
      grepl("longlat", sf::st_crs(target_crs)$wkt, ignore.case = TRUE)
    
    if (target_is_geo) {
      # For decimal degrees, we need to be careful with raster warping
      # Use project() for raster warping
      return(terra::project(r, target_crs, method = "bilinear"))
    } else {
      # For other CRS, use project()
      return(terra::project(r, target_crs))
    }
  }
  
  # Helper function to transform sf objects to target CRS
  .transform_sf_crs <- function(sf_obj, target_crs) {
    if (is.null(target_crs) || is.null(sf_obj)) return(sf_obj)
    
    current_crs <- sf::st_crs(sf_obj)
    
    # If already in target CRS or no CRS, return as-is
    if (is.na(current_crs) || 
        (!is.na(current_crs) && current_crs == sf::st_crs(target_crs))) {
      return(sf_obj)
    }
    
    return(sf::st_transform(sf_obj, target_crs))
  }
  
  .title <- function(label) {
    if (raw_ggplot) NULL else ggplot2::ggtitle(label)
  }
  
  .reclass_labels <- function(v) {
    labs <- c("0" = "None", "1" = "Low", "2" = "Medium", "3" = "High")
    
    factor(
      unname(labs[as.character(v)]),
      levels = risk_levels
    )
  }
  
  .get_raster <- function(obj, object_name = "raster", target_crs = NULL) {
    r <- NULL
    
    if (inherits(obj, "SpatRaster")) {
      r <- obj
    } else if (is.list(obj) &&
               !is.null(obj[[object_name]]) &&
               inherits(obj[[object_name]], "SpatRaster")) {
      r <- obj[[object_name]]
    } else {
      stop("Could not find a SpatRaster in this object.", call. = FALSE)
    }
    
    # Transform if target CRS is specified
    if (!is.null(target_crs)) {
      r <- .transform_raster_crs(r, target_crs)
    }
    
    return(r)
  }
  
  .rast_to_df <- function(r, value_name = "value", target_crs = NULL) {
    # If target CRS provided and r is a SpatRaster, transform it
    if (!is.null(target_crs) && inherits(r, "SpatRaster")) {
      r <- .transform_raster_crs(r, target_crs)
    }
    
    terra::as.data.frame(r, xy = TRUE, na.rm = FALSE) |>
      dplyr::rename(Longitude = x, Latitude = y) |>
      dplyr::rename(!!value_name := 3) |>
      dplyr::filter(!is.na(.data[[value_name]]))
  }
  
  .rast_list_to_df <- function(r_list,
                               object_name = "raster",
                               group_names = names(r_list),
                               value_name = "value",
                               group_col = "group",
                               target_crs = NULL) {
    
    r_stack <- terra::rast(
      purrr::map(group_names, \(nm) {
        .get_raster(r_list[[nm]], object_name, target_crs)
      })
    )
    
    names(r_stack) <- group_names
    
    terra::as.data.frame(r_stack, xy = TRUE, na.rm = FALSE) |>
      dplyr::rename(Longitude = x, Latitude = y) |>
      tidyr::pivot_longer(
        cols = -c(Longitude, Latitude),
        names_to = group_col,
        values_to = value_name
      ) |>
      dplyr::filter(!is.na(.data[[value_name]]))
  }
  
  .join_shps <- function(r_list, group_col = "group", target_crs = NULL) {
    purrr::imap(
      r_list,
      \(obj, nm) {
        if (is.null(obj$shp)) {
          stop("Could not find a `$shp` object in one of the list elements.", call. = FALSE)
        }
        
        shp <- obj$shp
        # Transform if target CRS is specified
        if (!is.null(target_crs)) {
          shp <- .transform_sf_crs(shp, target_crs)
        }
        
        dplyr::mutate(shp, !!group_col := nm)
      }
    ) |>
      dplyr::bind_rows()
  }
  
  .ensure_aoi_layer <- function(aoi, xy_df = NULL, crs_like = NULL) {
    if (inherits(aoi, c("sf", "sfc"))) {
      return(ggplot2::geom_sf(data = aoi, fill = "transparent", linewidth = 0.5))
    }
    
    if (is.null(xy_df) || is.null(crs_like)) {
      return(NULL)
    }
    
    bb <- sf::st_bbox(
      c(
        xmin = min(xy_df$Longitude, na.rm = TRUE),
        ymin = min(xy_df$Latitude, na.rm = TRUE),
        xmax = max(xy_df$Longitude, na.rm = TRUE),
        ymax = max(xy_df$Latitude, na.rm = TRUE)
      ),
      crs = sf::st_crs(crs_like)
    )
    
    ggplot2::geom_sf(
      data = sf::st_as_sf(sf::st_as_sfc(bb)),
      fill = "transparent",
      col = "transparent"
    )
  }
  
  .fill_scale <- function(type, fill_col = NULL) {
    is_reclass <- !is.null(fill_col) &&
      grepl("reclass", fill_col, ignore.case = TRUE)
    
    switch(
      type,
      
      species = ggplot2::scale_fill_gradientn(
        colors = grDevices::hcl.colors(120, "Greens", rev = TRUE)[31:120],
        na.value = "transparent"
      ),
      
      stressor = ggplot2::scale_fill_gradientn(
        colors = grDevices::hcl.colors(120, "Blues", rev = TRUE)[31:120],
        na.value = "transparent"
      ),
      
      overlap = ggplot2::scale_fill_gradientn(
        colors = grDevices::hcl.colors(100, "Viridis", rev = TRUE),
        na.value = "transparent"
      ),
      
      risk = {
        if (is_reclass) {
          ggplot2::scale_fill_manual(
            values = risk_cols,
            limits = risk_levels,
            breaks = risk_levels,
            drop = FALSE,
            na.value = "transparent",
            guide = ggplot2::guide_legend(
              override.aes = list(alpha = 1)
            )
          )
        } else {
          ggplot2::scale_fill_gradientn(
            colors = grDevices::hcl.colors(120, "YlOrRd", rev = TRUE)[21:120],
            na.value = "transparent"
          )
        }
      },
      
      stop("Unknown palette type.", call. = FALSE)
    )
  }
  
  .fill_scale_discrete <- function(type) {
    switch(
      type,
      
      species = ggplot2::scale_fill_brewer(
        palette = "Greens",
        drop = FALSE,
        na.value = "transparent"
      ),
      
      stressor = ggplot2::scale_fill_brewer(
        palette = "Blues",
        drop = FALSE,
        na.value = "transparent"
      ),
      
      overlap = ggplot2::scale_fill_viridis_d(
        option = "viridis",
        drop = FALSE,
        na.value = "transparent"
      ),
      
      risk = ggplot2::scale_fill_manual(
        values = risk_cols,
        limits = risk_levels,
        breaks = risk_levels,
        drop = FALSE,
        na.value = "transparent"
      ),
      
      stop("Unknown palette type.", call. = FALSE)
    )
  }
  
  .gg_raster <- function(df, fill_col, aoi_layer = NULL, palette_type = "risk") {
    palette_layer <- NULL
    if (!raw_ggplot) {
      palette_layer <- .fill_scale(palette_type, fill_col)
    }
    ggplot2::ggplot() +
      ggplot2::geom_tile(
        data = df,
        ggplot2::aes(x = Longitude, y = Latitude, fill = .data[[fill_col]]),
        show.legend = TRUE
      ) +
      palette_layer +
      aoi_layer
  }
  
  .gg_sf <- function(sf_obj, fill_col = "Rating", aoi_layer = NULL, palette_type = "risk") {
    palette_layer <- NULL
    if (!raw_ggplot) {
      palette_layer <- .fill_scale_discrete(palette_type)
    }
    ggplot2::ggplot() +
      ggplot2::geom_sf(
        data = sf_obj,
        ggplot2::aes(fill = .data[[fill_col]]),
        col = "transparent",
        show.legend = TRUE
      ) +
      palette_layer +
      aoi_layer
  }
  
  .is_raster_object <- function(obj, object_name = "raster") {
    inherits(obj, "SpatRaster") ||
      (is.list(obj) && !is.null(obj[[object_name]]) && inherits(obj[[object_name]], "SpatRaster"))
  }
  
  .is_raster_overlap <- function(overlap_one) {
    any(purrr::map_lgl(overlap_one, .is_raster_object))
  }
  
  .make_aoi <- function(aoi, xy_df, crs_like, target_crs = NULL) {
    if (is.null(aoi) || is.null(xy_df) || is.null(crs_like)) {
      return(NULL)
    }
    
    if (inherits(aoi, c("sf", "sfc"))) {
      # Transform AOI to target CRS if specified
      if (!is.null(target_crs)) {
        aoi <- .transform_sf_crs(aoi, target_crs)
      }
      return(ggplot2::geom_sf(data = aoi, fill = "transparent", linewidth = 0.5))
    }
    
    # Create bounding box using transformed coordinates
    bb <- sf::st_bbox(
      c(
        xmin = min(xy_df$Longitude, na.rm = TRUE),
        ymin = min(xy_df$Latitude, na.rm = TRUE),
        xmax = max(xy_df$Longitude, na.rm = TRUE),
        ymax = max(xy_df$Latitude, na.rm = TRUE)
      ),
      crs = sf::st_crs(target_crs %||% crs_like)
    )
    
    ggplot2::geom_sf(
      data = sf::st_as_sf(sf::st_as_sfc(bb)),
      fill = "transparent",
      col = "transparent"
    )
  }
  
  # risaMaps method
  if (inherits(x, "risaMaps")) {
    
    sp_list <- x$species_kernel_maps
    st_list <- x$stressor_kernel_maps
    overlap_list <- x$overlap_maps
    
    species_names <- names(sp_list)
    stressor_names <- names(st_list)
    
    raster_mode <- .is_raster_overlap(overlap_list[[1]])
    
    # SpatRaster mode
    if (raster_mode) {
      spp_df <- .rast_list_to_df(sp_list, object_name = "raster", 
                                 value_name = "value", target_crs = target_crs)
      stressor_df <- .rast_list_to_df(st_list, object_name = "raster", 
                                      value_name = "value", target_crs = target_crs)
      
      overlap_df <- purrr::map_dfr(
        species_names,
        \(sp) {
          .rast_list_to_df(
            overlap_list[[sp]],
            object_name = "raster",
            value_name = "value", 
            target_crs = target_crs
          ) |>
            dplyr::mutate(species = sp)
        }
      )
      
      aoi_layer <- .make_aoi(
        x$area_of_interest,
        xy_df = spp_df,
        crs_like = .get_raster(sp_list[[1]], "raster", target_crs = target_crs)
      )
      
      gg_spp <- .gg_raster(spp_df, "value", aoi_layer, "species") +
        ggplot2::facet_wrap(~ group, scales = "fixed") +
        base_theme +
        .title("Species KDE")
      
      gg_stress <- .gg_raster(stressor_df, "value", aoi_layer, "stressor") +
        ggplot2::facet_wrap(~ group, scales = "fixed") +
        base_theme +
        .title("Stressor KDE")
      
      gg_overlap <- .gg_raster(overlap_df, "value", aoi_layer, "overlap") +
        ggplot2::facet_grid(group ~ species, scales = "fixed") +
        base_theme +
        .title("Species-Stressor Overlaps")
      
      return(list(gg_spp, gg_stress, gg_overlap))
    }
    
    # sf mode
    spp_sf <- .join_shps(sp_list, target_crs = target_crs) |>
      dplyr::mutate(Rating = .reclass_labels(Rating))
    
    stressor_sf <- .join_shps(st_list, target_crs = target_crs) |>
      dplyr::mutate(Rating = .reclass_labels(Rating))
    
    all_overlap <- purrr::map_dfr(
      species_names,
      \(sp) {
        .join_shps(overlap_list[[sp]], target_crs = target_crs) |>
          dplyr::mutate(species = sp)
      }
    ) |>
      dplyr::mutate(Rating = .reclass_labels(Rating))
    
    aoi_layer <- if (area_of_interest) {
      .ensure_aoi_layer(x$area_of_interest)
    } else {
      NULL
    }
    
    gg_spp <- .gg_sf(spp_sf, "Rating", aoi_layer, "species") +
      ggplot2::facet_wrap(~ group, scales = "fixed") +
      base_theme +
      .title("Species KDE")
    
    gg_stress <- .gg_sf(stressor_sf, "Rating", aoi_layer, "stressor") +
      ggplot2::facet_wrap(~ group, scales = "fixed") +
      base_theme +
      .title("Stressor KDE")
    
    gg_overlap <- .gg_sf(all_overlap, "Rating", aoi_layer, "overlap") +
      ggplot2::facet_grid(group ~ species, scales = "fixed") +
      base_theme +
      .title("Species-Stressor Overlaps")
    
    return(list(gg_spp, gg_stress, gg_overlap))
  }
  
  # risaHRA method
  depth <- list_depth_base(x)
  
  stressor_names <- unique(x$summary_stats$STRESSOR)
  stressor_names <- stressor_names[stressor_names != "(FROM ALL STRESSORS)"]
  
  species_names <- unique(x$summary_stats$SPECIES)
  species_names <- species_names[species_names != "ECOSYSTEM"]
  
  # risaHRA depth == 2
  if (depth == 2) {
    risk_raw <- .rast_to_df(x$total_raw, value_name = "Risk", target_crs = target_crs)
    
    risk_reclassified <- .rast_to_df(
      x$total_hotspots_reclassified,
      value_name = "Risk (reclass.)",
      target_crs = target_crs
    ) |>
      dplyr::mutate(`Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`))
    
    aoi_layer <- .make_aoi(
      x$area_of_interest %||% NULL,
      xy_df = risk_raw,
      crs_like = x$total_raw, 
      target_crs = target_crs
    )
    
    if (length(stressor_names) == 1) {
      gg1 <- .gg_raster(risk_raw, "Risk", aoi_layer, "risk") +
        base_theme +
        .title("Total risk")
      
      gg2 <- .gg_raster(risk_reclassified, "Risk (reclass.)", aoi_layer, "risk") +
        base_theme +
        .title("Total risk (reclassified)")
      
      return(list(gg1, gg2))
    }
    
    stressor_risk_raw <- .rast_list_to_df(
      x,
      object_name = "Risk_map_raw",
      group_names = stressor_names,
      value_name = "Risk", 
      target_crs = target_crs
    )
    
    stressor_risk_reclass <- .rast_list_to_df(
      x,
      object_name = "Risk_map",
      group_names = stressor_names,
      value_name = "Risk (reclass.)", 
      target_crs = target_crs
    ) |>
      dplyr::mutate(`Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`))
    
    risk_max_ratings <- .rast_to_df(
      x$total_reclassified,
      value_name = "Highest str. risk (reclass.)", 
      target_crs = target_crs
    ) |>
      dplyr::mutate(
        `Highest str. risk (reclass.)` = .reclass_labels(`Highest str. risk (reclass.)`)
      )
    
    gg_str_risk_raw <- .gg_raster(stressor_risk_raw, "Risk", aoi_layer, "risk") +
      ggplot2::facet_wrap(~ group, scales = "fixed") +
      base_theme +
      .title("Stressor risk")
    
    gg_str_risk_reclass <- .gg_raster(stressor_risk_reclass, "Risk (reclass.)", aoi_layer, "risk") +
      ggplot2::facet_wrap(~ group, scales = "fixed") +
      base_theme +
      .title("Reclassified stressor risk")
    
    gg_str_risk_max <- .gg_raster(risk_max_ratings, "Highest str. risk (reclass.)", aoi_layer, "risk") +
      base_theme +
      .title("Highest stressor risk estimates")
    
    gg_risk_raw <- .gg_raster(risk_raw, "Risk", aoi_layer, "risk") +
      base_theme +
      .title("Total (combined) risk")
    
    gg_risk_reclass <- .gg_raster(risk_reclassified, "Risk (reclass.)", aoi_layer, "risk") +
      base_theme +
      .title("Total (combined) risk (reclassified)")
    
    return(list(
      gg_str_risk_raw,
      gg_str_risk_reclass,
      gg_str_risk_max,
      gg_risk_raw,
      gg_risk_reclass
    ))
  }
  
  # risaHRA depth > 2
  eco_risk_raw <- .rast_to_df(x$ecosys_risk_raw, value_name = "Risk", target_crs = target_crs)
  
  eco_risk_reclass <- .rast_to_df(
    x$ecosys_risk_classified,
    value_name = "Risk (reclass.)", 
    target_crs = target_crs
  ) |>
    dplyr::mutate(`Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`))
  
  aoi_layer <- .make_aoi(
    x$area_of_interest %||% NULL,
    xy_df = eco_risk_raw,
    crs_like = x$ecosys_risk_raw, 
    target_crs = target_crs
  )
  
  gg_ecosys_raw <- .gg_raster(eco_risk_raw, "Risk", aoi_layer, "risk") +
    base_theme +
    .title("Ecosystem risk")
  
  gg_ecosys_reclass <- .gg_raster(eco_risk_reclass, "Risk (reclass.)", aoi_layer, "risk") +
    base_theme +
    .title("Reclassified ecosystem risk")
  
  raw_risk_df <- purrr::map_dfr(
    species_names,
    \(sp) {
      .rast_to_df(x[[sp]]$total_raw, value_name = "Risk", target_crs = target_crs) |>
        dplyr::mutate(species = sp)
    }
  )
  
  recl_risk_df <- purrr::map_dfr(
    species_names,
    \(sp) {
      .rast_to_df(x[[sp]]$total_reclassified, 
                  value_name = "Risk (reclass.)", 
                  target_crs = target_crs) |>
        dplyr::mutate(
          species = sp,
          `Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`)
        )
    }
  )
  
  gg_raw_risk <- .gg_raster(raw_risk_df, "Risk", aoi_layer, "risk") +
    ggplot2::facet_grid(. ~ species) +
    base_theme +
    .title("Total (combined stress.) risk")
  
  gg_recl_risk <- .gg_raster(recl_risk_df, "Risk (reclass.)", aoi_layer, "risk") +
    ggplot2::facet_grid(. ~ species) +
    base_theme +
    .title("Reclassified total (combined stress.) risk")
  
  if (length(stressor_names) == 1) {
    return(list(
      gg_raw_risk,
      gg_recl_risk,
      gg_ecosys_raw,
      gg_ecosys_reclass
    ))
  }
  
  raw_str_risk_df <- purrr::map_dfr(
    species_names,
    \(sp) {
      .rast_list_to_df(
        x[[sp]],
        object_name = "Risk_map_raw",
        group_names = stressor_names,
        value_name = "Risk",
        group_col = "stressor", 
        target_crs = target_crs
      ) |>
        dplyr::mutate(species = sp)
    }
  )
  
  reclass_str_risk_df <- purrr::map_dfr(
    species_names,
    \(sp) {
      .rast_list_to_df(
        x[[sp]],
        object_name = "Risk_map",
        group_names = stressor_names,
        value_name = "Risk (reclass.)",
        group_col = "stressor", 
        target_crs = target_crs
      ) |>
        dplyr::mutate(
          species = sp,
          `Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`)
        )
    }
  )
  
  gg_str_raw_risk <- .gg_raster(raw_str_risk_df, "Risk", aoi_layer, "risk") +
    ggplot2::facet_grid(stressor ~ species) +
    base_theme +
    .title("Stressor risk")
  
  gg_str_recl_risk <- .gg_raster(reclass_str_risk_df, "Risk (reclass.)", aoi_layer, "risk") +
    ggplot2::facet_grid(stressor ~ species) +
    base_theme +
    .title("Reclassified stressor risk")
  
  return(list(
    gg_str_raw_risk,
    gg_str_recl_risk,
    gg_raw_risk,
    gg_recl_risk,
    gg_ecosys_raw,
    gg_ecosys_reclass
  ))
}