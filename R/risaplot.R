#' Plot RISA spatial outputs
#'
#' Generates publication-ready `ggplot2` maps from RISA output objects
#' (either `risaMaps` or `risaHRA` classes). It automatically detects whether
#' the inputs are rasters (`SpatRaster`) or vector (`sf`) objects and
#' produces a consistent set of plots for species, stressors, overlaps, and
#' risk outputs.
#'
#' @param x An object of class `risaMaps` or `risaHRA`.  
#'   - For `risaMaps`: the object must contain lists of species and stressor
#'     kernel density maps (`species_kernel_maps`, `stressor_kernel_maps`),
#'     and optionally `overlap_maps` and an `area_of_interest` polygon.  
#'   - For `risaHRA`: the object must contain total risk and reclassified
#'     maps (e.g., `total_raw`, `total_hotspots_reclassified`), and optionally
#'     per-stressor risk layers and `area_of_interest`.
#'
#' @return
#' A **list of `ggplot` objects**, typically including maps for species,
#' stressors, overlaps, and total risk (depending on input class).
#'
#' @section For `risaMaps` objects:
#' When `x` is of class `risaMaps`, the function returns a list of three plots:
#' \itemize{
#'   \item \strong{Species KDE:} Kernel density maps of each species'
#'     occurrence probability or spatial use.
#'   \item \strong{Stressor KDE:} Kernel density maps of each stressor’s
#'     spatial influence or exposure intensity.
#'   \item \strong{Species–Stressor Overlaps:} Combined maps showing
#'     the co-occurrence of species and stressor distributions.
#' }
#' All layers are faceted by species or stressor name, and overlaid with the
#' area of interest (AOI) when available.
#'
#' @section For `risaHRA` objects:
#' When `x` is of class `risaHRA`, the function returns a list of risk maps
#' depending on the hierarchical depth of the object (i.e., if one or more species
#' / habitats are considered.):
#' \describe{
#'   \item{Depth 2 (species-level HRA):}{
#'     \itemize{
#'       \item \strong{Stressor risk (raw):} Per-stressor total risk intensity.
#'       \item \strong{Stressor risk (reclassified):} Per-stressor risk ratings
#'         converted to ordinal classes ("None", "Low", "Medium", "High").
#'       \item \strong{Highest stressor risk:} Highest risk level across
#'         stressors (reclassified).
#'       \item \strong{Total (combined) risk:} Aggregated total risk across
#'         all stressors.
#'       \item \strong{Total (combined) risk, reclassified:} Reclassified
#'         ordinal version of the total risk.
#'     }
#'   }
#'   \item{Depth > 2 (ecosystem-level HRA):}{
#'     \itemize{
#'       \item \strong{Species-level risk maps:} Raw and reclassified risk
#'         maps for each modeled species.
#'       \item \strong{Stressor-specific risk maps:} Raw and reclassified
#'         species–stressor risk maps (faceted by species and stressor).
#'       \item \strong{Ecosystem risk (raw and reclassified):} Overall
#'         combined risk maps integrating all species and stressors.
#'     }
#'   }
#' }
#'
#' @details
#' The function handles both raster- and vector-based outputs:
#' \describe{
#'   \item{Raster inputs}{All `SpatRaster` objects are converted to data frames
#'     and plotted using `geom_tile()` for efficient rendering. The resulting
#'     plots are faceted by species, stressor, or overlap group.}
#'   \item{Vector inputs}{All `sf` objects are combined with a grouping column
#'     and plotted using `geom_sf()`. Faceting follows the same structure as for
#'     rasters.}
#' }
#'
#' Reclassified values are automatically converted to ordered factors with
#' levels: `"None"`, `"Low"`, `"Medium"`, and `"High"`.
#'
#' The area of interest (AOI) is overlaid as a transparent polygon when
#' available.
#'
#' @section Performance:
#' Internally, this function uses `purrr::map_dfr()` and `vctrs::vec_rbind()`
#' for fast concatenation of multiple raster or vector layers, avoiding costly
#' iterative `rbind()` calls. `geom_raster()` is used instead of `geom_tile()`
#' to optimize large raster plotting.
#'
#' @examples
#' \dontrun{
#' # Example with risaMaps
#' plots <- risaplot(risa_obj)
#' plots[[1]] # Plot first ggplot output
#'
#' # Example with risaHRA
#' hra_plots <- risaplot(hra_obj)
#' patchwork::wrap_plots(hra_plots)
#' }
#'
#' @seealso
#' [risa::risaMaps()], [risa::risaHRA()], [ggplot2::geom_raster()],
#' [ggplot2::geom_sf()], [vctrs::vec_rbind()]
#'
#' @export
risaplot <- function(x) {
  
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
      return(geom_sf(data = aoi, fill = "transparent", linewidth = 0.5))
    }

    stopifnot(!is.null(xy_df))
    stopifnot(!is.null(crs_like))
    bb <- st_bbox(
      c(xmin = min(xy_df$Longitude, na.rm = TRUE),
        ymin = min(xy_df$Latitude,  na.rm = TRUE),
        xmax = max(xy_df$Longitude, na.rm = TRUE),
        ymax = max(xy_df$Latitude,  na.rm = TRUE)),
      crs = st_crs(crs_like)
    )
    geom_sf(data = st_as_sf(st_as_sfc(bb)),
            fill = "transparent", col = "transparent")
  }
  
  # Convert a raster list into a data.frame
  .rast_list_to_df <- function(r_list, object_name = "raster",
                               group_names = NULL, value_name = "value") {
    if (is.null(group_names)) group_names <- names(r_list)
    
    r_stack <- terra::rast(map(group_names, ~ r_list[[.x]][[object_name]]))
    names(r_stack) <- group_names
    df <- terra::as.data.frame(r_stack, xy = TRUE)
    names(df)[1:2] <- c("Longitude","Latitude")
    df |>
      pivot_longer(-c(Longitude,Latitude),
                   names_to = "group",
                   values_to = value_name) |>
      filter(!is.na(.data[[value_name]]))
  }
  
  #Merge multiple sf objects
  .join_shps <- function(r_list, group_col = "group") {
    # bind_rows is fast and preserves sf; add group from names
    imap(r_list, ~ mutate(.x$shp, !!group_col := .y)) |>
      vctrs::vec_rbind()
  }
  
  # Create raster layer in ggplot
  .gg_raster <- function(df, fill_col, aoi_layer = NULL) {
    ggplot() +
      geom_tile(data = df, aes(x = Longitude, y = Latitude, fill = .data[[fill_col]])) +
      aoi_layer +
      coord_sf(expand = FALSE)
  }
  
  # checks if rasters
  .is_raster_overlap <- function(overlap_one) {
    any(map_lgl(overlap_one, ~ inherits(.x$raster, "SpatRaster")))
  }
  
  # risaMaps
  if (inherits(x, "risaMaps")) {
    sp_list       <- x$species_kernel_maps
    st_list       <- x$stressor_kernel_maps
    overlap_list  <- x$overlap_maps
    species_names <- names(sp_list)
    stressor_names<- names(st_list)
    aoi           <- x$area_of_interest
    
    # raster mode
    raster_mode <- .is_raster_overlap(overlap_list[[1]])
    
    if (raster_mode) {
      spp_df      <- .rast_list_to_df(sp_list,  "raster", value_name = "value")
      stressor_df <- .rast_list_to_df(st_list,  "raster", value_name = "value")
      
      # overlap: map over species, each returns long df with its stressors
      overlap_df <- map_dfr(species_names, function(sp) {
        .rast_list_to_df(overlap_list[[sp]], "raster", value_name = "value") |>
          mutate(species = sp)
      })
      
      aoi_layer <- .ensure_aoi_layer(aoi, xy_df = spp_df, crs_like = sp_list[[1]]$raster)
      
      gg_spp <- .gg_raster(spp_df, "value", aoi_layer) +
        facet_wrap(~ group, scales = "fixed") +
        ggtitle("Species KDE")
      
      gg_stress <- .gg_raster(stressor_df, "value", aoi_layer) +
        facet_wrap(~ group, scales = "fixed") +
        ggtitle("Stressor KDE")
      
      gg_overlap <- .gg_raster(overlap_df, "value", aoi_layer) +
        facet_grid(group ~ species, scales = "fixed") +
        ggtitle("Species-Stressor Overlaps")
      
      return(list(gg_spp, gg_stress, gg_overlap))
    } else {
      # shp mode
      spp_sf <- if (length(sp_list) == 1) sp_list[[1]]$shp else .join_shps(sp_list)
      stressor_sf <- if (length(st_list) == 1) st_list[[1]]$shp else .join_shps(st_list)
      all_overlap <- map_dfr(species_names, ~ mutate(.join_shps(overlap_list[[.x]]), species = .x))
      
      aoi_layer <- geom_sf(data = aoi, fill = "transparent", linewidth = 0.5)
      
      gg_spp <- ggplot() +
        geom_sf(data = spp_sf, aes(fill = Rating), col = "transparent") +
        aoi_layer + facet_wrap(~ group, scales = "fixed") +
        ggtitle("Species KDE")
      
      gg_stress <- ggplot() +
        geom_sf(data = stressor_sf, aes(fill = Rating), col = "transparent") +
        aoi_layer + facet_wrap(~ group, scales = "fixed") +
        ggtitle("Stressor KDE")
      
      gg_overlap <- ggplot() +
        geom_sf(data = all_overlap, aes(fill = Rating), col = "transparent") +
        aoi_layer + facet_grid(group ~ species, scales = "fixed") +
        ggtitle("Species-Stressor Overlaps")
      
      return(list(gg_spp, gg_stress, gg_overlap))
    }
    
  } else {
    # risaHRA
    depth          <- list_depth_base(x)
    str_vecs       <- unique(x$summary_stats$STRESSOR)
    stressor_names <- str_vecs[!str_vecs %in% "(FROM ALL STRESSORS)"]
    spp_vecs       <- unique(x$summary_stats$SPECIES)
    species_names  <- spp_vecs[!spp_vecs %in% "ECOSYSTEM"]
    
    # Prepare AOI layer lazily when needed
    aoi_layer <- NULL
    
    if (depth == 2) {
      risk_raw <- terra::as.data.frame(x$total_raw, xy = TRUE) |>
        rename(Longitude = x, Latitude = y) |>
        rename(Risk = 3)
      risk_reclassified <- terra::as.data.frame(x$total_hotspots_reclassified, xy = TRUE) |>
        rename(Longitude = x, Latitude = y) |>
        rename(`Risk (reclass.)` = 3) |>
        mutate(`Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`))
      
      if (is.null(aoi_layer)) {
        aoi_layer <- .ensure_aoi_layer(x$area_of_interest %||% NULL,
                                       xy_df = risk_raw, crs_like = x$total_raw)
      }
      
      if (length(stressor_names) == 1) {
        gg1 <- .gg_raster(risk_raw, "Risk", aoi_layer) + ggtitle("Total risk")
        gg2 <- .gg_raster(risk_reclassified, "Risk (reclass.)", aoi_layer) +
          ggtitle("Total risk (reclassified)")
        return(list(gg1, gg2))
      } else {
        stressor_risk_raw <- .rast_list_to_df(x, "Risk_map_raw",
                                              group_names = stressor_names,
                                              value_name = "Risk")
        stressor_risk_reclass <- .rast_list_to_df(x, "Risk_map",
                                                  group_names = stressor_names,
                                                  value_name = "Risk (reclass.)") |>
          mutate(`Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`))
        
        risk_max_ratings <- terra::as.data.frame(x$total_reclassified, xy = TRUE) |>
          rename(Longitude = x, Latitude = y) |>
          rename(`Highest str. risk (reclass.)` = 3) |>
          mutate(`Highest str. risk (reclass.)` =
                   .reclass_labels(`Highest str. risk (reclass.)`))
        
        gg_str_risk_raw <- .gg_raster(stressor_risk_raw, "Risk", aoi_layer) +
          facet_wrap(~ group, scales = "fixed") + ggtitle("Stressor risk")
        
        gg_str_risk_reclass <- .gg_raster(stressor_risk_reclass, "Risk (reclass.)", aoi_layer) +
          facet_wrap(~ group, scales = "fixed") + ggtitle("Reclassified stressor risk")
        
        gg_str_risk_max <- .gg_raster(risk_max_ratings, "Highest str. risk (reclass.)", aoi_layer) +
          ggtitle("Highest stressor risk estimates")
        
        gg_risk_raw <- .gg_raster(risk_raw, "Risk", aoi_layer) +
          ggtitle("Total (combined) risk")
        
        gg_risk_reclass <- .gg_raster(risk_reclassified, "Risk (reclass.)", aoi_layer) +
          ggtitle("Total (combined) risk (reclassified)")
        
        return(list(gg_str_risk_raw,
                    gg_str_risk_reclass,
                    gg_str_risk_max,
                    gg_risk_raw,
                    gg_risk_reclass))
      }
      
    } else {
      eco_risk_raw <- terra::as.data.frame(x$ecosys_risk_raw, xy = TRUE) |>
        rename(Longitude = x, Latitude = y) |>
        rename(Risk = 3)
      eco_risk_reclass <- terra::as.data.frame(x$ecosys_risk_classified, xy = TRUE) |>
        rename(Longitude = x, Latitude = y) |>
        rename(`Risk (reclass.)` = 3) |>
        mutate(`Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`))
      
      aoi_layer <- .ensure_aoi_layer(x$area_of_interest %||% NULL,
                                     xy_df = eco_risk_raw, crs_like = x$ecosys_risk_raw)
      
      gg_ecosys_raw     <- .gg_raster(eco_risk_raw, "Risk", aoi_layer) + ggtitle("Ecosystem risk")
      gg_ecosys_reclass <- .gg_raster(eco_risk_reclass, "Risk (reclass.)", aoi_layer) +
        ggtitle("Reclassified ecosystem risk")
      
      # Per-species totals
      raw_risk_df <- map_dfr(species_names, ~ {
        terra::as.data.frame(x[[.x]]$total_raw, xy = TRUE) |>
          rename(Longitude = x, Latitude = y) |>
          rename(Risk = 3) |>
          mutate(species = .x)
      })
      
      recl_risk_df <- map_dfr(species_names, ~ {
        terra::as.data.frame(x[[.x]]$total_reclassified, xy = TRUE) |>
          rename(Longitude = x, Latitude = y) |>
          rename(`Risk (reclass.)` = 3) |>
          mutate(species = .x,
                 `Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`))
      })
      
      gg_raw_risk  <- .gg_raster(raw_risk_df,  "Risk", aoi_layer) +
        facet_grid(. ~ species) + ggtitle("Total (combined stress.) risk")
      gg_recl_risk <- .gg_raster(recl_risk_df, "Risk (reclass.)", aoi_layer) +
        facet_grid(. ~ species) + ggtitle("Reclassified total (combined stress.) risk")
      
      if (length(stressor_names) == 1) {
        return(list(gg_raw_risk, gg_recl_risk, gg_ecosys_raw, gg_ecosys_reclass))
      }
      
      # Species × stressor maps
      raw_str_risk_df <- map_dfr(species_names, function(sp) {
        map_dfr(stressor_names, function(st) {
          terra::as.data.frame(x[[sp]][[st]]$Risk_map_raw, xy = TRUE) |>
            rename(Longitude = x, Latitude = y) |>
            rename(Risk = 3) |>
            mutate(species = sp, stressor = st)
        })
      })
      
      reclass_str_risk_df <- map_dfr(species_names, function(sp) {
        map_dfr(stressor_names, function(st) {
          terra::as.data.frame(x[[sp]][[st]]$Risk_map, xy = TRUE) |>
            rename(Longitude = x, Latitude = y) |>
            rename(`Risk (reclass.)` = 3) |>
            mutate(species = sp, stressor = st,
                   `Risk (reclass.)` = .reclass_labels(`Risk (reclass.)`))
        })
      })
      
      gg_str_raw_risk <- .gg_raster(raw_str_risk_df, "Risk", aoi_layer) +
        facet_grid(stressor ~ species) + ggtitle("Stressor risk")
      
      gg_str_recl_risk <- .gg_raster(reclass_str_risk_df, "Risk (reclass.)", aoi_layer) +
        facet_grid(stressor ~ species) + ggtitle("Reclassified stressor risk")
      
      return(list(gg_str_raw_risk,
                  gg_str_recl_risk,
                  gg_raw_risk,
                  gg_recl_risk,
                  gg_ecosys_raw,
                  gg_ecosys_reclass))
    }
  }
}
