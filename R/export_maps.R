#' Export nested list of maps to files
#'
#' Recursively writes `terra::SpatRaster` (.tif) and vector layers using the
#' chosen driver (ESRI Shapefile, GPKG, or GeoJSON) from a nested list to disk.
#' Optionally zips the top-level output directory.
#'
#' @param x A nested list containing `SpatRaster`, `sf`/`sfc`, or further lists.
#' @param out_dir Character path to output directory.
#' @param path Internal use: character vector of names under recursion.
#' @param zip_export Logical; if `TRUE`, zip the top-level output directory.
#' @param vector_driver One of `"ESRI Shapefile"`, `"GPKG"`, or `"GeoJSON"`.
#'   Controls the format for vector outputs. Default `"ESRI Shapefile"`.
#' @return Invisibly returns `NULL`. Files are written to `out_dir`.
#' @importFrom terra writeRaster
#' @importFrom sf st_write st_as_sf
#' @importFrom utils zip
#' @examples
#' # Creating test data
#' spp_df <- rbind(data.frame(long = rnorm(80, 0, 10),
#'                            lat = rnorm(80, 0, 10), species = "sp1"),
#'                 data.frame(long = rnorm(60, 0, 10),
#'                            lat = rnorm(60, 0, 10), species = "sp2"))
#' str_df <- rbind(data.frame(long = rnorm(100, 0, 5),
#'                            lat = rnorm(100, 0, 10), stressor = "trawling"),
#'                 data.frame(long = rnorm(50, 0, 10),
#'                            lat = rnorm(100, 0, 5), stressor = "gillnet"))
#'
#' # Create kernel maps of species and stressor distributions and overlap maps
#' risa_maps <- risa_prep(spp_df, str_df)
#'
#' # Export vector and raster objects as shapefiles and tiff to external folder
#' export_maps(risa_maps, "my_folder")
#' @export
export_maps <- function(
    x,
    out_dir,
    path = character(),
    zip_export = FALSE,
    vector_driver = c("ESRI Shapefile","GPKG","GeoJSON")
) {
  vector_driver <- match.arg(vector_driver)

  # Ensure destination exists
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Name helpers
  sanitize_name <- function(s) {
    s <- if (is.null(s) || !nzchar(s)) "layer" else as.character(s)
    s <- gsub("[^A-Za-z0-9_\\-]+", "_", s)
    s <- sub("^_+", "", s)
    s <- substr(s, 1L, 60L)
    if (!nzchar(s)) "layer" else s
  }
  base_name <- if (length(path)) sanitize_name(tail(path, 1)) else "layer"

  # Writers
  write_raster_safely <- function(r, fname_base) {
    out_file <- file.path(out_dir, paste0(fname_base, ".tif"))
    tryCatch(
      terra::writeRaster(
        r,
        filename  = out_file,
        filetype  = "GTiff",
        overwrite = TRUE
      ),
      error = function(e) warning("Failed to write raster: ", out_file, " (", e$message, ")", call. = FALSE)
    )
  }

  # Pick dsn + layer based on vector_driver
  vector_target <- function(folder, fname_base, driver) {
    if (driver == "ESRI Shapefile") {
      dsn <- folder
      layer <- fname_base
    } else if (driver == "GPKG") {
      dsn <- file.path(folder, paste0(fname_base, ".gpkg"))
      layer <- fname_base
    } else { # GeoJSON
      dsn <- file.path(folder, paste0(fname_base, ".geojson"))
      layer <- fname_base
    }
    list(dsn = dsn, layer = layer)
  }

  write_vector_safely <- function(v, layer_base) {
    if (inherits(v, "sfc")) v <- sf::st_as_sf(v)
    tgt <- vector_target(out_dir, layer_base, vector_driver)
    drv <- vector_driver
    tryCatch(
      sf::st_write(
        obj = v,
        dsn = tgt$dsn,
        layer = tgt$layer,
        driver = drv,
        delete_layer = TRUE,
        quiet = TRUE
      ),
      error = function(e) {
        dest <- if (drv == "ESRI Shapefile") file.path(out_dir, paste0(layer_base, ".shp")) else tgt$dsn
        warning("Failed to write vector: ", dest, " (", e$message, ")", call. = FALSE)
      }
    )
  }

  # Dispatch
  if (inherits(x, "SpatRaster")) {
    write_raster_safely(x, base_name)

  } else if (inherits(x, "sf") || inherits(x, "sfc")) {
    write_vector_safely(x, base_name)

  } else if (is.list(x)) {
    nms <- names(x)
    for (i in seq_along(x)) {
      nm <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else sprintf("item_%02d", i)
      nm <- sanitize_name(nm)
      sub_dir <- file.path(out_dir, nm)
      export_maps(
        x = x[[i]],
        out_dir = sub_dir,
        path = c(path, nm),
        zip_export = FALSE,
        vector_driver = vector_driver
      )
    }

  } else {
    message("Skipping unrecognized object at: ", paste(path, collapse = "$"))
  }

  # Zip only once at the top level
  if (isTRUE(zip_export) && length(path) == 0L) {
    zip_file <- paste0(normalizePath(out_dir, mustWork = FALSE), ".zip")
    message("Zipping export folder to: ", zip_file)
    utils::zip(
      zipfile = zip_file,
      files   = list.files(out_dir, full.names = TRUE, recursive = TRUE)
    )
  }

  invisible(NULL)
}
