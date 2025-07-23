#' Export nested list of maps to files
#'
#' Recursively writes `SpatRaster` and `sf` layers, respectively as `.tif` and shapefiles, from a nested list to disk, optionally zipping.
#'
#' @param x A nested list of `SpatRaster`, `sf`, or further lists.
#' @param out_dir Character path to output directory.
#' @param path Internal use: character vector of path under recursion.
#' @param zip_export Logical; if `TRUE`, zip the top-level output.
#' @return Invisibly returns `NULL`. Files are written to `out_dir`.
#' @importFrom terra writeRaster
#' @importFrom sf st_write
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
export_maps <- function(x, out_dir, path = character(), zip_export = FALSE) {
  # Ensure out_dir exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Helper to pick the filename from the current path
  get_base_name <- function(path) {
    if (length(path) > 0) {
      tail(path, 1)
    } else {
      "layer"
    }
  }
  base_name <- get_base_name(path)

  # Handle SpatRaster
  if (inherits(x, "SpatRaster")) {
    out_file <- file.path(out_dir, paste0(base_name, ".tif"))
    terra::writeRaster(x, filename = out_file, overwrite = TRUE)

    # Handle sf
  } else if (inherits(x, "sf")) {
    sf::st_write(
      obj = x,
      dsn = out_dir,
      layer = base_name,
      driver = "ESRI Shapefile",
      delete_layer = TRUE
    )

    # Recurse into nested lists
  } else if (is.list(x)) {
    for (nm in names(x)) {
      sub_dir <- file.path(out_dir, nm)
      export_maps(
        x = x[[nm]],
        out_dir = sub_dir,
        path = c(path, nm),
        zip_export = FALSE
      )
    }

    # Skip anything else
  } else {
    message("Skipping unrecognized object: ", paste(path, collapse = "$"))
  }

  # Zip the top-level folder
  if (zip_export && identical(sys.parent(), 0L)) {
    zip_file <- paste0(out_dir, ".zip")
    message("Zipping export folder to: ", zip_file)
    # I don't know if this will work in all systems, but let's prey for the best
    zip(zipfile = zip_file,
        files = list.files(out_dir, full.names = TRUE, recursive = TRUE))
  }
}
