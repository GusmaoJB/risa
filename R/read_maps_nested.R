#' Read spatial files into a nested list mirroring directory structure
#'
#' Scans a root directory (non-recursively at each level) for raster files
#' (`.tif`, `.tiff`) and vector files (`.shp`; also supports `.gpkg`, `.geojson`)
#' and reads them as `terra::SpatRaster` or `sf` objects. Subdirectories are
#' traversed recursively; the returned list mirrors the folder structure.
#'
#' @param dir_path character(1). Path to the root directory containing spatial data.
#' @return A nested `list`:
#'   - Intermediate elements are named by subdirectory names and contain further lists.
#'   - Leaf elements are named by file basenames (without extensions) and contain
#'     `SpatRaster` or `sf` objects.
#' @importFrom terra rast
#' @importFrom sf st_read
#' @importFrom tools file_ext file_path_sans_ext
#' @examples
#' \dontrun{
#' # Suppose your files are organized like:
#' # data/
#' # ├─ regionA/
#' # │  ├─ landcover.tif
#' # │  └─ roads.shp
#' # └─ regionB/
#' #    └─ zone1/
#' #       └─ elevation.tif
#'
#' maps <- read_maps_nested("data")
#' # Access the landcover raster:
#' maps$regionA$landcover
#' # Access the roads sf object:
#' maps$regionA$roads
#' # Access the elevation raster:
#' maps$regionB$zone1$elevation
#' }
#' @export
read_maps_nested <- function(dir_path) {
  if (!is.character(dir_path) || length(dir_path) != 1L) {
    stop("`dir_path` must be a single character path.")
  }
  if (!dir.exists(dir_path)) {
    stop("Directory not found: ", dir_path)
  }

  # internal recursive worker
  .read_dir <- function(d) {
    out <- list()

    # files in this directory (not recursive)
    files <- list.files(
      d,
      pattern = "(?i)\\.(tif|tiff|shp|gpkg|geojson)$",
      full.names = TRUE,
      recursive = FALSE
    )

    for (f in files) {
      ext  <- tolower(tools::file_ext(f))
      name <- tools::file_path_sans_ext(basename(f))

      # avoid name collisions within this folder
      if (name %in% names(out)) {
        name <- make.unique(c(names(out), name))[length(names(out)) + 1L]
      }

      obj <- switch(
        ext,
        tif     = tryCatch(terra::rast(f), error = function(e) { warning("Failed to read raster: ", f, " (", e$message, ")"); NULL }),
        tiff    = tryCatch(terra::rast(f), error = function(e) { warning("Failed to read raster: ", f, " (", e$message, ")"); NULL }),
        shp     = tryCatch(sf::st_read(f, quiet = TRUE), error = function(e) { warning("Failed to read vector: ", f, " (", e$message, ")"); NULL }),
        gpkg    = tryCatch(sf::st_read(f, quiet = TRUE), error = function(e) { warning("Failed to read vector: ", f, " (", e$message, ")"); NULL }),
        geojson = tryCatch(sf::st_read(f, quiet = TRUE), error = function(e) { warning("Failed to read vector: ", f, " (", e$message, ")"); NULL }),
        { warning("Unsupported extension: ", ext, " for file ", f); NULL }
      )

      if (!is.null(obj)) out[[name]] <- obj
    }

    # subdirectories (not recursive here; recurse manually)
    subdirs <- list.dirs(d, full.names = TRUE, recursive = FALSE)
    # guard against weird entries
    if (length(subdirs)) {
      is_dir <- dir.exists(subdirs)
      subdirs <- subdirs[is_dir]
    }

    for (sd in subdirs) {
      nm <- basename(sd)
      # avoid name collisions with files already added
      if (nm %in% names(out)) {
        nm <- make.unique(c(names(out), nm))[length(names(out)) + 1L]
      }
      out[[nm]] <- .read_dir(sd)  # <- correct recursive call
    }

    out
  }

  .read_dir(dir_path)
}
