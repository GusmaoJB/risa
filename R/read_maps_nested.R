#' Read spatial files into a nested list mirroring directory structure
#'
#' This function scans a root directory (non‐recursively at each level) for raster
#' files (`.tif`, `.tiff`) and shapefiles (`.shp`), reads them as `terra::SpatRaster`
#' or `sf` objects, respectively, and builds a nested R list whose names correspond
#' to folder and file basenames (without extensions).
#'
#' @param dir_path `character(1)`. Path to the root directory containing your spatial data.
#'   Subdirectories will be traversed recursively and represented as nested list elements.
#'
#' @return A `list` of arbitrary depth:
#'   - Intermediate elements are named by subdirectory names and contain further lists.
#'   - Leaf elements are named by file basenames and contain `SpatRaster` or `sf` objects.
#'
#' @importFrom terra rast
#' @importFrom sf st_read
#' @export
#'
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
read_maps_nested <- function(dir_path) {

  files <- list.files(dir_path,
                      pattern = "\\.(tif|tiff|shp)$",
                      full.names = TRUE,
                      recursive = FALSE,
                      ignore.case = TRUE)
  out <- list()

  for (f in files) {
    ext  <- tolower(tools::file_ext(f))
    name <- tools::file_path_sans_ext(basename(f))
    obj  <- switch(ext,
                   tif = terra::rast(f),
                   tiff = terra::rast(f),
                   shp = sf::st_read(f, quiet = TRUE),
                   stop("This extension is not supported: ", ext)
    )
    out[[name]] <- obj
  }

  subdirs <- list.dirs(dir_path, full.names = TRUE, recursive = FALSE)
  for (d in subdirs) {
    folder_name <- basename(d)
    out[[folder_name]] <- read_spatial_nested(d)
  }
  out
}
