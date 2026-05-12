#' Merge a list of sf layers
#'
#' Combines multiple `sf` objects stored in a list into a single object. By
#' default, the function performs a simple row-wise merge of the input layers,
#' keeping all geometries as independent features. Attribute columns that are not
#' shared among all input layers are preserved and filled with `NA` where absent.
#'
#' The function can also return a coordinate table or, if explicitly requested,
#' perform a spatial union of all geometries. The union option should be used with
#' care because `sf::st_union()` performs a topological dissolve and may be slow
#' or memory-intensive for large or complex geometries.
#'
#' All input layers are transformed to the CRS of the first element in
#' `shp_list`, when possible. If an input layer has no CRS, it is kept as-is and
#' a warning is issued.
#'
#' @param shp_list A non-empty list of `sf` objects.
#' @param output Character. Output type. One of `"sf"`, `"coords"`, or
#'   `"union"`. If `"sf"`, the input layers are combined into a single `sf`
#'   object using row binding. If `"coords"`, a `data.frame` with the coordinates
#'   of all geometries is returned. If `"union"`, all geometries are dissolved
#'   into a single union geometry using `sf::st_union()`.
#'
#' @return
#' Depending on `output`, returns:
#'
#' * `"sf"`: a single `sf` object containing all input features.
#' * `"coords"`: a `data.frame` with columns `X`, `Y`, and `group`.
#' * `"union"`: an `sf` object containing the dissolved union of all input
#'   geometries.
#'
#' @details
#' The `"sf"` option performs a simple merge and is usually the appropriate
#' choice when the goal is to combine species, stressor, or other spatial layers
#' into a single object without changing the geometries. Unlike `base::rbind()`,
#' this implementation uses `dplyr::bind_rows()`, allowing input layers to have
#' different attribute columns.
#'
#' The `"coords"` option extracts the coordinates of all input geometries using
#' `sf::st_coordinates()` and returns them as a regular `data.frame`. A `group`
#' column is added using the names of `shp_list`, or the list index when names
#' are not available.
#'
#' The `"union"` option applies `sf::st_union()` to the merged geometries. This
#' is a spatial/topological operation, not a simple merge, and can consume large
#' amounts of memory for large datasets.
#'
#' @importFrom sf st_crs st_transform st_geometry st_as_sf st_union st_coordinates
#' @importFrom dplyr bind_rows
#'
#' @examples
#' # Create test data
#' vec1 <- df_to_shp(data.frame(long = c(1, 2, 2, 4),
#'                              lat  = c(4, 4, 2, 2)))
#'
#' vec2 <- df_to_shp(data.frame(long = c(2, 5, 4, 6),
#'                              lat  = c(4, 4, 2, 2)))
#'
#' vec_list <- list(species1 = vec1, species2 = vec2)
#'
#' # Merge into a single sf object
#' shp <- merge_shp(vec_list, output = "sf")
#'
#' # Extract coordinates as a data.frame
#' coords <- merge_shp(vec_list, output = "coords")
#'
#' # Dissolve all geometries into a single union geometry
#' union_shp <- merge_shp(vec_list, output = "union")
#'
#' @export
merge_shp <- function(shp_list, output = c("sf", "coords", "union")) {
  output <- match.arg(output)

  if (!is.list(shp_list) || length(shp_list) == 0L) {
    stop("`shp_list` must be a non-empty list of sf objects.")
  }

  if (!all(vapply(shp_list, inherits, logical(1), what = "sf"))) {
    stop("All elements of `shp_list` must be sf objects.")
  }

  crs0 <- sf::st_crs(shp_list[[1]])

  shp_list <- lapply(shp_list, function(s) {
    crsS <- sf::st_crs(s)

    if (is.na(crs0)) return(s)

    if (is.na(crsS)) {
      warning("One input layer has no CRS. It was kept as-is.")
      return(s)
    }

    if (crsS != crs0) {
      s <- sf::st_transform(s, crs0)
    }

    s
  })

  nm <- names(shp_list)

  shp_list <- Map(function(s, i) {
    group_name <- if (!is.null(nm) && nzchar(nm[i])) nm[i] else as.character(i)
    s$group <- group_name
    s
  }, shp_list, seq_along(shp_list))

  if (output == "sf") {
    return(dplyr::bind_rows(shp_list))
  }

  if (output == "coords") {
    out_list <- lapply(shp_list, function(s) {
      coords <- sf::st_coordinates(s)

      if (nrow(coords) == nrow(s)) {
        group <- s$group
      } else if ("L1" %in% colnames(coords)) {
        group <- s$group[coords[, "L1"]]
      } else {
        group <- rep(unique(s$group), nrow(coords))
      }

      data.frame(
        X = coords[, 1],
        Y = coords[, 2],
        group = group,
        stringsAsFactors = FALSE
      )
    })

    out <- do.call(rbind, out_list)
    row.names(out) <- NULL
    return(out)
  }

  if (output == "union") {
    merged <- dplyr::bind_rows(shp_list)
    geom_union <- sf::st_union(sf::st_geometry(merged))
    return(sf::st_as_sf(geom_union))
  }
}


#' Convert a data.frame to an sf object
#'
#' Builds an `sf` POINT layer from longitude/latitude columns of a data frame.
#' If `crs` is not supplied, tries `guess_crs()`. If guessing fails,
#' falls back to EPSG:4326 only when coordinates look like degrees; otherwise uses `NA` CRS.
#'
#' @param df A data.frame with at least two numeric columns for coordinates, or an `sf` (returned as-is).
#' @param lon,lat Optional column names or indices for longitude and latitude.
#'   If omitted, will try common names (lon/long/longitude/x, lat/latitude/y), else the first two numeric columns.
#' @param crs Optional CRS (e.g., EPSG code like 4326, WKT, or an `sf::crs` object).
#' @param guess Logical; if `TRUE`, try `guess_crs()` when `crs` is not provided. Default `TRUE`.
#' @param keep_coords Logical; keep original lon/lat columns (`TRUE`, default) or drop them (`FALSE`).
#' @param drop_na Logical; drop rows with non-finite lon/lat (`TRUE`, default). If `FALSE`, will error on NA/Inf.
#' @param quiet Logical; suppress messages/warnings (default `TRUE`). Set to `FALSE` for informative messages.
#' @return An `sf` object with geometry column and original attributes.
#' @importFrom sf st_as_sf st_crs
#' @examples
#' df <- data.frame(long = c(1,2,2,4), lat = c(4,4,2,2))
#' vec <- df_to_shp(df)
#' @export
df_to_shp <- function(df,
                      lon = NULL, lat = NULL,
                      crs = NULL,
                      guess = TRUE,
                      keep_coords = TRUE,
                      drop_na = TRUE,
                      quiet = TRUE) {

  # Checks input
  if (inherits(df, "sf")) {
    if (!quiet) message("Input is already an 'sf' object; returning as-is.")
    return(df)
  }
  if (!is.data.frame(df)) stop("`df` must be a data.frame.")

  # Choose lon/lat columns
  nms <- names(df)
  if (is.null(lon) || is.null(lat)) {
    canon  <- tolower(nms)
    lon_ix <- which(canon %in% c("lon","long","longitude","x","longitud"))
    lat_ix <- which(canon %in% c("lat","latitude","y","latitud"))
    if (length(lon_ix) >= 1 && length(lat_ix) >= 1) {
      lon <- nms[lon_ix[1]]; lat <- nms[lat_ix[1]]
    } else {
      num_idx <- which(vapply(df, is.numeric, logical(1)))
      if (length(num_idx) < 2) stop("Could not find lon/lat: supply `lon`/`lat` or add two numeric columns.")
      lon <- nms[num_idx[1]]; lat <- nms[num_idx[2]]
      if (!quiet) message(sprintf("Using numeric columns '%s' (lon) and '%s' (lat); supply `lon`/`lat` to control.", lon, lat))
    }
  } else {
    if (is.numeric(lon)) lon <- nms[lon]
    if (is.numeric(lat)) lat <- nms[lat]
  }
  if (!all(c(lon, lat) %in% nms)) stop("`lon`/`lat` columns not found in `df`.")
  if (!is.numeric(df[[lon]]) || !is.numeric(df[[lat]])) stop("`lon` and `lat` must be numeric.")

  # Handle NA/Inf
  if (drop_na) {
    keep <- is.finite(df[[lon]]) & is.finite(df[[lat]])
    if (!all(keep)) df <- df[keep, , drop = FALSE]
  } else if (any(!is.finite(df[[lon]]) | !is.finite(df[[lat]]))) {
    stop("`lon`/`lat` contain NA/Inf; set `drop_na = TRUE` to drop them.")
  }

  # Determine CRS
  target_crs <- NULL
  if (!is.null(crs)) {
    target_crs <- sf::st_crs(crs)
  } else if (isTRUE(guess) && exists("guess_crs", mode = "function")) {
    guessed <- tryCatch(
      if (quiet) suppressMessages(guess_crs(df)) else guess_crs(df),
      error = function(e) NA
    )
    target_crs <- tryCatch(sf::st_crs(guessed), error = function(e) NA)
  }
  if (is.null(target_crs) || is.na(target_crs)) {
    rng_lon <- range(df[[lon]], na.rm = TRUE)
    rng_lat <- range(df[[lat]], na.rm = TRUE)
    if (all(rng_lon >= -180 & rng_lon <= 180) && all(rng_lat >= -90 & rng_lat <= 90)) {
      target_crs <- sf::st_crs(4326)
    } else {
      if (!quiet) message("CRS could not be determined; creating geometry with NA CRS.")
      target_crs <- NA
    }
  }

  # Build sf
  sf::st_as_sf(
    df,
    coords = c(lon, lat),
    crs = target_crs,
    remove = !keep_coords
  )
}


#' Split a data.frame or sf into a list of sf by group
#'
#' Splits rows into separate `sf` objects based on a grouping column.
#' If `df` is a data.frame, it's converted once with `df_to_shp()` and then split.
#'
#' @param df A data.frame (with lon/lat columns) or an `sf` object.
#' @param group Optional grouping column (name or index). If omitted:
#'   - for data.frame with ≥3 columns, uses the 3rd column;
#'   - for sf, returns a single-element list (no split).
#' @param drop_na_group Logical; if `FALSE`, rows with NA in `group` are returned
#'   as a `<NA>` element. Default `TRUE`.
#' @param ... Passed to `df_to_shp()` when `df` is a data.frame (e.g., `lon`, `lat`, `crs`).
#' @return A named list of `sf` objects, one per group level (or a single item if no grouping).
#' @examples
#' df <- data.frame(
#'   long = c(1,2,2,4, 2,5,4,6),
#'   lat  = c(4,4,2,2, 4,4,2,2),
#'   species = rep(c("sp1","sp2"), each = 4)
#' )
#' vec_list <- df_to_list(df, group = "species")
#' @export
df_to_list <- function(df, group = NULL, drop_na_group = TRUE, ...) {
  nm <- deparse(substitute(df))

  # Normalize to sf once
  if (inherits(df, "sf")) {
    sf_obj <- df
  } else if (is.data.frame(df)) {
    if (is.null(group) && ncol(df) >= 3L) group <- names(df)[3L]
    sf_obj <- df_to_shp(df, ...)
  } else {
    stop("`df` must be a data.frame or an sf object.")
  }

  # If no grouping requested/available, return single item
  if (is.null(group)) {
    lst <- list(sf_obj)
    names(lst) <- nm
    return(lst)
  }

  # Resolve group column name safely
  grp_name <- if (is.numeric(group)) names(sf_obj)[group] else group
  sf_col <- attr(sf_obj, "sf_column")
  if (!grp_name %in% names(sf_obj) || grp_name == sf_col) {
    stop("`group` must refer to a non-geometry column present in `df`/`sf`.")
  }

  g <- sf_obj[[grp_name]]

  # Split efficiently; control NA handling
  out <- split(sf_obj, g, drop = TRUE)
  if (!drop_na_group && anyNA(g)) {
    out[["<NA>"]] <- sf_obj[is.na(g), , drop = FALSE]
  }

  # If nothing to split (single level), still return a named list
  if (length(out) == 0L) {
    out <- list(sf_obj)
    names(out) <- nm
  }

  out
}


#' Compute the nesting depth of a list
#'
#' Determines the maximum number of nested list layers in an R object.
#'
#' @param x An R object. If it is a list (possibly nested), its maximum nesting depth is returned; otherwise 0.
#' @return An integer scalar: the depth of list nesting (0 for non‐lists, 1 for an empty list, etc.).
#' @examples
#' # Not a list
#' list_depth_base(42)
#'
#' # Empty list
#' list_depth_base(list())
#'
#' # Nested list
#' nested <- list(a = list(b = list(c = 1)))
#' list_depth_base(nested)
#'
#' @export
list_depth_base <- function(x) {
  if (!is.list(x) || length(x) == 0) {
    return(if (is.list(x)) 1L else 0L)
  }
  1L + max(vapply(x, list_depth_base, integer(1)))
}


#' Convert an object into a named list of `sf` objects
#'
#' This function standardizes different input types (`sf`, `data.frame`,
#' or list of those) into a named list of `sf` objects. It is useful for
#' preparing heterogeneous inputs before further spatial operations.
#'
#' If the input is:
#' - An `sf` object: it will be wrapped into a list with a name.
#' - A `data.frame`: it will be converted into an `sf` object using
#'   `df_to_shp()` or split by a grouping column if provided.
#' - A list: each element must be an `sf` object or a `data.frame`, which
#'   will be converted accordingly. Names are preserved when possible,
#'   otherwise they are auto-generated using `label_prefix`.
#'
#' @param obj An object of class `sf`, `data.frame`, or a list of those.
#' @param group Optional. A column name (string) used to split a
#'   `data.frame` into multiple `sf` layers. Ignored if `obj` is not a
#'   `data.frame`.
#' @param label_prefix Character string used as prefix for auto-generated
#'   names if input layers are unnamed. Default is `"layer"`.
#' @return A named list of `sf` objects.
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Example with a single sf object
#' nc <- st_read(system.file("shape/nc.shp", package="sf"), quiet = TRUE)
#' lst1 <- as_sf_list(nc)
#'
#' # Example with a data.frame (requires df_to_shp defined in package)
#' df <- data.frame(x = c(1,2), y = c(3,4), id = c("A","B"))
#' lst2 <- as_sf_list(df)
#'
#' # Example with a list of sf/data.frame
#' lst3 <- as_sf_list(list(nc1 = nc[1:10,], nc2 = nc[11:20,]))
#' }
#'
#' @export
as_sf_list <- function(obj, group = NULL, label_prefix = "layer") {
  if (inherits(obj, "sf")) {
    nm <- if (!is.null(attr(obj, "name")) && nzchar(attr(obj, "name"))) {
      attr(obj, "name")
    } else label_prefix
    return(setNames(list(obj), nm))
  }
  if (is.data.frame(obj)) {
    if (!is.null(group) && group %in% names(obj)) {
      sp <- split(obj, obj[[group]], drop = TRUE)
      lst <- lapply(sp, df_to_shp)
      return(lst)
    } else {
      return(df_to_list(obj))
    }
  }
  if (is.list(obj)) {
    nms <- names(obj)
    lst <- lapply(obj, function(el) {
      if (inherits(el, "sf")) return(el)
      if (is.data.frame(el))  return(df_to_shp(el))
      stop("List elements must be `sf` or data.frame.")
    })
    if (is.null(nms) || any(!nzchar(nms))) {
      if (!is.null(nms)) {
        nms[nchar(nms) == 0] <- paste0(label_prefix, seq_len(sum(nchar(nms) == 0)))
        names(lst) <- nms
      } else {
        names(lst) <- paste0(label_prefix, seq_along(lst))
      }
    }
    return(lst)
  }
  stop("Input must be an `sf`, data.frame, or a list of those.")
}


#' Rescale a SpatRaster to a New Numeric Range
#'
#' Linearly rescales the values of a `terra::SpatRaster` to a new user-defined
#' numeric range. The output is returned as a `SpatRaster` with values rounded
#' to four decimal places.
#'
#' @param x A `terra::SpatRaster` object containing the values to be rescaled.
#' @param new_range A numeric vector of length 2 giving the new minimum and
#'   maximum values for rescaling (`c(new_min, new_max)`).
#'
#' @return A `terra::SpatRaster` object with values rescaled to the specified
#'   range.
#'
#' @details
#' The function performs linear rescaling using:
#' \deqn{
#'   x' = a + \frac{(x - \min(x)) (b - a)}{\max(x) - \min(x)}
#' }
#' where \eqn{a} and \eqn{b} are the user-defined minimum and maximum.
#'
#' If the raster has no variation (`min == max`), the function throws an error.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r <- rast(matrix(runif(100), 10, 10))
#' r2 <- scale_vars(r, c(0, 1))
#' }
#'
#' @importFrom terra minmax
#' @export
scale_vars <- function(x, new_range, quiet = TRUE) {

  if (!is.numeric(new_range) || length(new_range) != 2L) {
    if (!quiet) stop("Error: new_range must be a numeric vector of length 2 (min, max).")
  }
  if (new_range[1] >= new_range[2]) {
    if (!quiet) stop("Error: minimum value for rescaling is not valid (must be < max).")
  }

  mm <- terra::minmax(x)
  old_min <- mm[1]
  old_max <- mm[2]

  if (old_min == old_max) {
    if (!quiet) stop("Error: cannot rescale a raster with constant values (min == max).")
  }

  scaled_x <- new_range[1] +
    ((x - old_min) * (new_range[2] - new_range[1])) / (old_max - old_min)

  return(round(scaled_x, 4))
}


#' Apply Value Rescaling to a List of SpatRasters
#'
#' Applies linear rescaling to each element of a list of `SpatRaster` objects,
#' using a user-defined numeric range. Optionally rounds integer rasters after
#' scaling.
#'
#' @param rlist A list of `terra::SpatRaster` objects to be rescaled.
#' @param new_range A numeric vector of length 2 defining the new value range
#'   for rescaling (`c(new_min, new_max)`). If `NULL`, no rescaling is applied.
#' @param cat Logical. If `TRUE`, any raster that is internally stored as an
#'   integer type (`terra::is.int()`) will be rounded after scaling.
#' @param quiet Logical. If `TRUE`, suppresses progress messages.
#'
#' @return A list of `terra::SpatRaster` objects, rescaled where appropriate.
#'
#' @details
#' For each raster:
#' * Its current min/max is compared with `new_range`.
#' * If identical, no rescaling is performed.
#' * Otherwise, the raster is passed to `scale_vars()`.
#'
#' The function preserves the structure of the input list.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r1 <- rast(matrix(runif(100), 10, 10))
#' r2 <- rast(matrix(runif(100), 10, 10))
#' lst <- list(a = r1, b = r2)
#'
#' scaled <- apply_scale(lst, c(0, 1))
#' }
#'
#' @importFrom terra minmax is.int
#' @export
apply_scale <- function(rlist, new_range, cat = FALSE, quiet = TRUE) {
  x <- rlist

  if (!is.null(new_range)) {
    for (i in seq_along(x)) {
      r <- x[[i]]
      mm <- terra::minmax(r)
      old_min <- mm[1]
      old_max <- mm[2]

      if (old_min == new_range[1] && old_max == new_range[2]) {
        if (!quiet) message(names(x)[i], " range and target range are equal. Nothing to do...")
      } else {
        x[[i]] <- scale_vars(r, new_range)
      }

      if (cat && terra::is.int(x[[i]])) {
        x[[i]] <- round(x[[i]])
      }
    }
  }
  return(x)
}


#' Check HRA criteria table format
#'
#' Checks whether a criteria table follows the expected format used by the
#' Habitat Risk Assessment (HRA) model implemented in *InVEST*. The function
#' verifies that all required columns are present in a single criteria
#' `data.frame`, or in each element of a list of criteria tables.
#'
#' @description
#' This function is used internally to validate the criteria table supplied to
#' HRA routines. In single-species applications, the input should be a
#' `data.frame`. In ecosystem-level applications, where criteria may differ
#' among species, the input should be a named list of `data.frame` objects.
#'
#' @param df_or_list A `data.frame` containing the HRA criteria table, or a
#'   list of `data.frame` objects. Each table must contain the standard HRA
#'   criteria columns: `STRESSOR`, `ATTRIBUTES`, `RATING`, `DQ`, `WEIGHT`, and
#'   `E/C`.
#'
#' @details
#' The expected input follows the standard criteria table structure used in the
#' *InVEST* Habitat Risk Assessment model. Each row represents a criterion used
#' to estimate either exposure (`E`) or consequence (`C`) scores for one or more
#' stressors.
#'
#' The required columns are:
#'
#' \itemize{
#'   \item `STRESSOR`: name of the stressor associated with the criterion.
#'   \item `ATTRIBUTES`: name of the criterion or spatial attribute.
#'   \item `RATING`: criterion rating. This may be a numeric value or `NA` when
#'   the rating is supplied by a spatial raster layer.
#'   \item `DQ`: data quality score.
#'   \item `WEIGHT`: criterion weight.
#'   \item `E/C`: whether the criterion contributes to exposure (`E`) or
#'   consequence (`C`).
#' }
#'
#' If `df_or_list` is a single `data.frame`, the function checks whether all
#' required columns are present and returns the original object. If `df_or_list`
#' is a list, the function checks each element independently and returns the
#' validated list.
#'
#' @return
#' Returns the original `data.frame` or list of `data.frame` objects if all
#' required columns are present. Otherwise, an error is raised.
#'
#' @examples
#' path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
#' df <- read.csv(path)
#' crit_list <- criteria_reshape(df)
#' crit_list
#'
#' @export
check_criteria <- function(df_or_list) {
  req <- c("STRESSOR","ATTRIBUTES","RATING","DQ","WEIGHT","E/C")
  if (is.data.frame(df_or_list)) {
    if (!all(req %in% names(df_or_list))) stop("Criteria must contain: ", paste(req, collapse=", "))
    return(df_or_list)
  }
  if (is.list(df_or_list)) {
    lapply(df_or_list, function(d) {
      if (!is.data.frame(d) || !all(req %in% names(d))) {
        stop("Each criteria element must contain: ", paste(req, collapse=", "))
      }
      d
    })
  } else stop("'criteria' must be a data.frame (single) or list (ecosystem).")
}



