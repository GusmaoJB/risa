#' Merge a list of sf layers
#'
#' @param shp_list list of sf objects
#' @param group_size optional column name to carry into points mode (old behavior)
#' @param return one of c("sf_union","points") — default "sf_union" for risa_prep()
#' @return if return = "sf_union": an sf with dissolved geometry (union of inputs)
#'         if return = "points":   a data.frame of XY (old behavior)
#' @examples
#' # Create test data
#' vec1 <- df_to_shp(data.frame(long = c(1,2,2,4), lat = c(4,4,2,2)))
#' vec2 <- df_to_shp(data.frame(long = c(2,5,4,6), lat = c(4,4,2,2)))
#' vec_list <- list(vec1, vec2)
#'
#' # Convert vector list into data.frame
#' df <- merge_shp(vec_list)
#' @export
merge_shp <- function(shp_list, group_size = NULL, return = c("sf_union","points")) {
  return <- match.arg(return)

  if (!is.list(shp_list) || length(shp_list) == 0L) {
    stop("`shp_list` must be a non-empty list of sf objects.")
  }
  if (!all(vapply(shp_list, inherits, logical(1), what = "sf"))) {
    stop("All elements of `shp_list` must be sf objects.")
  }

  # --------- New default: union to a single sf geometry (for risa_prep) -------
  if (return == "sf_union") {
    # Harmonize CRS to the first layer (safer than mixing)
    crs0 <- sf::st_crs(shp_list[[1]])
    shp_list <- lapply(shp_list, function(s) {
      crsS <- sf::st_crs(s)
      same <- !is.na(crsS) && !is.na(crs0) && identical(crsS$wkt, crs0$wkt)
      if (!same) sf::st_transform(s, crs0) else s
    })

    # Union all geometries; handle single-element list
    geoms <- lapply(shp_list, sf::st_geometry)
    geom_union <- if (length(geoms) == 1L) geoms[[1]] else Reduce(sf::st_union, geoms)

    # Wrap back to sf (CRS preserved)
    return(sf::st_as_sf(geom_union))
  }

  # Warn if CRSs differ; coordinates will be mixed as-is
  crs_keys <- vapply(shp_list, function(s) {
    crs <- sf::st_crs(s)
    if (!is.null(crs$epsg)) paste0("epsg:", crs$epsg) else if (!is.null(crs$wkt)) crs$wkt else "NA"
  }, character(1))
  if (length(unique(crs_keys)) > 1L) {
    warning("Input layers have different CRS; coordinates are kept as-is. Consider transforming beforehand.")
  }

  out_list <- vector("list", length(shp_list))
  nm_list  <- names(shp_list)

  for (i in seq_along(shp_list)) {
    shp <- shp_list[[i]]
    group_name <- if (!is.null(nm_list) && nzchar(nm_list[i])) nm_list[i] else as.character(i)

    # Coordinates for all vertices
    coords <- sf::st_coordinates(shp)

    df <- data.frame(
      X = coords[, 1],
      Y = coords[, 2],
      group = rep(group_name, nrow(coords)),
      stringsAsFactors = FALSE
    )

    if (!is.null(group_size) && group_size %in% names(shp)) {
      n_per_feat <- sf::st_npoints(sf::st_geometry(shp), by_feature = TRUE)
      df[[group_size]] <- rep(shp[[group_size]], times = n_per_feat)
    }

    out_list[[i]] <- df
  }

  out <- do.call(rbind, out_list)
  row.names(out) <- NULL
  out
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

  # Early return for sf
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
    # Default grouping for data.frame: 3rd column if not provided
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
  sf_col   <- attr(sf_obj, "sf_column")
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
