#' Merge a list of 'sf' objects into a data frame of coordinates
#'
#' Combines multiple vector objects of class `sf` into a single data frame with X, Y, and group labels.
#'
#' @param shp_list A list of `sf` objects.
#' @param group_size Optional column name to include in output.
#' @returns A data.frame with columns `X`, `Y`, `group`, and optionally `group_size`.
#' @examples
#' # Create test data
#' vec1 <- df_to_shp(data.frame(long = c(1,2,2,4), lat = c(4,4,2,2)))
#' vec2 <- df_to_shp(data.frame(long = c(2,5,4,6), lat = c(4,4,2,2)))
#' vec_list <- list(vec1, vec2)
#'
#' # Convert vector list into data.frame
#' df <- merge_shp(vec_list)
#' df
#'
#' @export
merge_shp <- function(shp_list, group_size = NULL) {
  if (!is.list(shp_list)) {
    stop("Error: Input must be a list of shapefiles ('sf' objects).")
  }

  df <- data.frame(matrix(ncol = 3, nrow = 0))
  names(df) <- c("X", "Y", "group")

  for (i in seq_along(shp_list)) {
    shp <- shp_list[[i]]
    if (!inherits(shp, "sf")) {
      stop("Error: All list elements must be shapefiles ('sf' objects).")
    }

    group_name <- if (is.null(names(shp_list)) || names(shp_list)[i] == "") paste(i) else names(shp_list)[i]

    coords <- st_coordinates(shp)
    data_out <- data.frame(X = coords[, 1], Y = coords[, 2], group = group_name)

    if (!is.null(group_size) && group_size %in% colnames(shp)) {
      data_out[[group_size]] <- shp[[group_size]]
    }

    df <- rbind(df, data_out)
  }

  return(df)
}


#' Convert a data.frame to an sf object
#'
#' Builds an `sf` point layer from lon/lat columns of a data frame, guessing CRS via `guess_crs()`.
#'
#' @param df A data.frame with at least two numeric columns (lon, lat).
#' @return An `sf` object with geometry column and original attributes.
#' @importFrom sf st_as_sf
#' @examples
#' # example code
#' df <- data.frame(long = c(1,2,2,4), lat = c(4,4,2,2))
#' vec <- df_to_shp(df)
#' class(vec)
#' @export
df_to_shp <- function(df) {
  # 1) Input checks
  if (inherits(df, "sf")) {
    message("Warning: Input is a 'sf' (shapefile) objects. Nothing to do...")
    return(df)
  }
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame.")
  }
  if (ncol(df) < 2) {
    stop("Data frame must have at least two columns: longitude (1st) and latitude (2nd).")
  }

  lonlat <- as.matrix(df[, 1:2])
  if (!is.numeric(lonlat)) {
    stop("The first two columns must be numeric (lon/lat).")
  }

  # 2) Guess CRS from those coords
  guessed_epsg <- guess_crs(df)
  if (is.na(guessed_epsg)) {
    warning("Could not guess a valid geographic CRS. Setting CRS to NA.")
  }

  # 3) Build the sf object
  sf_obj <- st_as_sf(
    df,
    coords = names(df)[1:2],
    crs = guessed_epsg,
    remove = FALSE
  )

  return(sf_obj)
}


#' Split a data.frame or sf into a list of sf by group
#'
#' Splits data.frame rows into separate `sf` objects based on a grouping factor (column 3).
#'
#' @param df A data.frame or `sf` with at least two coordinate columns and optional group column.
#' @return A named list of `sf` objects, one per group level.
#' @examples
#' # Create test data
#' df <- data.frame(long = c(1,2,2,4,2,5,4,6),
#'                  lat = c(4,4,2,2,4,4,2,2),
#'                  species=rep(c("sp1", "sp2"), each=4))

#' #create vector list
#' vec_list <- df_to_list(df)
#' vec_list
#' @export
df_to_list <- function(df) {
  # Input validation
  if (inherits(df, "sf")) {
    shp_list <- list(df)
    names(shp_list) <- deparse(substitute(df))
    return(shp_list)
  }
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame.")
  }
  if (ncol(df) < 2) {
    stop("Data frame must have at least two columns: longitude, latitude")
  }

  # Extract the species/grouping variable (3rd column) and its levels
  group_var <- rep(deparse(substitute(df)), dim(df)[1])
  if (ncol(df) > 2) {
    group_var <- df[[3]]
  }
  grp_factor <- as.factor(group_var)
  levels_vec <- levels(grp_factor)

  # Split and convert
  shp_list <- lapply(levels_vec, function(level_value) {
    df_sub <- df[grp_factor == level_value, , drop = FALSE]
    df_to_shp(df_sub)
  })

  # Name list elements by the factor levels
  names(shp_list) <- levels_vec

  return(shp_list)
}
