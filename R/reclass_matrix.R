#' Generate reclassification matrix for raster values
#'
#' Builds a matrix mapping value ranges to discrete class codes,
#' optionally excluding lowest 5 percent.
#'
#' @param raster A `SpatRaster` object from `terra`.
#' @param n_classes Integer 1–10 number of classes.
#' @param exclude_lowest Logical; if `TRUE`, values below 5 percent of max are mapped to `NA`.
#' @returns A numeric matrix with columns `from`, `to`, and `class` for use with `terra::classify()`.
#' @importFrom terra values
#' @examples
#' # Creating example data
#' m <- matrix(rep(c(0.01,2,3,4,5), 5), ncol = 5, byrow = TRUE)
#' coords <- data.frame(long = c(1,2,3,4,5), lat = c(1,2,3,4,5))
#' xres <- diff(coords$long)[1]
#' yres <- diff(coords$lat)[1]
#' e <- ext(min(coords$long) - xres/2,
#'          max(coords$long) + xres/2,
#'          min(coords$lat) - yres/2,
#'          max(coords$lat) + yres/2)
#' r <- rast(nrows = nrow(m),
#'           ncols = ncol(m),
#'           ext = e, crs = "")
#' values(r) <- as.vector(m)
#' plot(r)
#'
#' # Create reclassification matrix
#' rec_mt <- reclass_matrix(r) # exclude lowest values
#' rec_mt
#' rec_mt_all <- reclass_matrix(r, exclude_lowest = FALSE) # with all values
#' rec_mt_all
#' @export
reclass_matrix <- function(raster, n_classes = 3, exclude_lowest = TRUE, custom_max = NULL) {
  # 1. input checks
  if (!is.numeric(n_classes) ||
      n_classes %% 1 != 0 ||
      n_classes < 1 ||
      n_classes > 10) {
    stop("n_classes must be an integer between 1 and 10.")
  }

  # 2. pull out the values
  vals    <- values(raster, na.rm = TRUE)
  min_val <- min(vals, na.rm = TRUE)
  max_val <- ifelse(custom_max = NULL, max(vals, na.rm = TRUE), custom_max)
  val_05  <- 0.05 * max_val

  # 3. set the break bounds
  if (exclude_lowest) {
    lower_break <- val_05
  } else {
    lower_break <- min_val
  }
  breaks <- seq(lower_break, max_val, length.out = n_classes + 1)

  # 4. build the matrix
  mat <- NULL

  # If we want to exclude the lowest 5%, map them to NA
  if (exclude_lowest) {
    mat <- matrix(c(-Inf, val_05, NA), ncol = 3, byrow = TRUE)
  }

  # 4b. now append the n_classes bins
  for (i in seq_len(n_classes)) {
    from  <- breaks[i]
    to    <- if (i < length(breaks)) breaks[i + 1] else Inf
    cls   <- i
    mat   <- rbind(mat, c(from, to, cls))
  }

  colnames(mat) <- c("from", "to", "class")
  return(mat)
}
