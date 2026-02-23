#' Generate reclassification matrix for raster values
#'
#' Builds a matrix mapping value ranges to discrete class codes,
#' optionally excluding lowest 5 percent.
#'
#' @param raster A `terra::SpatRaster` object.
#' @param n_classes Integer 1–10 number of classes.
#' @param exclude_lowest Logical; if `TRUE`, values below 5 percent of max are mapped to `NA`.
#' @param lowest_prop Numeric; Defines the cutoff for the lowest proportion  (default `0.05`).
#' @returns A numeric matrix with columns `from`, `to`, and `class` for use with `terra::classify()`.
#' @importFrom terra global
#' @examples
#' # Example
#' r <- terra::rast(nrows=5, ncols=5, vals=c(0.01,2,3,4,5))
#' reclass_matrix(r)
#' @export
reclass_matrix <- function(raster, n_classes = 3, exclude_lowest = TRUE, lowest_prop = 0.05) {
  # Input checks
  if (!inherits(raster, "SpatRaster")) stop("`raster` must be a terra::SpatRaster.")
  if (!is.numeric(n_classes) || n_classes %% 1 != 0 || n_classes < 1 || n_classes > 10) {
    stop("`n_classes` must be an integer between 1 and 10.")
  }
  if (!is.logical(exclude_lowest) || length(exclude_lowest) != 1L) {
    stop("`exclude_lowest` must be a single logical (TRUE/FALSE).")
  }

  # Min/Max without loading all values
  mm <- terra::global(raster, fun = c("min", "max"), na.rm = TRUE)
  min_val <- as.numeric(mm[1, "min"])
  max_val <- as.numeric(mm[1, "max"])
  if (!is.finite(min_val) || !is.finite(max_val)) {
    stop("Raster contains no finite values (all NA/Inf).")
  }

  # Cutoff for the lowest proportion
  val_05 <- lowest_prop * max_val

  # Lower break for the classed bins
  lower_break <- if (exclude_lowest) val_05 else min_val
  if (lower_break > max_val) lower_break <- min_val

  # Equal-width breaks -> n_classes bins
  breaks <- seq(lower_break, max_val, length.out = n_classes + 1)

  # Build bins: make from/to the same length (n_classes)
  from <- breaks[-length(breaks)]
  to <- breaks[-1]
  to[length(to)] <- Inf  # last bin open-ended
  cls <- seq_len(n_classes)

  mat <- cbind(from = from, to = to, class = cls)

  # Prepend NA bin for excluded-lowest portion
  if (exclude_lowest) {
    mat <- rbind(c(from = -Inf, to = val_05, class = NA_real_), mat)
  }

  storage.mode(mat) <- "double"
  mat
}
