#' Compute Summary Statistics for a Single HRA Run
#'
#' This function summarizes exposure (E), consequence (C), and risk (R) maps
#' for each stressor in a single Habitat Risk Assessment (HRA) run. It extracts
#' global minimum, maximum, and mean values for E, C, and raw risk layers, and
#' calculates the percentage of raster cells falling into each risk class
#' (None, Low, Medium, High). It also computes overall aggregated statistics
#' across all stressors.
#'
#' @param rr A named list of stressor results, where each element must contain
#'   'SpatRaster' objects for 'E_criteria', 'C_criteria',
#'   'Risk_map_raw', and 'Risk_map'.
#' @param total_risk_raw A 'SpatRaster' with overall raw risk values
#'   aggregated across stressors.
#' @param total_class A 'SpatRaster' with overall classified risk values
#'   (0 = None, 1 = Low, 2 = Medium, 3 = High).
#' @return A 'data.frame' with one row per stressor and an additional row
#'   with the overall statistics. The columns include:
#'   \itemize{
#'     \item 'STRESSOR': Stressor name.
#'     \item 'E_min, 'E_max', 'E_mean': Exposure statistics.
#'     \item 'C_min, 'C_max', 'C_mean': Consequence statistics.
#'     \item 'R_min, 'R_max', 'R_mean': Raw risk statistics.
#'     \item 'R%high, 'R%medium', 'R%low', 'R%None':
#'           Percentages of cells in each risk class.
#'   }
#' @details
#' This function is primarily intended as an internal utility for summarizing
#' outputs of the HRA workflow. It uses \code{terra::global()} to extract
#' summary values and \code{terra::freq()} to tabulate class frequencies.
#' @seealso \code{\link{hra}} for running a full HRA analysis.
#' @examples
#' \dontrun{
#' # Suppose rr is a list of results from hra4()
#' stats <- compute_summary_stats_single(
#'   rr = hra_result,
#'   total_risk_raw = hra_result$total_raw,
#'   total_class = hra_result$total_reclassified
#' )
#' head(stats)
#' }
#' @export
compute_summary_stats_single <- function(rr, total_risk_raw, total_class) {
  stressor_names <- setdiff(names(rr), c("total_raw", "total_reclassified", "total_hotspots_reclassified", "summary_stats"))
  rows <- lapply(stressor_names, function(st) {
    E <- rr[[st]]$E_criteria
    C <- rr[[st]]$C_criteria
    Rr <- rr[[st]]$Risk_map_raw
    Rc <- rr[[st]]$Risk_map
    gE <- terra::global(E,  c("min","max","mean"), na.rm=TRUE)
    gC <- terra::global(C,  c("min","max","mean"), na.rm=TRUE)
    gR <- terra::global(Rr, c("min","max","mean"), na.rm=TRUE)
    fq <- terra::freq(Rc)
    cnt <- c(`0`=0,`1`=0,`2`=0,`3`=0)
    if (!is.null(fq) && nrow(fq)) {
      vv <- as.character(fq[, "value"]); cnt[vv] <- fq[, "count"]
    }
    tot <- sum(cnt)

    data.frame(
      STRESSOR = st,
      E_min = as.numeric(gE[1,"min"]), E_max = as.numeric(gE[1,"max"]), E_mean = as.numeric(gE[1,"mean"]),
      C_min = as.numeric(gC[1,"min"]), C_max = as.numeric(gC[1,"max"]), C_mean = as.numeric(gC[1,"mean"]),
      R_min = as.numeric(gR[1,"min"]), R_max = as.numeric(gR[1,"max"]), R_mean = as.numeric(gR[1,"mean"]),
      `R%high` = if (tot) 100*cnt["3"]/tot else NA_real_,
      `R%medium` = if (tot) 100*cnt["2"]/tot else NA_real_,
      `R%low` = if (tot) 100*cnt["1"]/tot else NA_real_,
      `R%None` = if (tot) 100*cnt["0"]/tot else NA_real_
    )
  })
  per <- do.call(rbind, rows)

  fqt <- terra::freq(total_class)
  cntt <- c(`0`=0,`1`=0,`2`=0,`3`=0)
  if (!is.null(fqt) && nrow(fqt)) {
    vv <- as.character(fqt[, "value"]); cntt[vv] <- fqt[, "count"]
  }
  tot <- sum(cntt)

  overall <- data.frame(
    STRESSOR="(FROM ALL STRESSORS)",
    E_min=min(per$E_min,na.rm=TRUE), E_max=max(per$E_max,na.rm=TRUE), E_mean=mean(per$E_mean,na.rm=TRUE),
    C_min=min(per$C_min,na.rm=TRUE), C_max=max(per$C_max,na.rm=TRUE), C_mean=mean(per$C_mean,na.rm=TRUE),
    R_min=min(per$R_min,na.rm=TRUE), R_max=max(per$R_max,na.rm=TRUE), R_mean=mean(per$R_mean,na.rm=TRUE),
    `R%high` = if (tot) 100*cntt["3"]/tot else NA_real_,
    `R%medium` = if (tot) 100*cntt["2"]/tot else NA_real_,
    `R%low` = if (tot) 100*cntt["1"]/tot else NA_real_,
    `R%None` = if (tot) 100*cntt["0"]/tot else NA_real_
  )
  rbind(overall, per)
}
