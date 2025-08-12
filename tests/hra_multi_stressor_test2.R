#' Generate summary statistics for risaHRA outputs
#' @param list an object of class "risaHRA"
#' @param verbose logical; print progress messages?
#' @return data.frame with per-stressor (and ecosystem) stats
#' @importFrom terra global freq
get_stats <- function(list, verbose = FALSE) {
  if (!inherits(list, "risaHRA")) {
    stop("Input must be a 'risaHRA' object.")
  }

  log_msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  # helper: percentages of classes 0/1/2/3 (None/low/medium/high)
  class_pcts <- function(r) {
    if (is.null(r)) {
      return(setNames(rep(NA_real_, 4), c("R%high","R%medium","R%low","R%None")))
    }
    f <- terra::freq(r, digits = 0)
    f <- as.data.frame(f)
    f <- f[!is.na(f$value), , drop = FALSE]
    tot <- sum(f$count)
    getc <- function(v) {
      i <- which(f$value == v)
      if (length(i)) f$count[i] else 0L
    }
    if (tot == 0) {
      setNames(rep(NA_real_, 4), c("R%high","R%medium","R%low","R%None"))
    } else {
      setNames(100 * c(getc(3), getc(2), getc(1), getc(0)) / tot,
               c("R%high","R%medium","R%low","R%None"))
    }
  }

  # helper: stats for a single stressor sublist
  stressor_stats <- function(st_name, s) {
    # Expect single-layer rasters for each metric
    E  <- s$E_criteria
    C  <- s$C_criteria
    R  <- s$Risk_map_raw
    RR <- s$Risk_map

    e <- as.numeric(terra::global(E,  c("min","max","mean"), na.rm = TRUE)[1, ])
    c <- as.numeric(terra::global(C,  c("min","max","mean"), na.rm = TRUE)[1, ])
    r <- as.numeric(terra::global(R,  c("min","max","mean"), na.rm = TRUE)[1, ])
    p <- class_pcts(RR)

    data.frame(
      STRESSOR = st_name,
      E_min = e[1], E_max = e[2], E_mean = e[3],
      C_min = c[1], C_max = c[2], C_mean = c[3],
      R_min = r[1], R_max = r[2], R_mean = r[3],
      t(as.data.frame.list(p)),
      check.names = FALSE
    )
  }

  # helper: all stressors for a species-level list
  get_stressor_stats <- function(sublist) {
    # stressors are the list elements that themselves are lists
    stressor_names <- names(sublist)[vapply(sublist, is.list, logical(1))]
    log_msg("The list has %d stressors.", length(stressor_names))

    out <- do.call(
      rbind,
      lapply(stressor_names, function(nm) {
        log_msg("Calculating stats for stressor %s", nm)
        stressor_stats(nm, sublist[[nm]])
      })
    )

    # "FROM ALL STRESSORS" row: E/C/R from per-stressor rows; R% from total reclass raster
    total_reclass <- sublist[["total"]]  # expected to be the reclassified total risk
    p_tot <- class_pcts(total_reclass)

    all_row <- data.frame(
      STRESSOR = "(FROM ALL STRESSORS)",
      E_min  = min(out$E_min,  na.rm = TRUE),
      E_max  = max(out$E_max,  na.rm = TRUE),
      E_mean = mean(out$E_mean, na.rm = TRUE),
      C_min  = min(out$C_min,  na.rm = TRUE),
      C_max  = max(out$C_max,  na.rm = TRUE),
      C_mean = mean(out$C_mean, na.rm = TRUE),
      R_min  = min(out$R_min,  na.rm = TRUE),
      R_max  = max(out$R_max,  na.rm = TRUE),
      R_mean = mean(out$R_mean, na.rm = TRUE),
      t(as.data.frame.list(p_tot)),
      check.names = FALSE
    )

    rbind(all_row, out)
  }

  if ("ecosys_risk_raw" %in% names(list)) {
    eco_risk_index <- which(names(list) == "ecosys_risk_raw")
    log_msg("Calculating summary statistics for %d species.", eco_risk_index - 1L)

    # species blocks are before ecosystem layers
    sp_idxs <- seq_len(eco_risk_index - 1L)
    sp_out <- do.call(
      rbind,
      lapply(sp_idxs, function(i) {
        sp_name <- names(list)[i]
        cbind.data.frame(SPECIES = sp_name, get_stressor_stats(list[[i]]),
                         check.names = FALSE)
      })
    )

    log_msg("Calculating general ecosystem risk stats.")
    eco_raw    <- list[[eco_risk_index]]
    eco_recl   <- list[[eco_risk_index + 1L]]

    r <- as.numeric(terra::global(eco_raw, c("min","max","mean"), na.rm = TRUE)[1, ])
    p <- class_pcts(eco_recl)

    eco_row <- data.frame(
      SPECIES = "ECOSYSTEM",
      STRESSOR = "(FROM ALL STRESSORS)",
      E_min = NA_real_, E_max = NA_real_, E_mean = NA_real_,
      C_min = NA_real_, C_max = NA_real_, C_mean = NA_real_,
      R_min = r[1], R_max = r[2], R_mean = r[3],
      t(as.data.frame.list(p)),
      check.names = FALSE
    )

    rbind(sp_out, eco_row)
  } else {
    log_msg("Calculating summary statistics for one species.")
    get_stressor_stats(list)
  }
}
