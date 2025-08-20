

sum_weighted <- function(df_map) {
  if (!nrow(df_map)) return(zero_r)
  parts <- vector("list", nrow(df_map))
  E_coefs <- list()
  C_coefs <- list()
  if (decay %in% c("linear","exponential","polynomial_2nd","polynomial_3rd","polynomial_4th") &&
      !is.null(buffer_m[stressor])) {
    for (i in 1:nrow(df_map)) {
      att <- df_map$ATTRIBUTES[i]
      crit_type <-  df_map[["E/C"]]
      r   <- rlist[[stressor]][[att]]
      if (is.null(r) || !inherits(r,"SpatRaster")) {
        stop("Missing/invalid raster for '", stressor, "' / attribute '", att, "'.")
      }
      decay_coef <- decay_coeffs(r, decay, buffer_m[stressor])

      if (identical(ec_col[i], "E")) {
        E_coefs[[att]] <- decay_coef
      } else if (identical(ec_col[i], "C")) {
        C_coefs[[att]] <- decay_coef
      }

      start_val <- as.numeric(terra::global(r, "min", na.rm = TRUE)[[1]])
      r <- get_decay_map(r, decay_coef, start_val)
      r <- .align_to(r, sp_distr, categorical = TRUE)
      r <- .align_to(r, sp_distr, categorical = TRUE)
      parts[[i]] <- r / (df_map$DQ[i] * df_map$WEIGHT[i])
    }
    list(weighted_result = Reduce(`+`, parts))

  }
  for (i in 1:nrow(df_map)) {
    att <- df_map$ATTRIBUTES[i]
    r   <- rlist[[stressor]][[att]]
    if (is.null(r) || !inherits(r,"SpatRaster")) {
      stop("Missing/invalid raster for '", stressor, "' / attribute '", att, "'.")
    }
    r <- .align_to(r, sp_distr, categorical = TRUE)
    parts[[i]] <- r / (df_map$DQ[i] * df_map$WEIGHT[i])
  }
  Reduce(`+`, parts)
}



