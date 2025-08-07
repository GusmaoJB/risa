hra <- function(raster_list, species_distr, criteria, equation = c("euclidean", "multiplicative"), r_max = 3, n_overlap = NULL) {

  if (all(!names(criteria) %in% c("STRESSOR", "ATTRIBUTES", "RATING", "DQ", "WEIGHT", "E/C"))){
    stop("Criteria table must have the following columns: 'STRESSOR', 'ATTRIBUTES', 'RATING', 'DQ', 'WEIGHT', and 'E/C'")
  }

  if (r_max > 10 | r_max < 1) {
    stop("Maximum criteria score must be between 1 and 10")
  }

  equation <- match.arg(equation)
  stressors <- unique(criteria$STRESSOR)
  stressors <- stressors[!stressors %in% c(NA, "", "NA")]

  if (is.null(n_overlap)) {
    n_overlap <- length(stressors)
  }

  if (all(!names(raster_list) == stressors)) {
    stop("The list names in 'raster_list' must match the stressor names described in the 'criteria' table")
  }

  # Initializing list of risk results
  risk_results <- list()

  # Reformatting species distribution raster
  sp_distr <- species_distr
  sp_presence <- terra::ifel(!is.na(sp_distr), 1, NA)

  for (stressor in stressors) {
    crit_stressor <- criteria[criteria$STRESSOR %in% c(NA, "NA", "", stressor),]

    # Split consequence (C) and exposure (E)
    C_df <- crit_stressor[crit_stressor$`E/C` == "C",]
    E_df <- crit_stressor[crit_stressor$`E/C` == "E",]

    C_df$DQ <- as.integer(C_df$DQ)
    C_df$WEIGHT <- as.integer(C_df$WEIGHT)

    E_df$DQ <- as.integer(E_df$DQ)
    E_df$WEIGHT <- as.integer(E_df$WEIGHT)

    # Select constant criteria (with rating), and divide value by DQ * Weight
    E_const <- E_df[!is.na(E_df$RATING),]
    C_const <- C_df[!is.na(C_df$RATING),]

    E_numer_const <- sum(E_const$RATING / (E_const$DQ * E_const$WEIGHT))
    C_numer_const <- sum(C_const$RATING / (C_const$DQ * C_const$WEIGHT))

    # Select mapped criteria (RATING == NA), then sum over raster/(DQ*WEIGHT)
    E_mapped <- E_df[is.na(E_df$RATING),]
    C_mapped <- C_df[is.na(C_df$RATING),]

    # Now compute E and C criteria blocks
    E_numer_rast <- 0
    for (crit in seq_len(nrow(E_mapped))) {
      raster_mapping_E <- raster_list[[stressor]][[E_mapped$ATTRIBUTES[crit]]]
      E_numer_rast <- E_numer_rast + (raster_mapping_E / (E_mapped$DQ[crit] * E_mapped$WEIGHT[crit]))
    }

    C_numer_rast <- 0
    for (crit in seq_len(nrow(C_mapped))) {
      raster_mapping_C <- raster_list[[stressor]][[C_mapped$ATTRIBUTES[crit]]]
      C_numer_rast <- C_numer_rast + (raster_mapping_C / (C_mapped$DQ[crit] * C_mapped$WEIGHT[crit]))
    }

    # Define denominator (is the same for all cells)
    E_denom <- sum(1/(E_df$DQ * E_df$WEIGHT))
    C_denom <- sum(1/(C_df$DQ * C_df$WEIGHT))

    # Final score raster
    E_score_raster <- (E_numer_const + E_numer_rast) / E_denom
    C_score_raster <- (C_numer_const + C_numer_rast) / C_denom

    # Converting species distributions with zeros (necessary for outputs)
    sp_distr_zeros <- terra::ifel(!is.na(sp_distr), 0, NA)

    # Setting layer with zero E and C criteria
    E_score_raster_map <- terra::mosaic(E_score_raster, sp_distr_zeros, fun="first")
    C_score_raster_map <- terra::mosaic(terra::mask(C_score_raster, sp_distr), sp_distr_zeros, fun="first")

    # Compute risk scores
    if (equation == "multiplicative") {
      m_jkl <- r_max^2
      risk_scores <- C_score_raster * E_score_raster * sp_presence
    } else if (equation == "euclidean") {
      m_jkl <- sqrt(2 * (r_max - 1)^2)
      risk_scores <- sqrt((E_score_raster - 1)^2 + (C_score_raster - 1)^2) * sp_presence
    }

    # Include layer with risk zero into the risk score raster
    risk_scores <- terra::mosaic(risk_scores, sp_distr_zeros, fun="first")

    # Classify using InVEST criteria
    risk_classified <- terra::ifel(risk_scores == 0, 0,
                                   terra::ifel(risk_scores < (1/3)*m_jkl, 1,
                                               terra::ifel(risk_scores < (2/3)*m_jkl, 2, 3)))

    risk_results[[stressor]] <- list(
      E_criteria = E_score_raster_map,
      C_criteria = C_score_raster_map,
      Risk_map_raw = risk_scores,
      Risk_map = risk_classified
    )
  }

  # Calculate cumulative risk (sum across all stressors)
  total_risk <- Reduce("+", lapply(risk_results, function(x) x$Risk_map_raw))
  total_risk_classified <- terra::ifel(total_risk == 0, 0,
                                    terra::ifel(total_risk < (1/3)*m_jkl*n_overlap, 1,
                                                terra::ifel(total_risk < (2/3)*m_jkl*n_overlap, 2, 3)))

  risk_results$total_raw <- total_risk
  risk_results$total <- total_risk_classified

  return(risk_results)
}


# Generate summary statistics
get_stats <- function(list) {
  output_df <- data.frame()
  for (stressor in names(list)) {
    if (is.list(list[[stressor]])){
      rasters <- c(list[[stressor]]$E_criteria, list[[stressor]]$C_criteria,
                   list[[stressor]]$Risk_map_raw, list[[stressor]]$Risk_map)
      names(rasters) <- c("E", "C", "R", "R_reclass")
      df <- as.data.frame(rasters, na.rm = FALSE)
      total_cells <- sum(df$R_reclass >= 0, na.rm=TRUE)
      stats <- data.frame(
        STRESSOR = stressor,
        E_min   = min(df$E,   na.rm = TRUE),
        E_max   = max(df$E,   na.rm = TRUE),
        E_mean  = mean(df$E,  na.rm = TRUE),
        C_min   = min(df$C,   na.rm = TRUE),
        C_max   = max(df$C,   na.rm = TRUE),
        C_mean  = mean(df$C,  na.rm = TRUE),
        R_min   = min(df$R,   na.rm = TRUE),
        R_max   = max(df$R,   na.rm = TRUE),
        R_mean  = mean(df$R,  na.rm = TRUE),
        `R%high` = sum(df$R_reclass == 3, na.rm = TRUE) / total_cells * 100,
        `R%medium` = sum(df$R_reclass == 2, na.rm = TRUE) / total_cells * 100,
        `R%low` = sum(df$R_reclass == 1, na.rm = TRUE) / total_cells * 100,
        `R%None` = sum(df$R_reclass == 0, na.rm = TRUE) / total_cells * 100)
      output_df <- rbind.data.frame(output_df, stats)
    }
  }
  total_df <- as.data.frame(list$total, na.rm = FALSE)
  all_stressors_df <- rbind.data.frame(
    data.frame(
      STRESSOR = "(FROM ALL STRESSORS)",
      E_min   = min(output_df$E_min,   na.rm = TRUE),
      E_max   = max(output_df$E_max,   na.rm = TRUE),
      E_mean  = mean(output_df$E_mean,  na.rm = TRUE),
      C_min   = min(output_df$C_min,   na.rm = TRUE),
      C_max   = max(output_df$C_max,   na.rm = TRUE),
      C_mean  = mean(output_df$C_mean,  na.rm = TRUE),
      R_min   = min(output_df$R_min,   na.rm = TRUE),
      R_max   = max(output_df$R_max,   na.rm = TRUE),
      R_mean  = mean(output_df$R_mean,  na.rm = TRUE),
      `R%high` = sum(total_df == 3, na.rm = TRUE) / total_cells * 100,
      `R%medium` = sum(total_df == 2, na.rm = TRUE) / total_cells * 100,
      `R%low` = sum(total_df == 1, na.rm = TRUE) / total_cells * 100,
      `R%None` = sum(total_df == 0, na.rm = TRUE) / total_cells * 100
    )
  )
  output_df <- rbind.data.frame(all_stressors_df, output_df)
  return(output_df)
}


many_hra <- function(raster_list, dist_list, criteria, equation = c("euclidean", "multiplicative"), r_max = 3, n_overlap = NULL) {
  if (!is.list(raster_list)) {
    stop("'raster_list' must be a list of list named after each species/habitat")
  }
  if (!is.list(dist_list)) {
    stop("'dist_list' must' be a list of list named after each species/habitat")
  }
  if (all(!names(criteria) %in% c("SPECIES", "STRESSOR", "ATTRIBUTES", "RATING", "DQ", "WEIGHT", "E/C"))){
    stop("Criteria table must have the following columns: 'SPECIES', 'STRESSOR', 'ATTRIBUTES', 'RATING', 'DQ', 'WEIGHT', and 'E/C'")
  }
  equation <- match.arg(equation)
  stressors <- unique(criteria$STRESSOR)
  stressors <- stressors[!stressors %in% c(NA, "", "NA")]

  if (is.null(n_overlap)) {
    n_overlap <- length(stressors)
  }

  results <- list()

  for (species in names(raster_list)) {
    results[[species]] <- hra(raster_list[[species]], dist_list[[species]], criteria[criteria$SPECIES == species,], equation, r_max, n_overlap)
  }

}
