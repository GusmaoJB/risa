habitat_risk_assessment <- function(raster_list, criteria, equation = c("euclidean", "multiplicative"), n_classes = 3) {

  if (!inherits(raster_list, "risa")) {
    stop("Input 'raster_list' must be a 'risa' object.")
  }

  equation <- match.arg(equation)
  stressors <- unique(criteria$STRESSOR)
  stressors <- stressors[!stressors %in% c(NA, "", "NA")]

  if (!all(stressors %in% names(raster_list$stressor_kernel_maps))) {
    stop("The stressor names in 'criteria' must match the names in 'raster_list$stressor_kernel_maps'.")
  }

  risk_results <- list()

  sp_distr <- raster_list$species_distributions[[1]]$raster
  sp_presence <- terra::ifel(!is.na(sp_distr), 1, NA)

  for (stressor in stressors) {
    crit_stressor <- criteria[criteria$STRESSOR == stressor,]

    # Split consequence (C) and exposure (E)
    C_df <- crit_stressor[crit_stressor$`E/C` == "C",]
    E_df <- crit_stressor[crit_stressor$`E/C` == "E",]

    C_df$DQ <- as.integer(C_df$DQ)
    C_df$WEIGHT <- as.integer(C_df$WEIGHT)
    E_df$DQ <- as.integer(E_df$DQ)
    E_df$WEIGHT <- as.integer(E_df$WEIGHT)

    # Define mapped rasters
    raster_mapping_C <- raster_list$stressor_kernel_maps[[stressor]]$raster
    raster_mapping_E <- raster_list$overlap_maps[[1]][[stressor]]$raster

    # Normalize and aggregate C and E scores
    C_score_raster <- raster_mapping_C / terra::global(raster_mapping_C, "max", na.rm=TRUE)[[1]]
    E_score_raster <- raster_mapping_E / terra::global(raster_mapping_E, "max", na.rm=TRUE)[[1]]

    # Apply weighting
    C_weighted <- sum(C_df$WEIGHT * C_df$DQ, na.rm=TRUE)
    E_weighted <- sum(E_df$WEIGHT * E_df$DQ, na.rm=TRUE)

    C_score_raster <- (C_score_raster * C_weighted) / sum(C_df$WEIGHT, na.rm=TRUE)
    E_score_raster <- (E_score_raster * E_weighted) / sum(E_df$WEIGHT, na.rm=TRUE)

    # Converting species distributions into a raster with zeros
    sp_distr_zeros <- terra::ifel(!is.na(sp_distr), 0, NA)

    # Setting nozes with zero E and C criteria
    C_score_raster_map <- terra::mosaic(C_score_raster, sp_distr_zeros, fun="first")
    E_score_raster_map <- terra::mosaic(E_score_raster, sp_distr_zeros, fun="first")

    # Compute risk scores
    r_max <- n_classes
    if (equation == "multiplicative") {
      m_jkl <- r_max^2
      risk_scores <- C_score_raster * E_score_raster * sp_presence
      risk_scores <- terra::mosaic(risk_scores, sp_distr_zeros, fun="first")
    } else if (equation == "euclidean") {
      m_jkl <- sqrt(2 * (r_max - 1)^2)
      risk_scores <- sqrt((E_score_raster - 1)^2 + (C_score_raster - 1)^2) * sp_presence
      risk_scores <- terra::mosaic(risk_scores, sp_distr_zeros, fun="first")
    }

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

  # Combine stressor risk rasters
  total_risk <- Reduce("+", lapply(risk_results, function(x) x$Risk_map_raw))

  risk_results$total <- total_risk

  return(risk_results)
}





# Creating test data
spp_df <- rbind(data.frame(long = rnorm(80, 0, 10),
                           lat = rnorm(80, 0, 10), species = "species1"))
str_df <- rbind(data.frame(long = rnorm(100, 0, 5),
                           lat = rnorm(100, 0, 10), stressor = "stressor1"),
                data.frame(long = rnorm(50, 0, 10),
                           lat = rnorm(100, 0, 5), stressor = "stressor2"))

# Create kernel maps of species and stressor distributions and overlap maps
risa_maps <- risa_prep(spp_df, str_df)

#Load example data
path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
df <- read.csv(path)

#Inspect dataframe
df[,1:8]

#Reshape criteria table
crit_list <- criteria_reshape(df)

test <- habitat_risk_assessment(risa_maps, crit_list[[1]], equation = "euclidean", n_classes = 3)


terra::plot(test$stressor1$E_criteria)
terra::plot(test$stressor1$C_criteria)
terra::plot(test$stressor2$E_criteria)
terra::plot(test$stressor1$Risk_map_raw)
terra::plot(test$stressor1$Risk_map)
terra::plot(test$stressor2$Risk_map_raw)
terra::plot(test$stressor2$Risk_map)

terra::plot(test$total)

crit_list
