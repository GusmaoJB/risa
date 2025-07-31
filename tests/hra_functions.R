# Function to perform habitat risk assessment (HRA)

# Creating test data
species1 <- data.frame(long = rnorm(80, 0, 2),
                       lat = rnorm(80, 0, 2), species = "species1")

stressor1 <- data.frame(long = rnorm(100, 0, 1),
                        lat = rnorm(100, 0, 2), stressor = "stressor1")

# Create kernel maps of species and stressor distributions and overlap maps
raster_list <- risa_prep(species1, stressor1, output_layer_type = "both")

# Loand example data and reshape for analysis
path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
criteria_csv <- read.csv(path)

#Reshape criteria table
crit_list <- criteria_reshape(criteria_csv)
criteria <- crit_list[[1]]
criteria <- criteria[!criteria$STRESSOR == "stressor2",]
criteria <- cbind.data.frame(SPECIES = c(rep("species1", dim(criteria)[1])), criteria)

criteria[,c("RATING", "DQ", "WEIGHT")] <- ifelse(criteria[,c("RATING", "DQ", "WEIGHT")] == 2, 3, NA)


# Split consequence (C) and exposure (E)
C_df <- criteria[criteria$`E/C` == "C",]
E_df <- criteria[criteria$`E/C` == "E",]

C_df$DQ <- as.integer(C_df$DQ)
C_df$WEIGHT <- as.integer(C_df$WEIGHT)

E_df$DQ <- as.integer(E_df$DQ)
E_df$WEIGHT <- as.integer(E_df$WEIGHT)

# Define which rasters correspond to your mapped‐criteria names
raster_mapping_C <- list(
  intensity = raster_list$stressor_kernel_maps$stressor1$raster
)
raster_mapping_E <- list(
  likelihood_of_interaction = raster_list$overlap_maps$species1$stressor1$raster
)

# Select constant criteria (with rating), and divide value by DQ * Weight
E_const <- E_df[!is.na(E_df$RATING),]
C_const <- C_df[!is.na(C_df$RATING),]

E_numer_const <- sum(E_const$RATING / (E_const$DQ * E_const$WEIGHT))
C_numer_const <- sum(C_const$RATING / (C_const$DQ * C_const$WEIGHT))

# Select mapped criteria (RATING == NA), then sum over raster/(DQ*WEIGHT)
E_mapped <- E_df[is.na(E_df$RATING),]
C_mapped <- C_df[is.na(C_df$RATING),]

C_numer_rast <- 0
for(i in seq_len(nrow(C_mapped))){
  crit <- E_mapped[i,]
  attr_name <- names(crit)[3]
  r <- raster_mapping_C[[1]]
  C_numer_rast <- C_numer_rast + (r / (crit$DQ * crit$WEIGHT))
}

E_numer_rast <- 0
for(i in seq_len(nrow(E_mapped))){
  crit <- E_mapped[i,]
  attr_name <- names(crit)[3]
  r <- raster_mapping_E[[1]]
  E_numer_rast <- E_numer_rast + (r / (crit$DQ * crit$WEIGHT))
}

terra::plot(C_numer_rast)
terra::plot(E_numer_rast)

sp1_distr <- raster_list$species_distributions$species1$raster
sp1_distr <- terra::ifel(!is.na(sp1_distr),0,NA)

#C_numer_rast <- terra::mosaic(sp1_distr, C_numer_rast, fun="last")
#E_numer_rast <- terra::mosaic(sp1_distr, E_numer_rast, fun="last")

# Define denominator (is the same for all cells)
E_denom <- sum(1/(E_df$DQ * E_df$WEIGHT))
C_denom <- sum(1/(C_df$DQ * C_df$WEIGHT))

# Final score raster
E_score_raster <- (E_numer_const + E_numer_rast) / E_denom
C_score_raster <- (C_numer_const + C_numer_rast) / C_denom

# Mask to where the species is actually present (>0) and compute risk
presence_mask <- raster_list$species_kernel_maps$species1$raster > 0

# For multiplicative risk estimates = E*C
risk_HRA_multi <- terra::mosaic((C_score_raster * E_score_raster), sp1_distr, fun="first")

terra::plot(risk_HRA_multi, main="Bycatch risk (Multiplicative)")

# For Euclidean risk estimattes
eucl_scores <- (sqrt((E_score_raster  - 1)^2 + (C_score_raster  - 1)^2))
risk_HRA_eucl <- terra::mosaic(eucl_scores, sp1_distr, fun="first")
terra::plot(risk_HRA_eucl, main="Bycatch risk (Euclidean)")

# Now we need to reclassify our Euclidean raster to 1-3 scores
re_mat_HRA_eucl <- reclass_matrix(eucl_scores, n_classes = 3, exclude_lowest = FALSE, custom_max=2.828427)
risk_HRA_eucl_reclass <- terra::classify(risk_HRA_eucl, re_mat_HRA_eucl, include.lowest=TRUE)

# plot to check
terra::plot(risk_HRA_eucl_reclass, main="Bycatch risk (1–3)")

# Generate summary statistics
library(dplyr)
summarize_hra <- function(E_raster,
                          C_raster,
                          R_raster,
                          R_reclass_raster) {

  s <- c(E_raster, C_raster, R_raster, R_reclass_raster)
  names(s) <- c("E", "C", "R", "R_reclass")
  df <- as.data.frame(s, na.rm = FALSE)

  summary_df <- df %>%
    summarise(
      E_min   = min(E,   na.rm = TRUE),
      E_max   = max(E,   na.rm = TRUE),
      E_mean  = mean(E,  na.rm = TRUE),
      C_min   = min(C,   na.rm = TRUE),
      C_max   = max(C,   na.rm = TRUE),
      C_mean  = mean(C,  na.rm = TRUE),
      R_min   = min(R,   na.rm = TRUE),
      R_max   = max(R,   na.rm = TRUE),
      R_mean  = mean(R,  na.rm = TRUE),
      total_cells = n(),
      `R%high` = sum(R_reclass == 3, na.rm = TRUE) / total_cells * 100,
      `R%medium` = sum(R_reclass == 2, na.rm = TRUE) / total_cells * 100,
      `R%low` = sum(R_reclass == 1, na.rm = TRUE) / total_cells * 100,
      `R%None` = sum(R_reclass == 0) / total_cells * 100
    )
  return(as.data.frame(summary_df, check.names = FALSE))
}

# Include areas where risk is zero (species occurr but no risk) in E and C score rasters
E_score_raster_z <- terra::mosaic(E_score_raster, sp1_distr, fun="first")
C_score_raster_z <- terra::mosaic(C_score_raster, sp1_distr, fun="first")

summarize_hra(E_score_raster_z, C_score_raster_z, risk_HRA_eucl, risk_HRA_eucl_reclass)

export_maps(raster_list, "C:/Users/gusma/Documents/research/test_hra/maps")


