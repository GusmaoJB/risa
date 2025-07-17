# Function to perform habitat risk assessment (HRA)

# Loading data and dummy criteria table
pont_shp <- st_read(paste(path_spps, "pont.shp", sep="/"))
traw_shp <- st_read(paste(path_spps, "traw.shp", sep="/"))
criteria <- read.csv(paste(dir, "byra_docs/criteria_scores_test.csv", sep="/"))
criteria

input_maps <- risa_prep(pont_shp, traw_shp, output_layer_type = "raster")

# Split consequence (C) and exposure (E)
C_df <- criteria %>% filter(E_C == "C")
E_df <- criteria %>% filter(E_C == "E")

# Function to compute a score raster for one block (C or E)
compute_block <- function(df_block, raster_mapping){

  # Select constant criteria (with rating), and divide value by DQ * Weight
  const <- df_block %>% filter(!is.na(RATING))
  numer_const <- sum(const$RATING / (const$DQ * const$WEIGHT))

  # Select mapped criteria (RATING == NA), then sum over raster/(DQ*WEIGHT)
  mapped <- df_block %>% filter(is.na(RATING))
  numer_rast <- 0
  for(i in seq_len(nrow(mapped))){
    crit <- mapped[i,]
    attr_name <- crit$HABITAT_STRESSOR_ATTRIBUTES
    r <- raster_mapping[[attr_name]]
    numer_rast <- numer_rast + (r / (crit$DQ * crit$WEIGHT))
  }

  # Define denominator (is the same for all cells)
  denom <- sum(1/(df_block$DQ * df_block$WEIGHT))

  # Assemble final score raster
  score_raster <- (numer_const + numer_rast) / denom
  return(score_raster)
}

# Define which rasters correspond to your mapped‐criteria names ---
raster_mapping_C <- list(
  intensity = raster_list$traw
)
raster_mapping_E <- list(
  likelihood_of_interaction = raster_list$overlap
)

# Compute overall C and E rasters
C_score <- compute_block(C_df, raster_mapping_C)
E_score <- compute_block(E_df, raster_mapping_E)

# Mask to where the dolphin is actually present (>0) and compute risk
presence_mask <- raster_list$pont > 0

# For multiplicative risk estimates = E*C
risk_HRA_multi <- mask(C_score * E_score, presence_mask)
plot(risk_HRA_multi, main="Bycatch risk (Multiplicative)")

# For Euclidean risk estimattes
risk_HRA_eucl <- mask(sqrt((E_score - 1)^2 + (C_score - 1)^2), presence_mask)
plot(risk_HRA_eucl, main="Bycatch risk (Euclidean)")

# Now we need to reclassify our raster to 1-3 scores
re_mat_HRA <- reclass_matrix(risk_HRA_eucl, n_classes = 3, exclude_lowest = FALSE)
risk_HRA_reclassified <- classify(risk_HRA_eucl, re_mat_HRA, include.lowest=TRUE)

# plot to check
plot(risk_HRA_reclassified, main="Bycatch risk (1–3)")

# Generate summary statistics
summarize_hra <- function(E_raster,
                          C_raster,
                          R_raster) {

  s <- c(E_raster, C_raster, R_raster)
  names(s) <- c("E", "C", "R")
  df <- as.data.frame(s, na.rm = FALSE)
  df[is.na(df)] <- 0

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
      `R%low` = sum(R == 1, na.rm = TRUE) / total_cells * 100,
      `R%medium` = sum(R == 2, na.rm = TRUE) / total_cells * 100,
      `R%high` = sum(R == 3, na.rm = TRUE) / total_cells * 100,
      `R%None` = sum(R == 0) / total_cells * 100
    )
  return(as.data.frame(summary_df, check.names = FALSE))
}

summarize_hra(E_score, C_score, risk_HRA_reclassified)
