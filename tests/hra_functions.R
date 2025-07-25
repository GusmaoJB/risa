# Function to perform habitat risk assessment (HRA)

# Creating test data
species1 <- data.frame(long = rnorm(80, 0, 10),
                   lat = rnorm(80, 0, 10), species = "species1")

stressor <- data.frame(long = rnorm(100, 0, 5),
                   lat = rnorm(100, 0, 10), stressor = "stressor1")

# Create kernel maps of species and stressor distributions and overlap maps
input_maps <- risa_prep(species, stressor, output_layer_type = "raster")

library(risa)




df_long <- reshape_habitat_data("C:/Users/gusma/Documents/risa_maps_test/criteria/test_criteria.csv")
df_long
dim(df_long)

lines <- readLines("C:/Users/gusma/Documents/risa_maps_test/criteria/test_criteria.csv")
test <- read.csv("C:/Users/gusma/Documents/risa_maps_test/criteria/test_criteria.csv")
test
lines <- lines[nzchar(lines)]
lines


tes <- split_by_habitat("C:/Users/gusma/Documents/risa_maps_test/criteria/multi_species_criteria.csv")
tes

split_by_habitat <- function(path) {
  # 1) read & drop truly blank lines
  lines <- readLines(path)
  lines <- lines[nzchar(lines)]

  # 2) sniff delimiter on the very first line
  delim <- if (any(grepl("\t", lines[1], fixed = TRUE))) "\t" else ","

  # 3) parse out the two header rows
  h1 <- strsplit(lines[1], delim, fixed = TRUE)[[1]]
  h2 <- strsplit(lines[2], delim, fixed = TRUE)[[1]]

  # 4) locate the metric‐blocks and the E/C column
  rating_cols <- which(h2 == "RATING")
  dq_cols     <- which(h2 == "DQ")
  weight_cols <- which(h2 == "WEIGHT")
  ec_col      <- which(h2 == "E/C")[1]
  habitats    <- h1[rating_cols]

  # 5) bulk‐read the rest of the table (skip the two header rows)
  df <- read.table(
    path,
    header         = FALSE,
    sep            = delim,
    skip           = 2,
    stringsAsFactors = FALSE,
    fill           = TRUE,
    comment.char   = ""
  )
  # pad to the same width as the header row
  if (ncol(df) < length(h2)) {
    df[ , (ncol(df)+1):length(h2)] <- NA
  }

  # 6) drop
  #    • any row that’s entirely blank
  #    • the “HABITAT STRESSOR OVERLAP PROPERTIES” marker
  #    • any row whose second column is literally the word "RATING"
  df <- df[!apply(df, 1, function(r) all(is.na(r) | r == "")), ]
  df <- df[df[[1]] != "HABITAT STRESSOR OVERLAP PROPERTIES", ]
  df <- df[ df[[2]] != "RATING", ]

  # 7) slice out one data.frame per habitat
  out <- setNames(vector("list", length(habitats)), habitats)
  for (i in seq_along(habitats)) {
    out[[i]] <- data.frame(
      # first column: attribute/property name
      ATTRIBUTES_AND_PROPERTIES = df[[1]],
      # that habitat’s RATING/DQ/WEIGHT
      RATING =   suppressWarnings(as.numeric(df[[ rating_cols[i] ] ])),
      DQ     =   suppressWarnings(as.numeric(df[[ dq_cols[i]     ] ])),
      WEIGHT =  suppressWarnings(as.numeric(df[[ weight_cols[i] ] ])),
      # the shared E/C column
      `E/C`  =   df[[ ec_col ]],
      stringsAsFactors = FALSE,
      check.names      = FALSE
    )
  }

  out
}








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
