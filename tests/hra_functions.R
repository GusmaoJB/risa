# Function to perform habitat risk assessment (HRA)

# Creating test data
species <- data.frame(long = rnorm(80, 0, 10),
                   lat = rnorm(80, 0, 10), species = "sp1")

stressor <- data.frame(long = rnorm(100, 0, 5),
                   lat = rnorm(100, 0, 10), stressor = "trawling")

# Create kernel maps of species and stressor distributions and overlap maps
input_maps <- risa_prep(species, stressor, output_layer_type = "raster")

# Loading criteria table
reshape_habitat_data <- function(path) {
  # 1. read & drop truly empty lines
  lines <- readLines(path)
  lines <- lines[nzchar(lines)]

  # 2. auto‐detect delimiter on the header row
  hdr     <- lines[1]
  n_comma <- sum(strsplit(hdr, "")[[1]] == ",")
  n_tab   <- sum(strsplit(hdr, "")[[1]] == "\t")
  delim   <- if (n_tab > n_comma) "\t" else ","

  # 3. split each line
  parts <- strsplit(lines, delim, fixed = TRUE)

  # 4. grab species name from the very first row, 2nd field
  species <- parts[[1]][2]

  out <- list()
  current_stressor <- NA_character_
  in_attr   <- FALSE
  in_stress <- FALSE

  # 5. loop over the rest
  for (row in parts[-1]) {
    # —— NEW: skip any row where every field is blank
    if (all(trimws(row) == "")) next

    # detect start of Attributes
    if (row[1] == "HABITAT RESILIENCE ATTRIBUTES") {
      in_attr   <- TRUE
      in_stress <- FALSE
      next
    }
    # detect the overall Stressor header
    if (row[1] == "HABITAT STRESSOR OVERLAP PROPERTIES") {
      in_attr   <- FALSE
      in_stress <- FALSE
      next
    }

    # detect the start of any stressor block:
    #  "row2 == 'RATING'" but not one of our two section headers
    if (length(row) >= 2 &&
        row[2] == "RATING" &&
        ! row[1] %in% c(
          "HABITAT RESILIENCE ATTRIBUTES",
          "HABITAT STRESSOR OVERLAP PROPERTIES"
        )
    ) {
      current_stressor <- row[1]
      in_attr   <- FALSE
      in_stress <- TRUE
      next
    }


    # skip any header‐rows (those that say “RATING” in col 2)
    if (length(row) >= 2 && row[2] == "RATING") next

    # pad short rows so we can safely index 1:5
    if (length(row) < 5) {
      row <- c(row, rep(NA_character_, 5 - length(row)))
    }

    # build a data.frame row
    out[[length(out) + 1]] <- data.frame(
      HABITAT_NAME               = species,
      ATRIBUTES_AND_PROPERTIES   = row[1],
      STRESSOR                   = if (in_stress) current_stressor else NA_character_,
      RATING                     = as.numeric(row[2]),
      DQ                         = as.numeric(row[3]),
      WEIGHT                     = as.numeric(row[4]),
      `E/C`                      = row[5],
      stringsAsFactors           = FALSE,
      check.names                = FALSE
    )
  }

  # bind & return
  do.call(rbind, out)
}



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
