#' Round half up
#'
#' Rounds numbers to the specified number of digits, rounding halves up
#' rather than using base R's default "round half to even" behavior.
#'
#' @param x numeric vector.
#' @param digits integer, number of decimal digits to round to (default 0).
#'
#' @return A numeric vector of the same length as \code{x}, rounded.
#' @keywords internal
round_half_up <- function(x, digits = 0) {
  posneg <- sign(x)
  z <- abs(x) * 10^digits
  z <- z + 0.5 + sqrt(.Machine$double.eps)
  z <- trunc(z)
  z <- z / 10^digits
  z * posneg
}

#' Classify a numeric quality score into a stoplight color
#'
#' Rounds each value to the nearest whole number (0-3) and maps it to a
#' colorblind-friendly stoplight color: dark red (1), saturated yellow (2),
#' or teal green (3). A rounded value of 0, or \code{NA}, maps to black.
#'
#' @param x numeric vector with values ranging from 0 to 3.
#'
#' @return A character vector of hex color codes, the same length as \code{x}.
#' @keywords internal
stoplight_cols <- function(x) {
  if (!is.numeric(x)) {
    stop("Error: input must be a number ranging from 0 to 3")
  }
  if (any(x > 3, na.rm = TRUE) || any(x < 0, na.rm = TRUE)) {
    stop("Error: input must be a number ranging from 0 to 3")
  }
  r_x <- round_half_up(x)
  ifelse(r_x == 0 | is.na(r_x), "#000000",
         ifelse(r_x == 1, "#8B0000",
                ifelse(r_x == 2, "#FFFF00", "#009E73")))
}

#' Convert a data-quality score table into stoplight color tables
#'
#' Takes a quality/uncertainty scoring data.frame and converts the numeric
#' scores into colorblind-friendly stoplight colors at three levels:
#' ecosystem (per data type, pooled across everything), species (per
#' species x stressor x data type), and species with stressors combined
#' (per species x data type, averaged across stressors).
#'
#' @param quality.data.frame A data.frame with (in this column order)
#'   \code{level}, \code{data}, \code{species}, \code{stressor}, and
#'   \code{uncertainty}. \code{species} and \code{stressor} may be \code{NA}
#'   for rows not specific to a species or stressor (e.g. species-level or
#'   stressor-level rows). \code{uncertainty} must be numeric, from 0 to 3.
#'
#' @return A named list of three data.frames:
#' \describe{
#'   \item{ecosystem}{\code{light} (data type) x \code{value} (color),
#'     averaged across all rows for that data type.}
#'   \item{species_per_stressor}{\code{species} x \code{stressor} x
#'     \code{light} (data type) x \code{value} (color).}
#'   \item{species_comb_stressor}{\code{species} x \code{light} (data type)
#'     x \code{value} (color), averaged across stressors.}
#' }
#'
#' @examples
#' quality_data <- cbind.data.frame(
#'   level = c("species", "species", "species", "species",
#'             "stressor", "stressor",
#'             "interaction", "interaction", "interaction", "interaction"),
#'   data = c("Animal sightings distribution", "Habitat suitability",
#'            "Animal sightings distribution", "Habitat suitability",
#'            "Fishing effort", "Fishing effort",
#'            "Bycatch/stranding data", "Bycatch/stranding data",
#'            "Bycatch/stranding data", "Bycatch/stranding data"),
#'   species = c("dolphin", "dolphin", "turtle", "turtle", NA, NA,
#'               "dolphin", "turtle", "dolphin", "turtle"),
#'   stressor = c(NA, NA, NA, NA, "trawling", "gillnet",
#'                "trawling", "trawling", "gillnet", "gillnet"),
#'   uncertainty = c(3, 2, 2, 2, 2, 3, 1, 2, 1, 2)
#' )
#' scores_to_colors(quality_data)
#'
#' @export
scores_to_colors <- function(quality.data.frame) {
  
  quality_data <- quality.data.frame
  names(quality_data) <- c("level", "data", "species", "stressor", "uncertainty")
  
  spp <- na.omit(unique(quality_data$species))
  str <- na.omit(unique(quality_data$stressor))
  
  # Ecosystem stoplight
  ecosys_risk_df <- stats::aggregate(uncertainty ~ data, data = quality_data, FUN = mean)
  ecosys_risk_cols <- ecosys_risk_df
  names(ecosys_risk_cols) <- c("light", "value")
  ecosys_risk_cols$value <- stoplight_cols(ecosys_risk_df$uncertainty)
  
  # Species stoplight
  # Species per stressor
  spp_risk_df <- data.frame()
  for (sp in spp) {
    sp_subset <- quality_data[quality_data$species == sp | is.na(quality_data$species), ]
    sp_subset$species <- sp
    spp_risk_df <- rbind.data.frame(spp_risk_df, sp_subset)
  }
  spp_risk_df_2 <- na.omit(spp_risk_df)
  
  for (sp in spp) {
    sp_subset <- quality_data[!is.na(quality_data$species), ]
    sp_subset <- sp_subset[sp_subset$species == sp, ]
    for (st in str) {
      st_subset <- sp_subset[is.na(sp_subset$stressor), ]
      st_subset$stressor <- st
      spp_risk_df_2 <- rbind.data.frame(spp_risk_df_2, st_subset)
    }
  }
  
  spp_risk_vals <- spp_risk_df_2[, c("species", "stressor", "data", "uncertainty")]
  names(spp_risk_vals) <- c("species", "stressor", "light", "value")
  spp_risk_cols <- spp_risk_vals
  spp_risk_cols$value <- stoplight_cols(spp_risk_vals$value)
  
  # Species with combined stressors
  spp_comb_risk_vals <- stats::aggregate(value ~ species + light, data = spp_risk_vals, FUN = mean)
  spp_comb_risk_cols <- spp_comb_risk_vals
  spp_comb_risk_cols$value <- stoplight_cols(spp_comb_risk_vals$value)
  
  list(
    ecosystem = ecosys_risk_cols,
    species_per_stressor = spp_risk_cols,
    species_comb_stressor = spp_comb_risk_cols
  )
}