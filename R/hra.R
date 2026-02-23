#' Habitat Risk Assessment (HRA)
#'
#' Implements the Natural Capital Project’s InVEST Habitat Risk Assessment (HRA;
#' Sharp et al., 2016) model to estimate spatial risk from multiple stressors to
#' a species/habitat (single-species mode) or across multiple species (ecosystem
#' mode). For each stressor, the function combines Exposure (E) and Consequence (C)
#' criteria—mixing mapped rasters and constant scores—using either the Euclidean
#' or Multiplicative InVEST risk equations (Arkema et al., 2014). Inputs can be
#' nested lists of rasters (stressor → attributes; or species → stressor → attributes),
#' a species distribution raster (or list), and a criteria table (or list). Rasters
#' are auto-aligned to the species grid; outputs include per-stressor E/C maps,
#' raw and classed risk maps, ecosystem risk (when applicable), and a `summary_stats`
#' data frame. Optional reprojection to EPSG:4326 is supported.
#'
#' @section Decay functions:
#' A distance-based decay can be applied to E and C maps (and a general stressor
#' occurrence decay) when a named buffer (in meters) is supplied for each stressor.
#' Let \eqn{d} be the distance (in meters) from a cell to the nearest mapped
#' stressor presence, \eqn{B} the buffer distance for that stressor, and
#' \eqn{I(d \le B)} an indicator for being within the buffer. The implemented
#' options are:
#'
#' \itemize{
#'   \item \code{"none"} (InVEST): \deqn{f(d) = I(d \le B)}
#'   \item \code{"linear"} (InVEST): \deqn{f(d) = I(d \le B)\,\left(1 - \frac{d}{B}\right)}
#'   \item \code{"exponential"} (InVEST): \deqn{f(d) = I(d \le B)\,e^{-k d},\quad k = -\frac{\ln(10^{-6})}{B}}
#'   \item \code{"polynomial_2nd"} (extension): \deqn{f(d) = I(d \le B)\,\left(1 - \frac{d}{B}\right)^{2}}
#'   \item \code{"polynomial_3rd"} (extension): \deqn{f(d) = I(d \le B)\,\left(1 - \frac{d}{B}\right)^{3}}
#'   \item \code{"complementary_decay_2nd"} (extension): \deqn{f(d) = I(d \le B)\,\left(1 - \left(\frac{d}{B}\right)^{2}\right)}
#'   \item \code{"complementary_decay_3rd"} (extension): \deqn{f(d) = I(d \le B)\,\left(1 - \left(\frac{d}{B}\right)^{3}\right)}
#' }
#'
#' \strong{Note:} Only the first three options (\code{"none"}, \code{"linear"},
#' and \code{"exponential"}) are present in InVEST’s reference implementation.
#' The polynomial and complementary variants are smooth extensions that allow
#' faster-than-linear attenuation while remaining bounded in \eqn{[0,1]} on \eqn{[0,B]}.
#'
#' @references
#' Arkema KK, Verutes G, Bernhardt JR, Clarke C, Rosado S, Canto M, et al. (2014)
#' Assessing habitat risk from human activities to inform coastal and marine spatial
#' planning: a demonstration in Belize. \emph{Environmental Research Letters} 9:114016.
#'
#' Sharp R, Tallis H, Ricketts T, Guerry A, Wood SA, Chaplin-Kramer R, et al. (2016)
#' InVEST [Internet]. The Natural Capital Project, Stanford University, University of
#' Minnesota, The Nature Conservancy, and World Wildlife Fund. Available at
#' \url{http://www.naturalcapitalproject.org/software/}.
#'
#' @param raster_list A named list of rasters. In single-species mode (depth 2),
#'   it must be \code{list(stressor -> list(attribute -> SpatRaster))}. In ecosystem
#'   mode (depth 3), it must be \code{list(species -> list(stressor -> list(attribute -> SpatRaster)))}.
#' @param species_distr A \code{SpatRaster} (single-species) OR a named list of
#'   \code{SpatRaster} objects (ecosystem). If an element is a list, the first
#'   \code{SpatRaster} within it will be used.
#' @param criteria A \code{data.frame} (single) OR named list of \code{data.frame}s
#'   (ecosystem) with columns \code{STRESSOR}, \code{ATTRIBUTES}, \code{RATING},
#'   \code{DQ}, \code{WEIGHT}, and \code{E/C}.
#' @param equation \code{c("euclidean","multiplicative")}. Selects the InVEST risk
#'   equation used to combine E and C.
#' @param r_max Integer in \code{1..10}; maximum score used by the chosen equation.
#' @param n_overlap Optional integer: number of stressors assumed to overlap for
#'   classification scaling. By default it is inferred from the criteria.
#' @param output_decimal_crs Logical; if \code{TRUE}, reprojects outputs to EPSG:4326.
#' @param decay \code{c("none","linear","exponential","polynomial_2nd","polynomial_3rd",
#'   "complementary_decay_2nd","complementary_decay_3rd")}. See \emph{Decay functions}
#'   above. Only the first three are provided by InVEST; the others are extensions.
#' @param buffer_m Named numeric vector giving the buffer distance (meters) for each
#'   stressor (e.g., \code{c(stressor1 = 500000, stressor2 = 1000000)}). Required
#'   when \code{decay != "none"} for any stressor you wish to decay.
#'
#' @return An object of class \code{"risaHRA"}: a named list that, for each stressor,
#'   contains \code{E_criteria}, \code{C_criteria}, \code{Risk_map_raw}, and
#'   \code{Risk_map}; plus \code{total_raw}, \code{total}, and \code{summary_stats}
#'   in single-species mode. In ecosystem mode, per-species results are returned,
#'   along with \code{ecosys_risk_raw}, \code{ecosys_risk_classified}, and a
#'   combined \code{summary_stats} table.
#'
#' @details
#' Rasters are aligned to the species grid using nearest-neighbor for categorical
#' data and bilinear interpolation for continuous data. When a decay is used and
#' a buffer is supplied for a given stressor, E and C rasters are attenuated by
#' the chosen \eqn{f(d)} and an additional "general" decay is applied to overall
#' stressor occurrence. Cells outside species presence are masked.
#'
#' @importFrom terra compareGeom project mosaic mask ifel freq global app rast sprc cover
#'
#' @examples
#' # Creating test data
#' set.seed(12)
#' spp_df <- rbind(data.frame(long = rnorm(80, 0, 10),
#' lat = rnorm(80, 0, 10), species = "species1"),
#' data.frame(long = rnorm(60, 0, 10),
#' lat = rnorm(60, 0, 10), species = "species2"))
#' str_df <- rbind(data.frame(long = rnorm(100, 0, 5),
#' lat = rnorm(100, 0, 10), stressor = "stressor1"),
#' data.frame(long = rnorm(50, 0, 10),
#' lat = rnorm(100, 0, 5), stressor = "stressor2"))
#' # Create kernel maps of species and stressor distributions and overlap maps
#' risa_maps <- risa_prep(spp_df, str_df)
#' #Load example data
#' path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
#' df <- read.csv(path)
#' #Reshape criteria table
#' crit_list <- criteria_reshape(df)
#' # Selecting spatially explicit criteria ratings
#' # Note that the rasters in the stressors's list are named after the respective attribute in the criteria table.
#' rast_list <- list(
#' species1 = list(
#' stressor1 = list(
#' intensity = risa_maps$stressor_kernel_maps$stressor1$raster,
#' `likelihood of interaction`=risa_maps$overlap_maps$species1$stressor1$raster),
#' stressor2 = list(
#' intensity = risa_maps$stressor_kernel_maps$stressor2$raster,
#' `likelihood of interaction`=risa_maps$overlap_maps$species1$stressor2$raster)
#' ),
#' species2 = list(
#' stressor1 = list(
#' intensity = risa_maps$stressor_kernel_maps$stressor1$raster,
#' `likelihood of interaction`=risa_maps$overlap_maps$species2$stressor1$raster),
#' stressor2 = list(
#' intensity = risa_maps$stressor_kernel_maps$stressor2$raster,
#' `likelihood of interaction`=risa_maps$overlap_maps$species2$stressor2$raster)
#' )
#' )
#'
#' # Species' distributions rasters
#' spp_dist <- list(species1 = risa_maps$species_distributions$species1$raster,
#' species2 = risa_maps$species_distributions$species2$raster)
#'
#' # Simple example with one species and one stressor
#' test1 <- hra(rast_list[[1]], spp_dist[[1]], crit_list[[1]], equation = "euclidean")
#' terra::plot(test1$total_raw)
#' test1$summary_stats
#'
#' # Now with two species and two stressors
#' many_test <- hra(rast_list, spp_dist, crit_list, equation = "euclidean")
#' many_test$summary_stats
#' terra::plot(many_test$species1$total_raw)
#' terra::plot(many_test$species2$total_raw)
#' @export
hra <- function(
    raster_list, species_distr, criteria,
    equation = c("euclidean","multiplicative"),
    r_max = 3, n_overlap = NULL, output_decimal_crs = FALSE,
    decay = c("none", "linear", "exponential",
              "polynomial_2nd", "polynomial_3rd",
              "complementary_decay_2nd", "complementary_decay_3rd"),
    buffer_m = NULL) {

  depth <- list_depth_base(raster_list)
  equation <- match.arg(equation)
  decay <- match.arg(decay)
  m_jkl <- if (equation=="multiplicative") r_max^2 else sqrt(2*(r_max-1)^2)

  # Helpers
  # Check criteria input format
  .check_criteria <- function(df_or_list) {
    req <- c("STRESSOR","ATTRIBUTES","RATING","DQ","WEIGHT","E/C")
    if (is.data.frame(df_or_list)) {
      if (!all(req %in% names(df_or_list))) stop("Criteria must contain: ", paste(req, collapse=", "))
      return(df_or_list)
    }
    if (is.list(df_or_list)) {
      lapply(df_or_list, function(d) {
        if (!is.data.frame(d) || !all(req %in% names(d))) {
          stop("Each criteria element must contain: ", paste(req, collapse=", "))
        }
        d
      })
    } else stop("'criteria' must be a data.frame (single) or list (ecosystem).")
  }

  # Find/extract a SpatRaster if user passed a list container
  .as_raster <- function(x) {
    if (inherits(x, "SpatRaster")) return(x)
    if (is.list(x)) {
      for (el in x) {
        r <- .as_raster(el)
        if (inherits(r, "SpatRaster")) return(r)
      }
    }
    NULL
  }

  # Calculates C and E scores, as well as risk estimates, for a single species
  .single <- function(rlist, sp_distr, crit, equation,
                      r_max, n_overlap, output_decimal_crs,
                      decay, buffer_m) {
    if (list_depth_base(rlist) != 2L) stop("Single-species mode requires depth-2 raster_list.")

    # Accept list containers for sp_distr
    sp_distr <- .as_raster(sp_distr)
    if (!inherits(sp_distr, "SpatRaster")) stop("'species_distr' must be or contain a SpatRaster.")

    crit <- .check_criteria(crit)

    stressors <- unique(crit$STRESSOR)
    stressors <- stressors[!(is.na(stressors) | stressors %in% c("", "NA"))]
    if (length(stressors) == 0L) stop("No stressors found in criteria (STRESSOR column).")
    if (is.null(n_overlap)) n_overlap <- length(stressors)

    if (is.null(names(rlist)) || any(!nzchar(names(rlist)))) {
      stop("Top-level 'raster_list' must be a named list of stressors.")
    }
    if (!all(names(rlist) %in% stressors)) {
      stop("Names in 'raster_list' must be a subset of criteria STRESSOR values.")
    }

    sp_presence <- terra::ifel(!is.na(sp_distr), 1, NA)
    sp_distr_zeros <- terra::ifel(!is.na(sp_distr), 0, NA)
    zero_r <- sp_distr * 0

    res <- setNames(vector("list", length(names(rlist))), names(rlist))

    for (stressor in names(rlist)) {
      crit_stressor <- crit[is.na(crit$STRESSOR) | crit$STRESSOR %in% c("", "NA", stressor), , drop=FALSE]
      C_df <- crit_stressor[crit_stressor$`E/C` == "C", , drop=FALSE]
      E_df <- crit_stressor[crit_stressor$`E/C` == "E", , drop=FALSE]
      if (nrow(E_df) == 0L || nrow(C_df) == 0L) stop("Both E and C must have at least one criterion for '", stressor, "'.")

      suppressWarnings({
        E_df$DQ <- as.numeric(E_df$DQ); E_df$WEIGHT <- as.numeric(E_df$WEIGHT); E_df$RATING <- as.numeric(E_df$RATING)
        C_df$DQ <- as.numeric(C_df$DQ); C_df$WEIGHT <- as.numeric(C_df$WEIGHT); C_df$RATING <- as.numeric(C_df$RATING)
      })
      if (any(E_df$DQ<=0 | E_df$WEIGHT<=0, na.rm=TRUE) || any(C_df$DQ<=0 | C_df$WEIGHT<=0, na.rm=TRUE)) {
        stop("DQ/WEIGHT must be > 0 for stressor '", stressor, "'.")
      }

      E_const <- E_df[!is.na(E_df$RATING), , drop=FALSE]
      C_const <- C_df[!is.na(C_df$RATING), , drop=FALSE]
      E_mapped <- E_df[ is.na(E_df$RATING), , drop=FALSE]
      C_mapped <- C_df[ is.na(C_df$RATING), , drop=FALSE]

      sum_weighted <- function(df_map) {
        if (!nrow(df_map)) return(zero_r)
        parts <- vector("list", nrow(df_map))
        for (i in seq_len(nrow(df_map))) {
          att <- df_map$ATTRIBUTES[i]
          r   <- rlist[[stressor]][[att]]
          if (is.null(r) || !inherits(r,"SpatRaster")) {
            stop("Missing/invalid raster for '", stressor, "' / attribute '", att, "'.")
          }
          r <- align_to(r, sp_distr, categorical = TRUE)
          parts[[i]] <- r / (df_map$DQ[i] * df_map$WEIGHT[i])
        }
        Reduce(`+`, parts)
      }

      # Build numerators/denominators
      E_numer_const <- if (nrow(E_const)) sum(E_const$RATING / (E_const$DQ * E_const$WEIGHT)) else 0
      C_numer_const <- if (nrow(C_const)) sum(C_const$RATING / (C_const$DQ * C_const$WEIGHT)) else 0
      E_numer_rast <- sum_weighted(E_mapped)
      C_numer_rast <- sum_weighted(C_mapped)
      E_denom <- sum(1 / (E_df$DQ * E_df$WEIGHT))
      C_denom <- sum(1 / (C_df$DQ * C_df$WEIGHT))

      # Union of this stressor’s attributes
      collection <- terra::sprc(rlist[[stressor]])
      r_mos <- terra::mosaic(collection, fun = max)
      stress_occ <- terra::ifel(is.na(r_mos), NA, 1)

      dec_w <- NULL
      E_dec_w <- NULL
      C_dec_w <- NULL

      if (decay %in% c("none","linear","exponential",
                       "polynomial_2nd","polynomial_3rd",
                       "complementary_decay_2nd","complementary_decay_3rd") &&
          !is.null(buffer_m[stressor])) {
        # ensure alignment with species grid
        stress_occ_aligned <- align_to(stress_occ, sp_distr, categorical = TRUE)
        dec_w <- make_stressor_decay(stress_occ = stress_occ_aligned,
                                     buffer_m   = buffer_m[[stressor]],
                                     decay      = decay,
                                     habitat_mask = sp_presence)
        E_dec_w <- dec_w * (E_numer_const / E_denom)
        C_dec_w <- dec_w * (C_numer_const / C_denom)
      }

      # Combine constants + rasters into numerators (rasters)
      E_num <- (E_numer_const + E_numer_rast)
      C_num <- (C_numer_const + C_numer_rast)

      # Final E and C rasters (restricted to habitat pixels)
      E_score_raster <- (E_num / E_denom)
      C_score_raster <- (C_num / C_denom)

      # Apply decay weights and zeros outside buffer
      if (inherits(dec_w, "SpatRaster")) {
        E_score_raster <- terra::cover(E_score_raster, E_dec_w)
        C_score_raster <- terra::cover(C_score_raster, C_dec_w)
      }

      E_score_raster <- terra::mask(E_score_raster, sp_presence)
      C_score_raster <- terra::mask(C_score_raster, sp_presence)

      # Convenience maps with zeros inside presence
      E_map <- terra::cover(E_score_raster, sp_distr_zeros)
      C_map <- terra::cover(C_score_raster, sp_distr_zeros)

      # Risk calculations
      risk_raw <- if (equation=="multiplicative") {
        C_score_raster * E_score_raster
      } else {
        E_opp <- E_score_raster - 1
        C_opp <- C_score_raster - 1
        E_opp <- terra::ifel(E_opp < 0, 0, E_opp)
        C_opp <- terra::ifel(C_opp < 0, 0, C_opp)
        sqrt(E_opp^2 + C_opp^2)
      }

      risk_raw <- terra::cover(risk_raw, sp_distr_zeros)
      cut_off_decay <- 1e-6
      risk_raw <- terra::ifel(risk_raw < cut_off_decay, 0, risk_raw)
      risk_cls <- terra::ifel(
        risk_raw == 0, 0,
        terra::ifel(risk_raw < (1/3)*m_jkl, 1,
                    terra::ifel(risk_raw < (2/3)*m_jkl, 2, 3))
      )

      res[[stressor]] <- list(
        E_criteria = E_map,
        C_criteria = C_map,
        Risk_map_raw = risk_raw,
        Risk_map = risk_cls
      )
    }

    # Summaries for this species
    list_raw <- lapply(res, function(x) x$Risk_map_raw)
    stack_cls <- terra::rast(lapply(res, function(x) x$Risk_map))

    total_raw <- Reduce(`+`, list_raw)
    total_cls <- terra::app(stack_cls, fun = max, na.rm = TRUE)
    total_hotspots_cls <- terra::ifel(
      total_raw == 0, 0,
      terra::ifel(total_raw < (1/3)*m_jkl*n_overlap, 1,
                  terra::ifel(total_raw < (2/3)*m_jkl*n_overlap, 2, 3))
    )

    res$total_raw <- total_raw
    res$total_reclassified <- total_cls
    res$total_hotspots_reclassified <- total_hotspots_cls
    res$summary_stats <- compute_summary_stats_single(res, total_raw, total_hotspots_cls)

    if (isTRUE(output_decimal_crs)) {
      for (st in names(res)) if (is.list(res[[st]])) {
        for (nm in c("E_criteria","C_criteria","Risk_map_raw","Risk_map")) {
          res[[st]][[nm]] <- convert_to_decimal_degrees(res[[st]][[nm]])
        }
      }
      res$total_raw <- convert_to_decimal_degrees(total_raw)
      res$total_reclassified <- convert_to_decimal_degrees(total_cls)
      res$total_hotspots_reclassified <- convert_to_decimal_degrees(total_hotspots_cls)
    }

    class(res) <- c("risaHRA", class(res))
    res
  }

  # Dispatch
  if (depth == 2L) {
    return(.single(
      rlist = raster_list,
      sp_distr = species_distr,
      crit = .check_criteria(criteria),
      equation = equation,
      r_max = r_max,
      n_overlap = n_overlap,
      output_decimal_crs = output_decimal_crs,
      decay,
      buffer_m
    ))
  }

  # Ecosystem mode (depth 3)
  if (depth != 3L) stop("'raster_list' must be depth 2 (single) or 3 (ecosystem).")
  if (!is.list(species_distr) || !is.list(criteria)) {
    stop("In ecosystem mode, 'species_distr' and 'criteria' must be named lists matching raster_list species.")
  }
  species <- names(raster_list)
  if (!setequal(species, names(species_distr)) || !setequal(species, names(criteria))) {
    stop("Species names must match across raster_list, species_distr, and criteria.")
  }

  # Coerce species distributions to SpatRaster and pick a template grid
  sd <- lapply(species, function(sp) .as_raster(species_distr[[sp]]))
  names(sd) <- species
  if (any(!vapply(sd, inherits, logical(1), "SpatRaster"))) {
    stop("All 'species_distr' entries must be or contain a SpatRaster.")
  }

  template <- sd[[1]]

  # Infer n_overlap from union of stressors
  if (is.null(n_overlap)) {
    all_stressors <- unique(unlist(lapply(criteria, function(df) {
      df$STRESSOR[!(is.na(df$STRESSOR) | df$STRESSOR %in% c("", "NA"))]
    })))
    n_overlap <- length(all_stressors)
  }

  # Run single HRA per species
  results <- lapply(species, function(sp) {
    .single(raster_list[[sp]], sd[[sp]], criteria[[sp]],
            equation, r_max, n_overlap, FALSE,
            decay, buffer_m)
  })
  names(results) <- species

  # Calculating ecosystem risk
  # Ecosystem presence mask on template grid (union of species)
  presences <- lapply(sd, function(d) terra::ifel(!is.na(align_to(d, template, TRUE)), 1, 0))
  sum_pres <- Reduce(`+`, presences)
  eco_mask <- terra::ifel(sum_pres > 0, 1, NA)

  # Align per-species general risk (total_raw) to template (continuous)
  rlist <- lapply(species, function(sp) {
    r <- results[[sp]]$total_raw
    align_to(r, template, categorical = FALSE)
  })

  # Make a multilayer SpatRaster
  stk <- terra::rast(rlist)

  # Per-cell sum of risks (ignore NA) and per-cell count of overlapping species
  eco_sum <- terra::app(stk, sum, na.rm = TRUE)
  eco_cnt <- terra::app(stk, function(v) sum(!is.na(v)))

  # Average over overlapping species only (sum / count), keep NA where count==0
  eco_raw <- terra::ifel(eco_cnt > 0, eco_sum / eco_cnt, NA)

  # Ensure zeros inside the union mask
  eco_raw <- terra::cover(eco_raw, terra::ifel(!is.na(eco_mask), 0, NA))

  # Classify ecosystem risk
  eco_cls <- terra::ifel(
    eco_raw == 0, 0,
    terra::ifel(eco_raw < (1/3)*m_jkl*n_overlap, 1,
                terra::ifel(eco_raw < (2/3)*m_jkl*n_overlap, 2, 3))
  )

  # Summary table for ecosystem estimates
  per_species_stats <- do.call(rbind, lapply(names(results), function(sp) {
    cbind.data.frame(SPECIES = sp, results[[sp]]$summary_stats, row.names = NULL)
  }))
  gEco <- terra::global(eco_raw, c("min","max","mean"), na.rm=TRUE)
  fqE <- terra::freq(eco_cls)
  cntE <- c(`0`=0,`1`=0,`2`=0,`3`=0)
  if (!is.null(fqE) && nrow(fqE)) { vv <- as.character(fqE[, "value"]); cntE[vv] <- fqE[, "count"] }
  totE <- sum(cntE)
  eco_row <- data.frame(
    SPECIES="ECOSYSTEM", STRESSOR="(FROM ALL STRESSORS)",
    E_min=NA_real_, E_max=NA_real_, E_mean=NA_real_,
    C_min=NA_real_, C_max=NA_real_, C_mean=NA_real_,
    R_min=as.numeric(gEco[1,"min"]), R_max=as.numeric(gEco[1,"max"]), R_mean=as.numeric(gEco[1,"mean"]),
    `R%high` = if (totE) 100*cntE["3"]/totE else NA_real_,
    `R%medium` = if (totE) 100*cntE["2"]/totE else NA_real_,
    `R%low` = if (totE) 100*cntE["1"]/totE else NA_real_,
    `R%None` = if (totE) 100*cntE["0"]/totE else NA_real_
  )
  summary_stats <- rbind(eco_row, per_species_stats)

  out <- results
  out$ecosys_risk_raw <- eco_raw
  out$ecosys_risk_classified <- eco_cls
  out$summary_stats <- summary_stats

  if (isTRUE(output_decimal_crs)) {
    for (sp in species) {
      for (nm in names(out[[sp]])) if (is.list(out[[sp]][[nm]])) {
        for (k in c("E_criteria","C_criteria","Risk_map_raw","Risk_map")) {
          out[[sp]][[nm]][[k]] <- convert_to_decimal_degrees(out[[sp]][[nm]][[k]])
        }
      } else if (nm %in% c("total_raw", "total_reclassified", "total_hotspots_reclassified")) {
        out[[sp]][[nm]] <- convert_to_decimal_degrees(out[[sp]][[nm]])
      }
    }
    out$ecosys_risk_raw <- convert_to_decimal_degrees(out$ecosys_risk_raw)
    out$ecosys_risk_classified <- convert_to_decimal_degrees(out$ecosys_risk_classified)
  }

  class(out) <- c("risaHRA", class(out))
  out
}
