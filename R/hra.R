#' Habitat Risk Assessment (HRA): single species or ecosystem
#'
#' Single-species: `raster_list` is a depth-2 named list: stressor -> attribute rasters.
#' Ecosystem: `raster_list` is depth-3: species -> stressor -> attribute rasters.
#'
#' @param raster_list list
#' @param species_distr SpatRaster (single) OR named list of SpatRasters (ecosystem).
#'        If an element is a list, the first SpatRaster within it will be used.
#' @param criteria data.frame (single) OR named list of data.frames (ecosystem).
#' @param equation c("euclidean","multiplicative")
#' @param r_max integer in 1..10
#' @param n_overlap optional number of stressors (default inferred)
#' @param output_decimal_crs logical; reproject outputs to EPSG:4326 if TRUE
#' @return "risaHRA" list with per-stressor maps, totals, and summary_stats
#' @importFrom terra compareGeom project mosaic mask ifel freq global
#' @examples
#' # example code
#'
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
    r_max = 3, n_overlap = NULL, output_decimal_crs = FALSE
) {
  depth <- list_depth_base(raster_list)
  equation <- match.arg(equation)

  # Helpers
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

  .align_to <- function(src, template, categorical = TRUE) {
    if (terra::compareGeom(src, template, stopOnError = FALSE)) return(src)
    m <- if (categorical) "near" else "bilinear"
    terra::project(src, template, method = m)
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

  .compute_summary_stats_single <- function(rr, total_risk_raw, total_class) {
    stressor_names <- setdiff(names(rr), c("total_raw","total","summary_stats"))
    rows <- lapply(stressor_names, function(st) {
      E <- rr[[st]]$E_criteria
      C <- rr[[st]]$C_criteria
      Rr <- rr[[st]]$Risk_map_raw
      Rc <- rr[[st]]$Risk_map
      gE <- terra::global(E,  c("min","max","mean"), na.rm=TRUE)
      gC <- terra::global(C,  c("min","max","mean"), na.rm=TRUE)
      gR <- terra::global(Rr, c("min","max","mean"), na.rm=TRUE)

      fq <- terra::freq(Rc)
      cnt <- c(`0`=0,`1`=0,`2`=0,`3`=0)
      if (!is.null(fq) && nrow(fq)) {
        vv <- as.character(fq[, "value"]); cnt[vv] <- fq[, "count"]
      }
      tot <- sum(cnt)

      data.frame(
        STRESSOR = st,
        E_min  = as.numeric(gE[1,"min"]), E_max = as.numeric(gE[1,"max"]), E_mean = as.numeric(gE[1,"mean"]),
        C_min  = as.numeric(gC[1,"min"]), C_max = as.numeric(gC[1,"max"]), C_mean = as.numeric(gC[1,"mean"]),
        R_min  = as.numeric(gR[1,"min"]), R_max = as.numeric(gR[1,"max"]), R_mean = as.numeric(gR[1,"mean"]),
        `R%high`   = if (tot) 100*cnt["3"]/tot else NA_real_,
        `R%medium` = if (tot) 100*cnt["2"]/tot else NA_real_,
        `R%low`    = if (tot) 100*cnt["1"]/tot else NA_real_,
        `R%None`   = if (tot) 100*cnt["0"]/tot else NA_real_
      )
    })
    per <- do.call(rbind, rows)

    fqt <- terra::freq(total_class)
    cntt <- c(`0`=0,`1`=0,`2`=0,`3`=0)
    if (!is.null(fqt) && nrow(fqt)) {
      vv <- as.character(fqt[, "value"]); cntt[vv] <- fqt[, "count"]
    }
    tot <- sum(cntt)

    overall <- data.frame(
      STRESSOR="(FROM ALL STRESSORS)",
      E_min=min(per$E_min,na.rm=TRUE), E_max=max(per$E_max,na.rm=TRUE), E_mean=mean(per$E_mean,na.rm=TRUE),
      C_min=min(per$C_min,na.rm=TRUE), C_max=max(per$C_max,na.rm=TRUE), C_mean=mean(per$C_mean,na.rm=TRUE),
      R_min=min(per$R_min,na.rm=TRUE), R_max=max(per$R_max,na.rm=TRUE), R_mean=mean(per$R_mean,na.rm=TRUE),
      `R%high` = if (tot) 100*cntt["3"]/tot else NA_real_,
      `R%medium` = if (tot) 100*cntt["2"]/tot else NA_real_,
      `R%low` = if (tot) 100*cntt["1"]/tot else NA_real_,
      `R%None` = if (tot) 100*cntt["0"]/tot else NA_real_
    )
    rbind(overall, per)
  }

  .single <- function(rlist, sp_distr, crit, equation, r_max, n_overlap, output_decimal_crs) {
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

    sp_presence    <- terra::ifel(!is.na(sp_distr), 1, NA)
    sp_distr_zeros <- terra::ifel(!is.na(sp_distr), 0, NA)
    zero_r         <- sp_distr * 0

    res <- setNames(vector("list", length(names(rlist))), names(rlist))
    m_jkl <- if (equation=="multiplicative") r_max^2 else sqrt(2*(r_max-1)^2)

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

      E_const  <- E_df[!is.na(E_df$RATING), , drop=FALSE]
      C_const  <- C_df[!is.na(C_df$RATING), , drop=FALSE]
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
          r <- .align_to(r, sp_distr, categorical = TRUE)
          parts[[i]] <- r / (df_map$DQ[i] * df_map$WEIGHT[i])
        }
        Reduce(`+`, parts)
      }

      E_numer_const <- if (nrow(E_const)) sum(E_const$RATING / (E_const$DQ * E_const$WEIGHT)) else 0
      C_numer_const <- if (nrow(C_const)) sum(C_const$RATING / (C_const$DQ * C_const$WEIGHT)) else 0
      E_numer_rast  <- sum_weighted(E_mapped)
      C_numer_rast  <- sum_weighted(C_mapped)

      E_denom <- sum(1 / (E_df$DQ * E_df$WEIGHT))
      C_denom <- sum(1 / (C_df$DQ * C_df$WEIGHT))

      E_score_raster <- (E_numer_const + E_numer_rast) / E_denom
      C_score_raster <- (C_numer_const + C_numer_rast) / C_denom

      E_map <- terra::mosaic(E_score_raster, sp_distr_zeros, fun="first")
      C_map <- terra::mosaic(terra::mask(C_score_raster, sp_distr), sp_distr_zeros, fun="first")

      risk_raw <- if (equation=="multiplicative") {
        (C_score_raster * E_score_raster) * sp_presence
      } else {
        sqrt((E_score_raster - 1)^2 + (C_score_raster - 1)^2) * sp_presence
      }
      risk_raw <- terra::mosaic(risk_raw, sp_distr_zeros, fun="first")
      risk_cls <- terra::ifel(
        risk_raw == 0, 0,
        terra::ifel(risk_raw < (1/3)*m_jkl, 1,
                    terra::ifel(risk_raw < (2/3)*m_jkl, 2, 3))
      )

      res[[stressor]] <- list(
        E_criteria   = E_map,
        C_criteria   = C_map,
        Risk_map_raw = risk_raw,
        Risk_map     = risk_cls
      )
    }

    total_raw <- Reduce(`+`, lapply(res, function(x) x$Risk_map_raw))
    total_cls <- terra::ifel(
      total_raw == 0, 0,
      terra::ifel(total_raw < (1/3)*m_jkl*n_overlap, 1,
                  terra::ifel(total_raw < (2/3)*m_jkl*n_overlap, 2, 3))
    )

    res$total_raw     <- total_raw
    res$total         <- total_cls
    res$summary_stats <- .compute_summary_stats_single(res, total_raw, total_cls)

    if (isTRUE(output_decimal_crs)) {
      for (st in names(res)) if (is.list(res[[st]])) {
        for (nm in c("E_criteria","C_criteria","Risk_map_raw","Risk_map")) {
          res[[st]][[nm]] <- convert_to_decimal_degrees(res[[st]][[nm]])
        }
      }
      res$total_raw <- convert_to_decimal_degrees(total_raw)
      res$total     <- convert_to_decimal_degrees(total_cls)
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
      output_decimal_crs = output_decimal_crs
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
            equation, r_max, n_overlap, FALSE)
  })
  names(results) <- species

  # Ecosystem presence mask on template grid (union of species)
  presences <- lapply(sd, function(d) terra::ifel(!is.na(.align_to(d, template, TRUE)), 1, 0))
  sum_pres  <- Reduce(`+`, presences)
  eco_mask  <- terra::ifel(sum_pres > 0, 1, NA)

  # Ecosystem risk: average of species' total_raw, aligned to template
  m_jkl <- if (equation=="multiplicative") r_max^2 else sqrt(2*(r_max-1)^2)
  eco_raw <- template * 0
  for (sp in species) {
    r <- results[[sp]]$total_raw
    r <- .align_to(r, template, categorical = FALSE) # continuous
    r <- terra::ifel(is.na(r), 0, r)
    eco_raw <- eco_raw + r
  }
  eco_raw <- eco_raw / length(species)
  eco_raw <- terra::mosaic(eco_raw, terra::ifel(!is.na(eco_mask), 0, NA), fun="first")
  eco_cls <- terra::ifel(
    eco_raw == 0, 0,
    terra::ifel(eco_raw < (1/3)*m_jkl*n_overlap, 1,
                terra::ifel(eco_raw < (2/3)*m_jkl*n_overlap, 2, 3))
  )

  # Summary table
  per_species_stats <- do.call(rbind, lapply(names(results), function(sp) {
    cbind.data.frame(SPECIES = sp, results[[sp]]$summary_stats, row.names = NULL)
  }))
  gEco <- terra::global(eco_raw, c("min","max","mean"), na.rm=TRUE)
  fqE  <- terra::freq(eco_cls)
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
      } else if (nm %in% c("total_raw","total")) {
        out[[sp]][[nm]] <- convert_to_decimal_degrees(out[[sp]][[nm]])
      }
    }
    out$ecosys_risk_raw <- convert_to_decimal_degrees(out$ecosys_risk_raw)
    out$ecosys_risk_classified <- convert_to_decimal_degrees(out$ecosys_risk_classified)
  }

  class(out) <- c("risaHRA", class(out))
  out
}
