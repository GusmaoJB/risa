#' Reshape a criteria table into a list of data frames respective to each species
#'
#' Splits a standard criteria table for Habitat Risk Assessment
#' (a column with species and stressor attribute descriptors, columns with rating values for each species, and a criteria type column)
#' into a list of data frames (one for each species).  It extracts and recycles stressor names,
#' removes header/blank rows, adds a `STRESSOR` column.
#'
#' @param x A `data.frame` with a column with stressor names,
#' columns with species and stressor attributes (groups of three columns per species),
#' a last column with CRITERIA TYPE` values (E or C). It is separated in row blocks for
#' species attributes and stressor properties by "" or NAs.
#' @return A named `list` of `data.frame` objects, one per species. Each element has:
#' `STRESSOR`, which is a factor of stressor names, and the original attribute columns for that species.
#' @examples
#' #Load example data
#' path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
#' df <- read.csv(path)
#'
#' #Inspect dataframe
#' head(df)
#'
#' #Reshape criteria table
#' crit_list <- criteria_reshape(df)
#' crit_list
#' @export
criteria_reshape <- function(x) {
  if (!is.data.frame(x) || ncol(x) < 6L) {
    stop("`x` must be a data.frame with at least 6 columns.")
  }

  # Detect E/C column (by name or content)
  norm_name <- function(s) tolower(gsub("[^a-z]", "", s))
  name_hits <- which(norm_name(names(x)) %in% c("ec", "criteriatype"))
  crit_col  <- if (length(name_hits)) name_hits[length(name_hits)] else NULL

  if (is.null(crit_col)) {
    is_ec_like <- function(v) {
      tok <- toupper(trimws(as.character(v)))
      tok <- tok[!(is.na(tok) | tok == "")]
      if (!length(tok)) return(FALSE)
      tok <- tok[!tok %in% c("E/C","EC","E C","CRITERIA TYPE")]
      if (!length(tok)) return(TRUE)
      mean(tok %in% c("E","C")) >= 0.80
    }
    cand <- which(vapply(x, is_ec_like, logical(1)))
    if (length(cand)) crit_col <- cand[length(cand)]
  }
  if (is.null(crit_col)) stop("Could not find a column with only 'E'/'C' values (criteria type).")

  # Species triplets live between column 1 and E/C column
  mid_cols <- (crit_col - 1L) - 1L
  if (mid_cols < 3L) stop("Not enough species columns between the first column and the E/C column.")
  sp_n <- floor(mid_cols / 3L)
  sp_firsts <- 2L + 3L * (0:(sp_n - 1L))
  sp_names <- names(x)[sp_firsts]
  if (is.null(sp_names) || any(!nzchar(sp_names))) sp_names <- paste0("sp", seq_len(sp_n))

  # Build STRESSOR labels from blank-separated blocks in col 1
  lab <- as.character(x[[1L]])
  is_blank <- is.na(lab) | trimws(lab) == ""
  n <- nrow(x)
  blank_idx <- which(is_blank)

  # Vectorized section-title detector
  is_section_title <- function(s) {
    s0 <- tolower(trimws(as.character(s)))
    out <- grepl("resilience.*attribute", s0) |
      grepl("stressor.*overlap.*propert", s0) |
      grepl("^rating instruction", s0)
    out[is.na(out)] <- FALSE
    out
  }

  header_idx <- integer(0)
  header_name <- character(0)

  for (b in blank_idx) {
    j <- b + 1L
    # Advance to first non-blank
    while (j <= n && (is.na(lab[j]) || trimws(lab[j]) == "")) j <- j + 1L
    if (j > n) next
    # Skip ONLY section titles; DO NOT skip E/C-labeled stressor header rows
    while (j <= n && is_section_title(lab[j])) j <- j + 1L
    if (j <= n && !is_blank[j]) {
      header_idx <- c(header_idx, j)
      header_name <- c(header_name, lab[j])
    }
  }

  # Assign each row after a header to that stressor until next blank
  stressor_by_row <- rep(NA_character_, n)
  if (length(header_idx)) {
    for (k in seq_along(header_idx)) {
      start <- header_idx[k] + 1L
      nb <- blank_idx[blank_idx > header_idx[k]]
      end <- if (length(nb)) nb[1L] - 1L else n
      if (start <= end) stressor_by_row[start:end] <- header_name[k]
    }
  }

  # Helpers to detect internal header lines
  looks_like_internal_header <- function(row_vals) {
    vv <- tolower(trimws(as.character(row_vals)))
    any(grepl("^rating$|^dq$|^weight$|^e/?c$|^criteria", vv))
  }
  is_row_internal_header <- function(df_row) {
    ec_val <- tolower(trimws(as.character(df_row[[length(df_row)]])))
    has_eclabel <- ec_val %in% c("e/c","ec","e c","criteria type")
    has_tokens  <- looks_like_internal_header(df_row)
    has_eclabel || has_tokens
  }

  # Per-species reshape
  out <- vector("list", sp_n)

  for (i in seq_len(sp_n)) {
    start <- 2L + 3L * (i - 1L)
    cols <- c(1L, start:(start + 2L), crit_col)
    sub <- x[, cols, drop = FALSE]

    names(sub) <- c("ATTRIBUTES","RATING","DQ","WEIGHT","E/C")
    r1 <- as.character(sub[1L, , drop = TRUE])
    data_rows <- if (looks_like_internal_header(r1)) 2L:nrow(sub) else 1L:nrow(sub)

    # Drop blank rows, section titles, and any internal header-like rows
    attrv <- as.character(sub$ATTRIBUTES)
    sec_title_rows <- is_section_title(attrv)
    internal_header_rows <- vapply(seq_len(nrow(sub)), function(rr) {
      is_row_internal_header(sub[rr, , drop = FALSE])
    }, logical(1))

    drop_rows <- which(is.na(attrv) | trimws(attrv) == "" | sec_title_rows | internal_header_rows)
    keep <- setdiff(data_rows, drop_rows)

    if (!length(keep)) {
      warning("No data rows found for species ", sQuote(sp_names[i]), ". Returning empty table.")
      out[[i]] <- data.frame(STRESSOR=character(0), ATTRIBUTES=character(0),
                             RATING=integer(0), DQ=integer(0), WEIGHT=integer(0),
                             `E/C`=character(0))
      next
    }

    df <- sub[keep, , drop = FALSE]
    df$STRESSOR <- stressor_by_row[keep]

    suppressWarnings({
      df$RATING <- as.integer(df$RATING)
      df$DQ <- as.integer(df$DQ)
      df$WEIGHT <- as.integer(df$WEIGHT)
    })

    df <- df[, c("STRESSOR","ATTRIBUTES","RATING","DQ","WEIGHT","E/C")]
    rownames(df) <- NULL
    out[[i]] <- df
  }

  names(out) <- sp_names
  out
}
