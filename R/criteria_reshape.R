#’ Reshape a criteria table into a list of species‐specific data frames
#’
#’ Splits a criteria table, where the first column lists species attributes and stressors properties (with blank rows as separators)
#’ and the remaining columns are organized in triplets, with RATING, DQ, and WEIGHT per species plus a final CRITEIRA TYPE (E/C) column—
#’ into a list of data frames, one for each species.  It extracts and recycles stressor names,
#’ removes header/blank rows, adds a `STRESSOR` column.
#’
#’ @param x A `data.frame` with:
#’   * Column 1: stressor names, separated by blank (`""`) or `NA` rows
#’   * Columns with species and stressor attributes (groups of three columns per species)
#’   * Last column: `CRITERIA TYPE` values (E or C)
#’ @return A named `list` of `data.frame` objects, one per species. Each element has:
#’   * `STRESSOR`: factor of stressor names
#’   * the original attribute columns for that species
#’   * `RATING`: integer rating
#’ @examples
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
#’ @export
criteria_reshape <- function(x) {
  df_list <- list()
  sp_n <- (dim(x)[2]-2) %/% 3
  line_breaks <- which(x[,1] %in% c("", NA))
  stressor_n <- length(line_breaks)
  criteria_column <- dim(x)[2]
  criteria_column <- ifelse((criteria_column-2) %% 3 > 0, (criteria_column-1), criteria_column)

  # Getting stressor names
  line_breaks <- which(x[,1] %in% c("", NA))
  stressor_n <- length(line_breaks)
  stressor_names <- c()
  for(i in 1:stressor_n){
    if(i == 1) {
      stressor_names[i] <- x[line_breaks[i]+2,1]
    } else {
      stressor_names[i] <- x[line_breaks[i]+1,1]
    }
  }

  # Spliting dataframes
  counter <- 2
  for(i in 1:sp_n){
    sub_df <- x[,c(1,counter:(counter+2),criteria_column)]
    str_attr_n <- (dim(sub_df)[1] - line_breaks[1] - stressor_n*2)/stressor_n
    spe_attr <- rep(c("NA"), line_breaks[1]-2)
    str_attr <- rep(stressor_names, each = str_attr_n)
    names(sub_df) <- sub_df[1,]
    skip_rows <- c(1, line_breaks[1]:(line_breaks[1]+2))
    if (length(line_breaks) > 1) {
      skip_rows <- c(1, line_breaks[1]:(line_breaks[1]+2), line_breaks[-1], (line_breaks[-1]+1))
    }
    sub_df <- sub_df[-c(skip_rows),]
    sub_df <- cbind.data.frame(STRESSOR = c(spe_attr, str_attr), sub_df)
    sub_df$RATING <- as.integer(sub_df$RATING)
    df_list[[i]] <- sub_df
    counter <- counter+2
  }
  names(df_list) <- names(x)[seq(2,(sp_n*3), 3)]
  return(df_list)
}
