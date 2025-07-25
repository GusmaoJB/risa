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
#' #test
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
