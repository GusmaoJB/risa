#’ Reshape a criteria table into species‐specific data frames
#’
#’ Splits a criteria table—where the first column lists stressors (with blank rows as separators)
#’ and the remaining columns are organized in triplets per species plus a final rating column—
#’ into a list of data frames, one for each species. Extracts and recycles stressor names,
#’ removes header/blank rows, adds a `STRESSOR` column, and converts `RATING` to integer.
#’
#’ @param x A `data.frame` with:
#’   * Column 1: stressor names, separated by blank (`""`) or `NA` rows
#’   * Columns 2–(3*n+1): groups of three columns per species
#’   * Last column: `RATING` values
#’ @return A named `list` of `data.frame` objects, one per species. Each element has:
#’   * `STRESSOR`: factor of stressor names
#’   * the original attribute columns for that species
#’   * `RATING`: integer rating
#’ @examples
#’ \donttest{
#’ # Suppose you have a CSV at inst/extdata/criteria_example.csv matching the expected layout:
#’ path    <- system.file("extdata", "criteria_example.csv", package = "risa")
#’ criteria <- read.csv(path, stringsAsFactors = FALSE)
#’ result   <- criteria_reshape(criteria)
#’ str(result[[1]])  # inspect the first species
#’ }
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



path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
path

df   <- read.csv(path, stringsAsFactors = FALSE)

path <- "C:/Users/gusma/Documents/risa_maps_test/criteria"
criteria <- read.csv(paste(path, "/test_criteria.csv", sep=""))
multi_criteria <- read.csv(paste(path, "/multi_species_criteria.csv", sep=""))


crit_list <- criteria_reshape(multi_criteria)
crit_list
