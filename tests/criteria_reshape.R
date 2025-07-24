path <- "C:/Users/gusma/Documents/risa_maps_test/criteria"

criteria <- read.csv(paste(path, "/test_criteria.csv", sep=""))

df_list <- list()
multi_criteria <- read.csv(paste(path, "/multi_species_criteria.csv", sep=""))
multi_criteria
class(multi_criteria)
sp_n <- (dim(multi_criteria)[2]-2) %/% 3
line_breaks <- which(multi_criteria[,1] %in% c("", NA))
stressor_n <- length(line_breaks)
criteria_column <- dim(multi_criteria)[2]
criteria_column <- ifelse((criteria_column-2) %% 3 > 0, (criteria_column-1), criteria_column)
criteria_column

# Getting stressor names
line_breaks <- which(multi_criteria[,1] %in% c("", NA))
stressor_n <- length(line_breaks)
stressor_names <- c()
for(i in 1:stressor_n){
  if(i == 1) {
    stressor_names[i] <- multi_criteria[line_breaks[i]+2,1]
  } else {
    stressor_names[i] <- multi_criteria[line_breaks[i]+1,1]
  }
}

counter <- 2
for(i in 1:sp_n){
  sub_df <- multi_criteria[,c(1,counter:(counter+2),criteria_column)]
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
names(df_list) <- names(multi_criteria)[seq(2,(sp_n*3), 3)]


