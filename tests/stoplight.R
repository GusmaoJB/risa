library(risa)
library(ggplot2)

# Example data
quality_data <- cbind.data.frame(
  level	= c("species",	"species", "species", "species",
            "stressor",	"stressor",
            "interaction", "interaction",	"interaction",	"interaction"),
  data = c("Animal sightings distribution",	"Habitat suitability",
           "Animal sightings distribution",	"Habitat suitability", "Fishing effort",
           "Fishing effort",	"Bycatch/stranding data",	"Bycatch/stranding data",
           "Bycatch/stranding data",	"Bycatch/stranding data"),
  species	= c("dolphin", "dolphin", "turtle", "turtle", NA, NA, "dolphin",
              "turtle", "dolphin", "turtle"),
  stressor = c(NA,	NA,	NA,	NA,	"trawling",	"gillnet",	"trawling",	"trawling",	"gillnet",	"gillnet"),
  uncertainty	= c(3, 2, 2, 2, 2, 3, 1, 2, 1, 2)
)

# Construct a function that converts a quality data data.frame into three
# data.frames with color codes.

# score_to_colors <- function(quality.data.frame)

names(quality_data) <- c("level", "data", "species", "stressor", "uncertainty")

data_types <- na.omit(unique(quality_data$data))
spp <- na.omit(unique(quality_data$species))
str <- na.omit(unique(quality_data$stressor))

# Colorblind friendly stoplight colors: dark red (#8B0000), saturated yellow (#FFFF00),
# and teal green (#009E73). These colors translate to the gray scale, respectively,
# as very dark gray (~#2A2A2A), very light gray (~#E2E2E2), and medium gray (~#6A6A6A)

# Helpers

# Round half up
round_half_up <- function(x, digits = 0) {
  posneg <- sign(x)
  z <- abs(x) * 10^digits
  z <- z + 0.5 + sqrt(.Machine$double.eps)
  z <- trunc(z)
  z <- z / 10^digits
  z * posneg
}

# Color classifier
stoplight_cols <- function(x) {
  if(!is.numeric(x)) {
    stop("Error: input must be a number ranging from 0 to 3")
  }
  if(x > 3 | x < 0) {
    stop("Error: input must be a number ranging from 0 to 3")
  }
  r_x <- round_half_up(x)
  ifelse(r_x == 0 | is.na(r_x), "#000000",
         ifelse(r_x == 1, "#8B0000",
                ifelse(r_x == 2, "#FFFF00", "#009E73")))
}

# Null checker operator
`%||%` <- function(a, b) if (is.null(a) || is.na(a) || a == "") b else a

# Ecosystem stoplight
ecosys_risk_df <- aggregate(uncertainty ~ data, data=quality_data, FUN=mean)
ecosys_risk_cols <- ecosys_risk_df
names(ecosys_risk_cols) <- c("light", "value")

for (i in 1:length(ecosys_risk_cols$value)) {
  ecosys_risk_cols$value[i] <- stoplight_cols(ecosys_risk_df$uncertainty[i])
}

# Species stoplight
spp_risk_df <- data.frame()

for (sp in spp) {
  sp_subset <- quality_data[quality_data$species == sp |
                            is.na(quality_data$species),]
  sp_subset$species <- sp

  spp_risk_df <- rbind.data.frame(spp_risk_df, sp_subset)
}

spp_risk_df_2 <- na.omit(spp_risk_df)

for (sp in spp) {
  sp_subset <- quality_data[!is.na(quality_data$species),]
  sp_subset <- sp_subset[sp_subset$species == sp ,]
  for (st in str) {
    st_subset <- sp_subset[is.na(sp_subset$stressor),]
    st_subset$stressor <- st
    spp_risk_df_2 <- rbind.data.frame(spp_risk_df_2, st_subset)
  }
}

spp_risk_vals <- spp_risk_df_2[,c("species", "stressor", "data", "uncertainty")]
names(spp_risk_vals) <- c("species", "stressor", "light", "value")
spp_risk_cols <- spp_risk_vals

for (i in 1:length(spp_risk_cols$value)) {
  spp_risk_cols$value[i] <- stoplight_cols(spp_risk_vals$value[i])
}


# Species with combined stressors

spp_comb_risk_vals <- aggregate(value ~ species + light, data=spp_risk_vals, FUN=mean)
spp_comb_risk_cols <- spp_comb_risk_vals

for (i in 1:length(spp_comb_risk_cols$value)) {
  spp_comb_risk_cols$value[i] <- stoplight_cols(spp_comb_risk_vals$value[i])
}

# it must return a list with these dataframes:
list_color_dfs <- list(
  ecossystem = ecosys_risk_cols,
  species_per_stressor = spp_risk_cols,
  species_comb_stressor = spp_comb_risk_cols
)













spp_risk_spp <- aggregate(uncertainty ~ data + species, data = dat_comb_spp, FUN = mean)

spp_risk_cols <- list()

for (sp in spp) {
  subset_sp <- spp_risk_spp[spp_risk_spp$species %in% c(sp, "NA"),]
  sp_risk <- setNames(subset_sp$uncertainty, subset_sp$data)
  sp_colors <- sp_risk

  for (i in sp_risk) {
    sp_colors[i] <- stoplight_cols(sp_risk[i])
  }
  spp_risk_cols[[sp]] <- sp_colors

}

spp_risk_cols


cols

quality_data

# your real plot
p <- ggplot(mtcars, aes(wt, mpg)) +
  geom_point()

# data for the legend you want to show, not otherwise plotted
legend_0 <- data.frame(x = mean(mtcars$wt), y = mean(mtcars$mpg),
                        category = c("Low", "Medium", "High"))
legend_1 <- data.frame(x = mean(mtcars$wt), y = mean(mtcars$mpg),
                       category = c("Low", "Low", "Medium"))

p +
  geom_point(data = legend_0, aes(x, y, color = category), alpha = 0) +
  scale_color_manual(
    name = "Custom Legend",
    values = c("Low" = "blue", "Medium" = "orange", "High" = "red")
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme(legend.key = element_rect(fill = "gray20")) +
  facet_wrap(. ~ vs)
















