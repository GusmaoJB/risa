library(risa)
library(ggplot2)
library(cowplot)
library(patchwork)

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

names(quality_data) <- c("level", "data", "species", "stressor", "uncertainty")

data_types <- na.omit(unique(quality_data$data))
spp <- na.omit(unique(quality_data$species))
str <- na.omit(unique(quality_data$stressor))

# Colorblind friendly stoplight colors: dark red (#8B0000), saturated yellow (#FFFF00),
# and teal green (#009E73). These colors translate to the gray scale, respectively,
# as very dark gray (~#2A2A2A), very light gray (~#E2E2E2), and medium gray (~#6A6A6A)
# Helpers
# Color classifier
stoplight_cols <- function(x) {
  if(!is.numeric(x)) {
    stop("Error: input must be a number ranging from 0 to 3")
  }
  if(x > 3 | x < 0) {
    stop("Error: input must be a number ranging from 0 to 3")
  }
  r_x <- round(x)
  ifelse(r_x == 0 | is.na(r_x), "#000000",
         ifelse(r_x == 1, "#8B0000",
                ifelse(r_x == 2, "#FFFF00", "#009E73")))
}

# Null checker operator
`%||%` <- function(a, b) if (is.null(a) || is.na(a) || a == "") b else a

# Make ggplot legendts
make_legend <- function(color_map, legend_name = NULL, point_size = 3, key_fill = NA) {

  mock <- ggplot(df_unique, aes(x = 1, y = 1, color = category)) +
    geom_point() +
    scale_color_manual(
      name   = legend_name %||% "",
      values = color_map
    ) +
    guides(color = guide_legend(override.aes = list(size = point_size, alpha = 1)))

  if (!is.na(key_fill)) {
    mock <- mock + theme(legend.key = element_rect(fill = key_fill))
  }

  leg <- get_legend(mock)
  wrap_elements(leg)
}

# Ecosystem stoplight
ecosys_risk_df <- aggregate(uncertainty ~ data, data=quality_data, FUN=mean)
ecosys_risk <- setNames(ecosys_risk_df$uncertainty, ecosys_risk_df$data)
ecosys_risk_cols <- ecosys_risk

for (i in 1:length(ecosys_risk)) {
  ecosys_risk_cols[i] <- stoplight_cols(ecosys_risk[i])
}

# Species stoplight - combined stressors
head(quality_data)
dat_comb_spp <- quality_data
dat_comb_spp[is.na(dat_comb_spp)] <- "NA"
dat_comb_spp
spp_risk_spp <- aggregate(uncertainty ~ data + species, data = dat_comb_spp, FUN = mean)
spp_risk_spp <- spp_risk_spp[spp_risk_spp$species %in% c(spp, "NA"),]

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
















