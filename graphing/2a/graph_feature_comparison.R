library(ggplot2)
library(cowplot)
# library(patchwork)
# library(tidyr)
# library(dplyr)

data <- read.csv("./results/feature_set_comparison/summary.csv", header = TRUE, sep = ",")

data$model <- gsub("on_NA_lgbm-", "", data$model)


data$highlight <- ifelse(
  data$model == "ALL_FEATURES",
  "None",
  ifelse(
    data$model == "LL_P5_P3_CF_AAF_3mer_freq_5",
    "Best Seq Only",
    "None"
    # ifelse(
    #   # if model contains Biochem
    #   grepl("Biochem", data$model),
    #   "w/ Biochem",
    #   "None"
    # )
  )
)

features_included <- data.frame(model = unique(data$model))
# "Yes" if features_included$model contains LL 
features_included$LL <- ifelse(grepl("LL", features_included$model), "Yes", "No")
features_included$UTR_NT_PERCENT <- ifelse(grepl("P5_P3", features_included$model), "Yes", "No")
features_included$CODON_FREQ <- ifelse(grepl("CF", features_included$model), "Yes", "No")
features_included$AA_FREQ <- ifelse(grepl("AAF", features_included$model), "Yes", "No")
features_included$UTR5_3MER_FREQ <- ifelse(grepl("3mer_freq_5", features_included$model), "Yes", "No")
features_included$STRUCTURE <- ifelse(grepl("Struct", features_included$model), "Yes", "No")
features_included$BIOCHEM <- ifelse(grepl("Biochem", features_included$model), "Yes", "No")

features_included$OTHER <- "No"
# for every col except model col, make ALL_FEATURES "Yes"
for (col in colnames(features_included)) {
  if (col != "model") {
    features_included[features_included$model == "ALL_FEATURES", col] <- "Yes"
  }
}

# mean_with_se_0 <- function(x) {
#   m <- mean(x)
#   return(data.frame(y = m, ymin = m, ymax = m))
# }

median_with_se_0 <- function(x) {
  m <- median(x)
  return(data.frame(y = m, ymin = m, ymax = m))
}

# round down to nearest 0.1
min_r2 = floor(min(data$nan_r2) * 10) / 10
max_r2 = ceiling(max(data$nan_r2) * 10) / 10

p1 <- ggplot(data, aes(y = nan_r2, x = reorder(model, nan_r2, FUN = mean), fill = highlight)) +
  geom_errorbar(stat = "summary", fun.data = median_with_se_0, aes(color = highlight), width = 0.75, linewidth = 1.25) +
  geom_point(aes(y = nan_r2), color = "black", size = 0.75) +
  # scale_fill_manual(
  #   name = "Model",
  #   values = c("All" = "#0dff00", "Best Seq Only" = "#0091ff", "w/ Biochem" = "red")
  # ) +
  labs(y = "R2 Score", x = NULL, title = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 3, angle = 90, hjust = 1),
    # axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none",
  ) +
  coord_cartesian(ylim = c(min_r2, max_r2)) +
  scale_y_continuous(breaks = seq(min_r2, max_r2, by = 0.1))



legend <- get_legend(p1 + theme(legend.position = "bottom"))

p <- plot_grid(p1, legend, ncol = 1, rel_heights = c(1, 0.1))


# save data to csv
ggsave("./graphing/A/graph_feature_comparison.pdf", p, width = 5, height = 3.5, units = "in", bg = 'white')
# print(p)
