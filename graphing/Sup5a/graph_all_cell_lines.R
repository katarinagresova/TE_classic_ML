library(ggplot2)
library(cowplot)

data <- read.csv("./results/human/all_cell_lines/lgbm-LL_P5_P3_CF_AAF_3mer_freq_5/model_results.csv", 
                 header = TRUE, sep = ","
)
save_path <- "./graphing/B/graph_all_cell_lines_human.pdf"
# data <- read.csv("./results/mouse/all_cell_lines/lgbm-LL_P5_P3_CF_AAF_3mer_freq_5_Struct/model_results.csv", 
#                  header = TRUE, sep = ","
# )
# save_path <- "./graphing/B/graph_all_cell_lines_mouse.pdf"
# drop all rows that have mean in fold column
data <- data[!grepl("mean", data$fold),]
# in bio_source col remove bio_source_ prefix
data$bio_source <- gsub("bio_source_", "", data$bio_source)


data$highlight <- ifelse(
    data$bio_source == "mean_across_cell_lines", # data$bio_source == "mean_te" |
    "mean_te",
    "NoHighlight"
  )


# biochem_data <- read.csv("./results/all_cell_lines/lgbm-L_P5_P3_C_3mer5_Struct_Biochem/model_results.csv", 
#                  header = TRUE, sep = ","
# )
# # drop all rows that have mean in fold column
# biochem_data <- biochem_data[!grepl("mean", biochem_data$fold),]
# # in bio_source col remove bio_source_ prefix
# biochem_data$bio_source <- gsub("bio_source_", "", biochem_data$bio_source)


# biochem_data$highlight <- ifelse(
#     biochem_data$bio_source == "mean_te",
#     "mean_te",
#     "NoHighlight"
#   )

# round to nearest 0.1
# min_r2 <- floor(min(data$nan_r2) * 10) / 10
# max_r2 <- ceiling(max(data$nan_r2) * 10) / 10
# round to nearest 0.25
# min_r2 = floor(min(data$nan_r2) * 4) / 4
min_r2 <- 0
max_r2 <- ceiling(max(data$nan_r2) * 4) / 4
# min_r2 <- min(min(data$nan_r2), min(biochem_data$r2_score))
# max_r2 <- max(max(data$nan_r2), max(biochem_data$r2_score))

mean_with_se_0 <- function(x) {
  m <- mean(x)
  return(data.frame(y = m, ymin = m, ymax = m))
}

mean_min_max <- function(x) {
  return(data.frame(y = mean(x), ymin = min(x), ymax = max(x)))
}

p <- ggplot(data, aes(x = nan_r2, y = reorder(bio_source, nan_r2, FUN = mean), fill = highlight)) +
  geom_pointrange(stat = "summary", fun.data = mean_min_max, aes(color = highlight), size = 0.1) +
  # geom_point(aes(x = nan_r2, color = highlight), size = 0.1) +

  labs(x = expression("10-fold CV R"^2), y = NULL, title = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 4),
    axis.text.x = element_text(size = 6),
    # axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    title = element_text(size = 7)
  ) +
  coord_cartesian(xlim = c(min_r2, max_r2)) + 
  scale_x_continuous(breaks = seq(min_r2, max_r2, by = 0.25))

# save plot as png
# ggsave("./graphing/graph_all_cell_lines.pdf", p, width = 5, height = 5, units = "in", dpi = 300, bg = 'white')
ggsave(save_path, p, width = 3, height = 5, units = "in", dpi = 300, bg = 'white')
# print(p)

# Do the same for Biochem model


# biochem_p <- ggplot(biochem_data, aes(x = r2_score, y = reorder(bio_source, data$r2_score, FUN = mean), fill = highlight)) +
#   geom_bar(stat = "summary", fun = "mean") +
#   scale_fill_manual(
#     name = "Model",
#     values = c("mean_te" = "#0dff00")
#   ) +
#   geom_point(aes(x = r2_score), color = "black", size = 0.25) +
# #   geom_errorbar(stat = "summary", fun.data = mean_se, color = "black", width = 0.2) +

#   labs(x = "R2 Score", y = NULL, title = "L_P5_P3_C_3mer5_Struct_Biochem") +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_text(size = 4),
#     # axis.text.x = element_text(angle = 0, hjust = 1),
#     legend.position = "none",
#     title = element_text(size = 7)
#   ) +
#   coord_cartesian(xlim = c(min_r2, max_r2))

# # save plot as png
# ggsave("./graphing/graph_all_cell_lines_Biochem.pdf", biochem_p, width = 5, height = 5, units = "in", dpi = 300, bg = 'white')
# # print(p)

# legend <- get_legend(p + theme(legend.position = "bottom", legend.title = element_blank()))

# # combine the two plots and remove the y axis from the second plot
# combined_plot <- plot_grid(p, biochem_p, ncol = 2)
# combined_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))
# # save plot as png
# ggsave("./graphing/graph_all_cell_lines_combined.pdf", combined_plot, width = 10, height = 5, units = "in", dpi = 300, bg = 'white')