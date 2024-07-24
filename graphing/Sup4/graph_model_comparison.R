library(ggplot2)

data <- read.csv("./results/model_compare/summary.csv", 
                 header = TRUE, sep = ","
)
# from model name remove everything after -
data$model <- gsub("-.*", "", data$model)
data$model <- gsub("on_NA_", "", data$model)



data$highlight <- ifelse(
    data$model == "lgbm",
    "Best",
    ""
  )

# Custom function to calculate mean and set se to 0
mean_with_se_0 <- function(x) {
  m <- mean(x)
  return(data.frame(y = m, ymin = m, ymax = m))
}

# round down to nearest 0.5
min_r2 <- 0.0
max_r2 <- 0.65

p <- ggplot(data, aes(x = nan_r2, y = reorder(model, nan_r2, FUN = mean), fill = highlight)) +
  # geom_bar(stat = "summary", fun = "mean") +
  # scale_fill_manual(
  #   name = "Model",
  #   values = c("RedHighlight" = "red", "BlueHighlight" = "#0091ff", "NoHighlight" = "grey50")
  # ) +
  geom_point(aes(x = nan_r2), color = "black", size = 1) +
  geom_errorbar(stat = "summary", fun.data = mean_with_se_0, aes(color = 'red'), width = 0.5, linewidth = 1.25) +
  labs(x = expression("10-fold CV R"^2), y = NULL, title = NULL) +
  theme_minimal() +
  theme(
    # axis.text.y = element_text(size = 10),
    # axis.text.x = element_text(angle = 0, hjust = 1), 
    legend.position = "none",
    title = element_text(size = 7)
  ) + 
  coord_cartesian(xlim = c(min_r2, max_r2)) + 
  scale_x_continuous(breaks = seq(min_r2, max_r2, 0.1))

ggsave("./graphing/SupA/graph_model_comparison.pdf", p, width = 3, height = 3, units = "in", bg = 'white')
# print(p)


# biochem_data <- read.csv("./results/model_compare_Biochem/summary.csv", 
#                  header = TRUE, sep = ","
# )
# # from model name remove everything after -
# biochem_data$model <- gsub("-.*", "", biochem_data$model)



# biochem_data$highlight <- ifelse(
#     biochem_data$model == "lgbm",
#     "BlueHighlight",
#     "NoHighlight"
#   )

# biochem_p <- ggplot(biochem_data, aes(x = nan_r2, y = reorder(model, nan_r2, FUN = mean), fill = highlight)) +
#   geom_bar(stat = "summary", fun = "mean") +
#   scale_fill_manual(
#     name = "Model",
#     values = c("RedHighlight" = "red", "BlueHighlight" = "#0091ff", "NoHighlight" = "grey50")
#   ) +
#   geom_point(aes(x = nan_r2), color = "black", size = 1) +
#   labs(x = "R2 Score", y = NULL, title = "Model Comparison on L_P5_P3_C_3mer5_Struct_Biochem Features") +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_text(size = 10),
#     # axis.text.x = element_text(angle = 0, hjust = 1), 
#     legend.position = "none",
#     title = element_text(size = 7)
#   ) +
#   coord_cartesian(xlim = c(min(biochem_data$nan_r2), max(biochem_data$nan_r2)))

# ggsave("./graphing/graph_model_comparison_Biochem.pdf", biochem_p, width = 5, height = 5, units = "in", bg = 'white')
# # print(p)

# biochem_data$highlight <- "Biochem"

# data$highlight <- "no Biochem"

# combined_data <- data.frame(
#   model = c(data$model, biochem_data$model),
#   nan_r2 = c(data$nan_r2, biochem_data$nan_r2),
#   highlight = c(data$highlight, biochem_data$highlight)
# )


# combined_plot <- ggplot(combined_data, aes(x = model, y = nan_r2, fill = highlight)) +
#   geom_bar(stat = "summary", fun = "mean", position = position_dodge()) +
#   scale_fill_manual(
#     name = "Model",
#     values = c("Biochem" = "red", "no Biochem" = "#0091ff")
#   ) +
#   geom_point(aes(x = model, y = nan_r2), color = "black", size = 1, position = position_dodge(0.9)) +
#   labs(x = "R2 Score", y = NULL, title = "L_P5_P3_C_3mer5_Struct w/wo Biochem Features") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(size = 10, angle = 45),
#     # axis.text.x = element_text(angle = 0, hjust = 1), 
#     # legend.position = "none",
#     title = element_text(size = 7)
#   ) +
#   coord_cartesian(ylim = c(min(combined_data$nan_r2), max(combined_data$nan_r2)))

# ggsave("./graphing/graph_model_comparison_combined.pdf", combined_plot, width = 5, height = 5, units = "in", bg = 'white')