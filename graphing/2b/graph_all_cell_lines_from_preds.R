library(ggplot2)
library(ggpubr)
library(cowplot)


preds <- read.csv("./results/human/all_cell_lines/lgbm-LL_P5_P3_CF_AAF_3mer_freq_5/predictions.csv", 
  header = TRUE, sep = ","
)

colnames(preds) <- gsub("bio_source_", "", colnames(preds))

bio_sources <- colnames(preds)
bio_sources <- bio_sources[!bio_sources  == "SYMBOL"]
bio_sources <- bio_sources[!bio_sources  == "fold"]


bio_sources <- gsub("_pred", "", bio_sources)
bio_sources <- gsub("_true", "", bio_sources)
bio_sources <- unique(bio_sources)

pred_cols <- paste0(bio_sources, "_pred")
pred_cols <- pred_cols[!pred_cols %in% c("mean_te_pred", "mean_across_cell_lines_pred")]

# add new col to preds with mean of all bio_sources_no_mean preds called mean_across_cell_lines_pred
preds$mean_across_cell_lines_pred <- rowMeans(preds[,pred_cols], na.rm = TRUE)
preds$mean_across_cell_lines_true <- preds$mean_te_true

bio_sources <- c(bio_sources, "mean_across_cell_lines")


# 2d histogram of mean_accross_cell_lines_pred vs mean_accross_cell_lines_true
hist_data <- preds[,c("mean_across_cell_lines_pred", "mean_across_cell_lines_true")]
hist_plot1 <- ggplot(hist_data, aes(x = mean_across_cell_lines_true, y = mean_across_cell_lines_pred)) +
    geom_bin2d(bins = 100) +
    labs(x = "Mean TE", y = "Predicted Mean", title = NULL) +
    scale_fill_continuous(type = "viridis", limits = c(0,50)) +
    # stat_cor(method = "pearson", label.x = -2, label.y = 2, label.sep = "\n", size = 3) +
    stat_cor(method = "spearman", cor.coef.name = "rho", label.x = -2, label.y = 2, label.sep = "\n", size = 3) +
    theme_minimal() +
    coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5)) + 
    theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        # axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position = "none",
        title = element_text(size = 12)
    )

hist_data <- preds[,c("HEK293T_true", "mean_across_cell_lines_true")]
hist_plot2 <- ggplot(hist_data, aes(x = mean_across_cell_lines_true, y = HEK293T_true)) +
    geom_bin2d(bins = 100) +
    labs(x = "Mean TE", y = "HEK293T TE", title = NULL) +
    scale_fill_continuous(type = "viridis", limits = c(0,50)) +
    # stat_cor(method = "pearson", label.x = -2, label.y = 2, label.sep = "\n", size = 3) +
    stat_cor(method = "spearman", cor.coef.name = "rho", label.x = -2, label.y = 2, label.sep = "\n", size = 3) +
    theme_minimal() +
    coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5)) + 
    theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        # axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position = "none",
        title = element_text(size = 12)
    )

legend <- get_legend(hist_plot1 + theme(legend.position = "right"))
combo_hist <- plot_grid(plot_grid(hist_plot1, hist_plot2, ncol = 2), legend, ncol = 2, rel_widths = c(1, 0.1))
ggsave("./graphing/C/graph_pred_vs_true_mean_te_and_HEK293T.pdf", combo_hist, width = 10, height = 4.75, units = "in", dpi = 300, bg = 'white')