library(tidyverse)
library(LSD)
library(cowplot)

predictions <- read.csv("./results/mouse/all_cell_lines/lgbm-LL_P5_P3_CF_AAF_3mer_freq_5/predictions.csv", 
  header = TRUE, sep = ","
)


predictions <- predictions %>% 
  rename(ObservedMeanTE = mean_across_cell_lines_true, PredictedMeanTE = mean_across_cell_lines_pred)


pdf("./graphing/mouse/C/predicted_vs_observed_mean_TE_mouse.pdf", 5, 5) 
# 'predictions' below is a dataframe with at least two columns named PredictedMeanTE and ObservedMeanTE
predictions %>% 
  mutate(
    r = cor(PredictedMeanTE, ObservedMeanTE),
    rho = cor(PredictedMeanTE, ObservedMeanTE, method = "spearman"),
    ) %>%
  {
    heatscatter(
      .$PredictedMeanTE, .$ObservedMeanTE,
      xlim = c(-2, 2),
      xlab = "CV fitted model prediction",
      ylab = "Mean TE", bty='n', cex=0.3, las=1,
      main = str_glue("r = {.$r[1] %>% round(3)}, rho = {.$rho[1] %>% round(3)}")
    )
  }
dev.off()

pdf("./graphing/mouse/C/hippocampal_vs_observed_mean_TE_mouse.pdf", 5, 5) 
# 'predictions' below is a dataframe with at least two columns named PredictedMeanTE and ObservedMeanTE
predictions %>% 
  mutate(
    r = cor(bio_source_hippocampal_true, ObservedMeanTE, use = "complete.obs"),
    rho = cor(bio_source_hippocampal_true, ObservedMeanTE, method = "spearman", use = "complete.obs"),
    ) %>%
  {
    heatscatter(
      .$bio_source_hippocampal_true, .$ObservedMeanTE,
      xlim = c(-4, 4),
      xlab = "Hippocampal TE",
      ylab = "Mean TE", bty='n', cex=0.3, las=1,
      main = str_glue("r = {.$r[1] %>% round(3)}, rho = {.$rho[1] %>% round(3)}")
    )
  }
dev.off()