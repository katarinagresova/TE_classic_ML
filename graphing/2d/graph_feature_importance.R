library(tidyverse)
library(cowplot)


codon_to_aa_dict <- c(
    "TCA" = "S",    # Serine
    "TCC" = "S",    # Serine
    "TCG" = "S",    # Serine
    "TCT" = "S",    # Serine
    "TTC" = "F",    # Phenylalanine
    "TTT" = "F",    # Phenylalanine
    "TTA" = "L",    # Leucine
    "TTG" = "L",    # Leucine
    "TAC" = "Y",    # Tyrosine
    "TAT" = "Y",    # Tyrosine
    "TAA" = "X",    # Stop
    "TAG" = "X",    # Stop
    "TGC" = "C",    # Cysteine
    "TGT" = "C",    # Cysteine
    "TGA" = "X",    # Stop
    "TGG" = "W",    # Tryptophan
    "CTA" = "L",    # Leucine
    "CTC" = "L",    # Leucine
    "CTG" = "L",    # Leucine
    "CTT" = "L",    # Leucine
    "CCA" = "P",    # Proline
    "CCC" = "P",    # Proline
    "CCG" = "P",    # Proline
    "CCT" = "P",    # Proline
    "CAC" = "H",    # Histidine
    "CAT" = "H",    # Histidine
    "CAA" = "Q",    # Glutamine
    "CAG" = "Q",    # Glutamine
    "CGA" = "R",    # Arginine
    "CGC" = "R",    # Arginine
    "CGG" = "R",    # Arginine
    "CGT" = "R",    # Arginine
    "ATA" = "I",    # Isoleucine
    "ATC" = "I",    # Isoleucine
    "ATT" = "I",    # Isoleucine
    "ATG" = "M",    # Methionine (start)
    "ACA" = "T",    # Threonine
    "ACC" = "T",    # Threonine
    "ACG" = "T",    # Threonine
    "ACT" = "T",    # Threonine
    "AAC" = "N",    # Asparagine
    "AAT" = "N",    # Asparagine
    "AAA" = "K",    # Lysine
    "AAG" = "K",    # Lysine
    "AGC" = "S",    # Serine
    "AGT" = "S",    # Serine
    "AGA" = "R",    # Arginine
    "AGG" = "R",    # Arginine
    "GTA" = "V",    # Valine
    "GTC" = "V",    # Valine
    "GTG" = "V",    # Valine
    "GTT" = "V",    # Valine
    "GCA" = "A",    # Alanine
    "GCC" = "A",    # Alanine
    "GCG" = "A",    # Alanine
    "GCT" = "A",    # Alanine
    "GAC" = "D",    # Aspartic Acid
    "GAT" = "D",    # Aspartic Acid
    "GAA" = "E",    # Glutamic Acid
    "GAG" = "E",    # Glutamic Acid
    "GGA" = "G",    # Glycine
    "GGC" = "G",    # Glycine
    "GGG" = "G",    # Glycine
    "GGT" = "G"     # Glycine
)
amino_acid_to_color <- c(
    "S" = "#FF0000",    # Serine
    "F" = "#ffbb00",    # Phenylalanine
    "L" = "#ffd500",    # Leucine
    "Y" = "#c3ff00",    # Tyrosine
    "X" = "#6abc00",    # Stop
    "C" = "#00ffa2",    # Cysteine
    "W" = "#00b3ff",    # Tryptophan
    "P" = "#005eff",    # Proline
    "H" = "#a86eff",    # Histidine
    "Q" = "#9900ff",    # Glutamine
    "R" = "#fa93ff",    # Arginine
    "I" = "#ff00e6",    # Isoleucine
    "M" = "#9d004f",    # Methionine (start)
    "T" = "#950000",    # Threonine
    "N" = "#006f82",    # Asparagine
    "K" = "#004123",    # Lysine
    "V" = "#0b005c",    # Valine
    "A" = "#5f2800",    # Alanine
    "D" = "#ce800a",    # Aspartic Acid
    "E" = "#00ffe5",    # Glutamic Acid
    "G" = "#744c4c",     # Glycine
    "none" = "gray75"
)

# graph feature importances
data <- read.csv(
    "./results/human/all_cell_lines/lgbm-LL_P5_P3_CF_AAF_3mer_freq_5/feature_importance.csv",
    header = TRUE,
    row.names = 1
)

features <- read.csv(
    "./graphing/D/features-LL_P5_P3_CF_AAF_3mer_freq_5_Struct.csv",
    header = TRUE
)
features <- features[, -grep("fold", names(features))]
features <- features[, -grep("SYMBOL", names(features))]
features <- features[, -grep("mean_te_true", names(features))]
features <- features[, -grep("mean_across_cell_lines_true", names(features))]

bio_source_cols <- grep("_true", names(features), value = TRUE)

# print("bio_source_cols length")
# print(length(bio_source_cols))

# replace infs in cols struct_min_dG_UTR5 and struct_min_dG_CDS with NA
features$struct_min_dG_UTR5[features$struct_min_dG_UTR5 == Inf] <- NA
features$struct_min_dG_CDS[features$struct_min_dG_CDS == Inf] <- NA

# compute mean correlation of biosources with features and add to data
# iterate over each feature
for (feature in row.names(data)) {
    # get the correlation of the feature with each biosource
    temp <- grep(feature, names(features), value = TRUE)
    corrs <- cor(features[, temp], features[, bio_source_cols], use = "complete.obs", method = "spearman")
    stopifnot(length(corrs) == length(bio_source_cols))
    stopifnot(!any(is.na(corrs)))
    # get the mean correlation
    mean_corr <- mean(corrs)
    # add the mean correlation to the data
    data[feature, "mean_corr"] <- mean_corr
}

data$amino_acid <- "none"
# append the Amino Acid code to the end of the row name if it contains a codon
for (i in 1:length(row.names(data))) {
    # if the col name contains a codon
    if (grepl("_[ATGC]{3}", row.names(data)[i])) {
        # get the amino acid code
        # codon get using regex
        codon <- regmatches(row.names(data)[i], regexpr("[ATGC]{3}", row.names(data)[i]))
        aa <- codon_to_aa_dict[codon]
        # append the amino acid code to the end of the col name
        row.names(data)[i] <- paste(row.names(data)[i], aa, sep = "_")
        # set the highlight color
        data$amino_acid[i] <- aa
    }
}

# x = data rownames, y is mean_importance, plot with bar chart
# only show top 50 features
data <- data[order(data$mean_importance, decreasing = TRUE), ]
data$log_mean_importance <- log10(data$mean_importance)

first <- data[1:40, ]
# second <- data[51:100, ]


p <- ggplot(first, aes(x = log_mean_importance, y = reorder(row.names(first), log_mean_importance))) +
    geom_bar(stat = "identity", aes(fill = mean_corr)) +
    scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "Mean spearman corr.") +
    labs(x = expression("log"[10]*" of mean importance"), y = "Top 40 features") +
    theme(
        axis.text.y = element_text(size = 5),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        legend.position = "none",
    )
    # coord_cartesian(xlim = c(0, 6000))
    # scale_x_continuous(trans = "pseudo_log")

# p2 <- ggplot(second, aes(x = mean_importance, y = reorder(row.names(second), mean_importance), fill = amino_acid)) +
#     geom_bar(stat = "identity") +
#     scale_fill_manual(
#         name = "Amino Acid",
#         values = amino_acid_to_color,
#     ) +
#     labs(x = "Mean importance", y = "Feature") +
#     theme(
#         axis.text.y = element_text(size = 6),
#     ) +
#     coord_cartesian(xlim = c(0, 6000)) +
#     # scale_x_continuous(trans = "pseudo_log")

# combined <- plot_grid(p, p2, ncol = 2)

legend <- get_legend(p + theme(
    legend.position = "bottom",
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 7),
))
p <- plot_grid(p, legend, ncol = 1, nrow = 2, rel_heights = c(1, 0.1))

# save p
ggsave(
    # filename = "./graphing/graph_feature_importance.pdf",
    filename = "./graphing/D/graph_feature_importance.pdf",
    plot = p,
    width = 3,
    height = 5,
    units = "in"
)
