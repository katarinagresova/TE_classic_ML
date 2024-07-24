# From: https://github.com/vagarwal87/saluki_paper/blob/main/Fig3_S4/Fig3_S4e.R
# Fig 3 of Saluki Paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02811-x#Sec2

a=read.csv("./results/feature_set_comparison/summary_pivot.csv", header=T, sep=",", row.names=1)
round(a,3)
# print(a)

# remove the "on_NA_lgbm." prefix from the model names
colnames(a) <- gsub("on_NA_lgbm.", "", colnames(a))

results <- data.frame(tests = c(
    "3mer_freq_5 vs AAF",
    "AAF vs CF",
    "CF vs CF_3mer_freq_5",
    "CF_3mer_freq_5 vs LL_CF_3mer_freq_5",
    "LL_CF_3mer_freq_5 vs LL_P5_P3_CF_3mer_freq_5",
    "LL_CF_3mer_freq_5 vs LL_CF_AAF_3mer_freq_5",
    "LL_CF_AAF_3mer_freq_5 vs LL_P5_P3_CF_AAF_3mer_freq_5",
    "LL_P5_P3_CF_3mer_freq_5 vs LL_P5_P3_CF_AAF_3mer_freq_5",
    "LL_CF_3mer_freq_5 vs LL_P5_P3_CF_AAF_3mer_freq_5",
    "LL_P5_P3_CF_AAF_3mer_freq_5 vs LL_P5_P3_CF_AAF_3mer_freq_5_Struct",
    "LL_P5_P3_CF_AAF_3mer_freq_5 vs ALL_FEATURES",
    "LL_P5_P3_CF_AAF_3mer_freq_5_Struct vs ALL_FEATURES",
    "LL_P5_P3_CF_AAF_3mer_freq_5_Struct vs LL_P5_P3_CF_AAF_3mer_freq_5_Struct_Biochem",
    "Biochem vs LL_P5_P3_CF_AAF_3mer_freq_5_Struct_Biochem",
    "LL_P5_P3_CF_AAF_3mer_freq_5 vs LL_P5_P3_CF_AAF_DCF_3mer_freq_5"
    ),
    pvals = c(
        t.test(a$`3mer_freq_5`, a$AAF, paired=T,alternative='less')$p.value,
        t.test(a$AAF, a$CF, paired=T,alternative='less')$p.value,
        t.test(a$CF, a$CF_3mer_freq_5, paired=T,alternative='less')$p.value,
        t.test(a$CF_3mer_freq_5, a$LL_CF_3mer_freq_5, paired=T,alternative='less')$p.value,
        t.test(a$LL_CF_3mer_freq_5, a$LL_P5_P3_CF_3mer_freq_5, paired=T,alternative='less')$p.value,
        t.test(a$LL_CF_3mer_freq_5, a$LL_CF_AAF_3mer_freq_5, paired=T,alternative='less')$p.value,
        t.test(a$LL_CF_AAF_3mer_freq_5, a$LL_P5_P3_CF_AAF_3mer_freq_5, paired=T,alternative='less')$p.value,
        t.test(a$LL_P5_P3_CF_3mer_freq_5, a$LL_P5_P3_CF_AAF_3mer_freq_5, paired=T,alternative='less')$p.value,
        t.test(a$LL_CF_3mer_freq_5, a$LL_P5_P3_CF_AAF_3mer_freq_5, paired=T,alternative='less')$p.value,
        t.test(a$LL_P5_P3_CF_AAF_3mer_freq_5, a$LL_P5_P3_CF_AAF_3mer_freq_5_Struct, paired=T,alternative='less')$p.value,
        t.test(a$LL_P5_P3_CF_AAF_3mer_freq_5, a$ALL_FEATURES, paired=T,alternative='less')$p.value,
        t.test(a$LL_P5_P3_CF_AAF_3mer_freq_5_Struct, a$ALL_FEATURES, paired=T,alternative='less')$p.value,
        t.test(a$LL_P5_P3_CF_AAF_3mer_freq_5_Struct, a$LL_P5_P3_CF_AAF_3mer_freq_5_Struct_Biochem, paired=T,alternative='less')$p.value,
        t.test(a$Biochem, a$LL_P5_P3_CF_AAF_3mer_freq_5_Struct_Biochem, paired=T,alternative='less')$p.value,
        t.test(a$LL_P5_P3_CF_AAF_3mer_freq_5, a$LL_P5_P3_CF_AAF_DCF_3mer_freq_5, paired=T,alternative='less')$p.value
    )
)

results$p.corrected = pmin(results$pvals*nrow(results),1)
results$significant = results$p.corrected < 0.05
write.csv(results, "./graphing/A/t-tests.csv", row.names=F, quote=F)
