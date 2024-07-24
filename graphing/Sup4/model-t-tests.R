# From: https://github.com/vagarwal87/saluki_paper/blob/main/Fig3_S4/Fig3_S4e.R
# Fig 3 of Saluki Paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02811-x#Sec2

a=read.csv("./results/model_compare/summary_pivot.csv", header=T, sep=",", row.names=1)
round(a,3)
# print(a)

results = data.frame(tests = c(
    "randomforest vs lasso",
    "lasso vs elasticnet",
    "elasticnet vs lgbm"
    ), 
    pvals = c(
        t.test(a$randomforest, a$lasso, paired=T ,alternative='less')$p.value,
        t.test(a$lasso, a$elasticnet, paired=T ,alternative='less')$p.value,
        t.test(a$elasticnet, a$lgbm, paired=T ,alternative='less')$p.value
    )
)

results$p.corrected = pmin(results$pvals*nrow(results),1)
results$significant = results$p.corrected < 0.05
write.csv(results, "./graphing/SupA/model-t-tests.csv", row.names=F, quote=F)
