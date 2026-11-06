### Auto-added: locate repository root and load configuration --------
if (!exists("root", inherits = FALSE)) {
  .cfg_d <- normalizePath(getwd(), mustWork = FALSE)
  for (.cfg_i in seq_len(10)) {
    if (file.exists(file.path(.cfg_d, "config.R"))) {
      source(file.path(.cfg_d, "config.R")); break
    }
    .cfg_up <- dirname(.cfg_d)
    if (identical(.cfg_up, .cfg_d)) break
    .cfg_d <- .cfg_up
  }
  rm(.cfg_d, .cfg_i, .cfg_up)
}
# 24 ROC + 26 列线图（核心基因 ICAM1 HPN COL19A1）
library(pROC); library(ggplot2); library(rms)

lcpm0 <- as.matrix(read.table(file.path(ym,"05.bulk标准化质控/expr.logCPM.txt"),
                              header=TRUE, sep="\t", row.names=1, check.names=FALSE))
clinical <- read.csv(file.path(ym,"05.bulk标准化质控/clinical.csv"))
core <- readLines(file.path(ym,"23.交集特征基因/core_genes.txt"))
core <- intersect(core, rownames(lcpm0))
X <- t(lcpm0[core, clinical$sample])
df <- data.frame(clinical, X, check.names=FALSE)
df$PI <- as.integer(df$group=="Periimplantitis")
cat("核心基因:", paste(core,collapse=" "),"\n")

# ---------- 24 ROC ----------
resp <- df$group=="Periimplantitis"
roclist <- setNames(lapply(core, function(g){
  roc(resp, df[[g]], quiet=TRUE)
}), core)
aucs <- sapply(roclist, function(r) as.numeric(auc(r)))
# 三基因联合 logistic
glmfit <- glm(PI ~ ICAM1 + HPN + COL19A1, data=df, family=binomial)
prob <- predict(glmfit, type="response")
roclist$Combined <- roc(resp, prob, quiet=TRUE)
aucComb <- as.numeric(auc(roclist$Combined))
labs <- c(setNames(sprintf("%s AUC=%.3f", core, aucs), core),
          Combined=sprintf("Combined AUC=%.3f", aucComb))
cols <- c(ICAM1="#E64B35", HPN="#4DBBD5", COL19A1="#00A087", Combined="black")
p <- ggroc(roclist, legacy.axes=TRUE, linewidth=.8) +
  scale_color_manual(values=cols, labels=labs) +
  geom_abline(slope=1, intercept=0, lty=2, col="grey60") +
  labs(x="1 - Specificity", y="Sensitivity", color="") +
  theme_bw() + theme(legend.position=c(.72,.25))
ggsave(file.path(ym,"24.ROC/ROC.pdf"), p, width=7, height=6)
cat("AUC: "); print(round(c(aucs, Combined=aucComb),3))

# ---------- 26 列线图 ----------
dd <- datadist(df); options(datadist="dd")
lfit <- lrm(PI ~ ICAM1 + HPN + COL19A1, data=df, x=TRUE, y=TRUE, penalty=1)
pdf(file.path(ym,"26.列线图/nomogram.pdf"), width=9, height=7)
nom <- nomogram(lfit, fun=plogis, fun.at=c(.05,.1,.2,.5,.8,.9,.95),
                 funlabel="Risk of Peri-implantitis", lp=FALSE)
plot(nom)
dev.off()
# 校准曲线
pdf(file.path(ym,"26.列线图/calibrate.pdf"), width=6.5, height=6)
cal <- calibrate(lfit, method="boot", B=1000,
                 title="Calibration (bootstrap 1000)")
plot(cal, subtitles=FALSE, cex=.9)
dev.off()
cat("已输出 ROC.pdf / nomogram.pdf / calibrate.pdf\n")
