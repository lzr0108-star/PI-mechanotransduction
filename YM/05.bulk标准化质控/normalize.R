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
# 05 标准化 —— 只保留 Healthy + Periimplantitis（去除 Periodontitis）
library(edgeR); library(limma)

setwd(file.path(base, "05.bulk标准化质控"))

cnt0 <- as.matrix(read.csv(file.path(base,"03.bulk数据下载/GSE223924_raw_count.csv"),
                           row.names=1, check.names=FALSE))
prefix0 <- sub("[0-9]+$", "", colnames(cnt0))
keep.s <- prefix0 %in% c("H","PI")                 # ★ 去掉 PT(牙周炎)
cnt <- cnt0[, keep.s]
prefix <- prefix0[keep.s]
num <- sub("^[A-Z]+", "", colnames(cnt))
gname <- c(H="Healthy", PI="Periimplantitis")[prefix]
group <- factor(gname, levels=c("Healthy","Periimplantitis"))

clinical <- data.frame(sample=colnames(cnt), group=as.character(group),
                       patient=num, disease=ifelse(prefix=="H",0,1))
write.csv(clinical, "clinical.csv", row.names=FALSE)

y <- DGEList(counts=cnt, group=group)
design <- model.matrix(~0+group); colnames(design)<-levels(group)
keep <- filterByExpr(y, design)
cat("过滤前:", nrow(y$counts), "\n")
y <- y[keep,,keep.lib.sizes=FALSE]
cat("过滤后:", nrow(y$counts), "\n")
y <- calcNormFactors(y, method="TMM")

logcpm <- cpm(y, log=TRUE, prior.count=2)
out <- data.frame(geneNames=rownames(logcpm), logcpm, check.names=FALSE)
write.table(out, "expr.logCPM.txt", sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(y, "dge.rds")
cat("样本:", ncol(cnt), "(Healthy", sum(prefix=="H"), "+ PI", sum(prefix=="PI"), ")  基因:", nrow(logcpm), "\n")
cat("已输出 clinical.csv / expr.logCPM.txt / dge.rds\n")
