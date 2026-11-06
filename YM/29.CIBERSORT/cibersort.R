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
# 29-31 CIBERSORT 免疫浸润 + 组间差异 + 核心基因相关性
library(edgeR); library(limma); library(ggplot2); library(reshape2); library(pheatmap); library(cowplot)

ym <- file.path(root,"YM")
setwd(file.path(ym,"29.CIBERSORT"))

y <- readRDS(file.path(ym,"05.bulk标准化质控/dge.rds"))
clinical <- read.csv(file.path(ym,"05.bulk标准化质控/clinical.csv"))
core <- readLines(file.path(ym,"23.交集特征基因/core_genes.txt"))

# ---- 准备 CIBERSORT 输入：非 log 的 TMM-CPM（基因 x 样本）----
X <- cpm(y)                       # 非log，TMM
X <- X[, clinical$sample]
inputD <- data.frame(Gene=rownames(X), X, check.names=FALSE)
write.table(inputD, "cibersort_input.txt", sep="\t", quote=FALSE, row.names=FALSE)

# ---- 跑 CIBERSORT ----
ref <- file.path(root,"29.CIBERSORT/ref.txt")
source(file.path(root,"29.CIBERSORT/MRG29.CIBERSORT.R"))
set.seed(2026)
res <- CIBERSORT(ref, "cibersort_input.txt", perm=100)
write.table(data.frame(sample=rownames(res), res, check.names=FALSE),
            "CIBERSORT_raw.txt", sep="\t", quote=FALSE, row.names=FALSE)
cellNames <- colnames(res)[1:(ncol(res)-3)]
C <- as.data.frame(res[, cellNames])
C$sample <- rownames(C)
C$group <- clinical$group[match(C$sample, clinical$sample)]
write.table(C, "CIBERSORT_fraction.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("CIBERSORT 完成，样本:", nrow(C), " 细胞类型:", length(cellNames),"\n")

# ---- 30 组间免疫差异 ----
long <- melt(C[, c(cellNames,"group","sample")], id.vars=c("group","sample"),
             variable.name="cell", value.name="fraction")
# Wilcox 每细胞
pvals <- sapply(cellNames, function(cc){
  h <- C[C$group=="Healthy", cc]; p <- C[C$group=="Periimplantitis", cc]
  tryCatch(wilcox.test(h,p)$p.value, error=function(e) NA)
})
sigTab <- data.frame(cell=cellNames, p.value=pvals,
                     mean_H=sapply(cellNames, function(cc) mean(C[C$group=="Healthy",cc])),
                     mean_PI=sapply(cellNames, function(cc) mean(C[C$group=="Periimplantitis",cc])))
sigTab$sig <- ifelse(sigTab$p.value<.05,"*","")
write.table(sigTab, "CIBERSORT_diff.txt", sep="\t", quote=FALSE, row.names=FALSE)

# boxplot（22细胞，facet）
long$cell <- factor(long$cell, levels=cellNames)
p <- ggplot(long, aes(group, fraction, fill=group)) +
  geom_boxplot(outlier.size=.4, width=.7) +
  geom_jitter(width=.08, size=.7, alpha=.6) +
  scale_fill_manual(values=c("Healthy"="#4DBBD5","Periimplantitis"="#E64B35")) +
  facet_wrap(~cell, scales="free_y", ncol=5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
        panel.grid=element_blank(), legend.position="bottom")
ggsave("CIBERSORT_boxplot.pdf", p, width=13, height=10)
cat("显著差异免疫细胞:\n"); print(sigTab[sigTab$p.value<.05,c("cell","p.value","mean_H","mean_PI")], row.names=FALSE)

# ---- 31 核心基因-免疫细胞相关性（Spearman）----
logcpm <- as.matrix(read.table(file.path(ym,"05.bulk标准化质控/expr.logCPM.txt"),
                               header=TRUE, sep="\t", row.names=1, check.names=FALSE))
gmat <- t(logcpm[core, C$sample])    # samples x core
cmat <- as.matrix(C[, cellNames])    # samples x cells
CM <- matrix(NA, length(core), length(cellNames),
             dimnames=list(core, cellNames))
PM <- CM
for(i in seq_along(core)) for(j in seq_along(cellNames)){
  ct <- suppressWarnings(cor.test(gmat[,i], cmat[,j], method="spearman", exact=FALSE))
  CM[i,j] <- ct$estimate; PM[i,j] <- ct$p.value
}
write.table(data.frame(gene=rownames(CM), round(CM,3)),
            "core_immune_cor.txt", sep="\t", quote=FALSE, row.names=FALSE)
# 星号热图
star <- ifelse(PM<.05, ifelse(PM<.01,"**","*"),"")
pdf("core_immune_heatmap.pdf", width=13, height=4.5)
pheatmap(CM, display_numbers=star, number_color="black",
         color=colorRampPalette(c("#3C5488","white","#E64B35"))(50),
         border_color=NA, cluster_cols=FALSE, cluster_rows=FALSE,
         fontsize_row=10, fontsize_col=7)
dev.off()
cat("\n已输出 29-31 CIBERSORT 全部结果\n")
