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
# 32-34 ConsensusClusterPlus 分型 + PCA验证 + 亚型免疫差异
library(ConsensusClusterPlus); library(limma)
library(ggplot2); library(pheatmap); library(reshape2); library(cowplot)

ym <- file.path(root,"YM")
setwd(file.path(ym,"32.样品分型"))

clinical <- read.csv(file.path(ym,"05.bulk标准化质控/clinical.csv"))
lcpm <- as.matrix(read.table(file.path(ym,"05.bulk标准化质控/expr.logCPM.txt"),
                             header=TRUE, sep="\t", row.names=1, check.names=FALSE))
lcpm <- lcpm[, clinical$sample]                 # 基因 x 样本
C <- read.delim(file.path(ym,"29.CIBERSORT/CIBERSORT_fraction.txt"), check.names=FALSE)
cellNames <- setdiff(colnames(C), c("sample","group"))

# ---- 32 ConsensusClusterPlus：CCP 对「列(样本)」聚类，输入 基因(行) x 样本(列) ----
madg <- apply(lcpm,1,mad)
topG <- names(sort(madg, decreasing=TRUE))[1:min(5000,length(madg))]
dat <- lcpm[topG, , drop=FALSE]                # 5000基因 x 20样本
# 按基因(行) z-score
dat2 <- t(apply(dat, 1, function(z) (z-mean(z))/sd(z)))   # 仍是 基因 x 样本
dir.create("CCP_output", showWarnings=FALSE)
res <- ConsensusClusterPlus(dat2, maxK=6, reps=100, pItem=0.8, pFeature=1,
                            clusterAlg="km", distance="euclidean",
                            title=file.path(getwd(),"CCP_output"),
                            plot="pdf", seed=2026)
K <- 2
cl <- res[[K]]$consensusClass                  # 名字=样本（列）
cluster <- setNames(paste0("C",cl), names(cl))
cat("选定 K =",K," 亚型分布:\n"); print(table(cluster))
ct <- table(cluster, true=clinical$group[match(names(cluster),clinical$sample)])
print(ct)
write.table(data.frame(sample=names(cluster), cluster=as.character(cluster),
                       true=clinical$group[match(names(cluster),clinical$sample)]),
            "cluster_assignment.txt", sep="\t", quote=FALSE, row.names=FALSE)

# ---- 33 PCA 按 cluster（PCA：样本为列，转成行）----
pcaInput <- t(dat2)                            # 样本 x 基因
pca <- prcomp(pcaInput, scale.=FALSE)
ve <- round(100*pca$sdev^2/sum(pca$sdev^2),1)
pc <- as.data.frame(pca$x[,1:2]); colnames(pc)<-c("PC1","PC2")
pc$cluster <- cluster[rownames(pc)]
pc$true <- clinical$group[match(rownames(pc),clinical$sample)]
p <- ggplot(pc, aes(PC1,PC2, color=cluster, shape=true)) +
  geom_point(size=3) +
  stat_ellipse(aes(group=cluster), lty=2, level=.95) +
  ggrepel::geom_text_repel(aes(label=rownames(pc)), size=2.5, show.legend=FALSE, max.overlaps=20) +
  scale_color_manual(values=c(C1="#E64B35",C2="#4DBBD5",C3="#00A087")) +
  labs(x=paste0("PC1 (",ve[1],"%)"), y=paste0("PC2 (",ve[2],"%)")) +
  theme_bw() + theme(panel.grid=element_blank())
ggsave(file.path(ym,"33.分型PCA/PCA_cluster.pdf"), p, width=7.5, height=6)

# ---- 34 亚型间免疫差异 ----
C$cluster <- cluster[match(C$sample, names(cluster))]
C <- C[!is.na(C$cluster),]
long <- melt(C[, c(cellNames,"cluster","sample")], id.vars=c("cluster","sample"),
             variable.name="cell", value.name="fraction")
pvals <- sapply(cellNames, function(cc){
  a <- C[C$cluster=="C1",cc]; b <- C[C$cluster=="C2",cc]
  if(length(unique(c(a,b)))<2) return(NA)
  tryCatch(wilcox.test(a,b)$p.value, error=function(e) NA)
})
idt <- data.frame(cell=cellNames, p.value=pvals,
                  mean_C1=sapply(cellNames, function(cc) mean(C[C$cluster=="C1",cc])),
                  mean_C2=sapply(cellNames, function(cc) mean(C[C$cluster=="C2",cc])))
write.table(idt, file.path(ym,"34.分型免疫差异/cluster_immune_diff.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)
long$cell <- factor(long$cell, levels=cellNames)
pb <- ggplot(long, aes(cluster, fraction, fill=cluster)) +
  geom_boxplot(outlier.size=.4) + geom_jitter(width=.08,size=.7,alpha=.6) +
  scale_fill_manual(values=c(C1="#E64B35",C2="#4DBBD5")) +
  facet_wrap(~cell, scales="free_y", ncol=5) +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
                     panel.grid=element_blank(), legend.position="bottom")
ggsave(file.path(ym,"34.分型免疫差异/cluster_immune_boxplot.pdf"), pb, width=13, height=10)
cat("\n亚型间显著免疫细胞:\n"); print(idt[!is.na(idt$p.value)&idt$p.value<.05,], row.names=FALSE)
cat("\n已输出 32-34 分型结果\n")
