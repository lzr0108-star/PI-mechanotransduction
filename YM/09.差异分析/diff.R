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
# 09 差异分析（H+PI，limma voom，PI vs Healthy）
library(edgeR); library(limma); library(ggplot2); library(ggrepel); library(pheatmap)

setwd(file.path(base, "09.差异分析"))
y <- readRDS(file.path(base,"05.bulk标准化质控/dge.rds"))
clinical <- read.csv(file.path(base,"05.bulk标准化质控/clinical.csv"))

logFCcut <- 0.585; pcut <- 0.05
group <- factor(clinical$group, levels=c("Healthy","Periimplantitis"))
design <- model.matrix(~0+group); colnames(design)<-levels(group)
v <- voom(y, design, plot=FALSE)
fit <- lmFit(v, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(PI_vs_H="Periimplantitis-Healthy", levels=design)))

de <- topTable(fit2, coef="PI_vs_H", number=Inf, sort.by="P")
de$gene <- rownames(de)
de$sig <- ifelse(de$adj.P.Val<pcut & abs(de$logFC)>logFCcut,
                 ifelse(de$logFC>0,"Up","Down"),"NS")
write.table(de, "DE.PI_vs_H.all.txt", sep="\t", quote=FALSE, row.names=FALSE)
sig <- de[de$sig!="NS",]
write.table(sig, "DE.PI_vs_H.sig.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("PI vs Healthy: 上调 %d 下调 %d（|logFC|>%.3f, adj.P<%.2f）\n",
    sum(sig$sig=="Up"), sum(sig$sig=="Down"), logFCcut, pcut)

# 火山图
p <- ggplot(de, aes(logFC,-log10(adj.P.Val))) +
  geom_point(aes(color=sig), size=.8, alpha=.6) +
  scale_color_manual(values=c(Up="#E64B35",Down="#3C5488",NS="grey80")) +
  geom_vline(xintercept=c(-logFCcut,logFCcut), lty=2, col="grey50") +
  geom_hline(yintercept=-log10(pcut), lty=2, col="grey50") +
  geom_text_repel(data=subset(de,sig!="NS")[1:20,], aes(label=gene),
                  size=2.6, max.overlaps=25, segment.size=.3) +
  labs(title="Peri-implantitis vs Healthy", x="log2 Fold Change", y="-log10(adj.P)") +
  theme_bw() + theme(panel.grid=element_blank(), legend.position=c(.9,.9), legend.title=element_blank())
ggsave("volcano.pdf", p, width=7, height=6)

# 热图（top 上调/下调各50）
n <- 50
up <- head(rownames(subset(de,sig=="Up")), n)
dn <- head(rownames(subset(de,sig=="Down")), n)
mat <- v$E[c(up,dn),]
ann <- data.frame(group=group); rownames(ann)<-clinical$sample
pdf("heatmap.PI_vs_H.pdf", width=9, height=10)
pheatmap(mat, annotation_col=ann, color=colorRampPalette(c("navy","white","firebrick3"))(50),
         scale="row", cluster_cols=FALSE, fontsize_row=5, fontsize_col=8, border_color=NA)
dev.off()
cat("已输出 DE 表 / volcano.pdf / heatmap.PI_vs_H.pdf\n")
