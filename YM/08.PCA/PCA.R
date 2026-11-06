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
# 08 PCA（H+PI，二组着色 + 椭圆 + 标签）
library(ggplot2); library(ggrepel)

setwd(file.path(base, "08.PCA"))
logcpm <- as.matrix(read.table(file.path(base,"05.bulk标准化质控/expr.logCPM.txt"),
                               header=TRUE, sep="\t", row.names=1, check.names=FALSE))
clinical <- read.csv(file.path(base,"05.bulk标准化质控/clinical.csv"))

pca <- prcomp(t(logcpm), scale.=FALSE)
ve <- round(100*pca$sdev^2/sum(pca$sdev^2), 1)
pc <- as.data.frame(pca$x[,1:2]); colnames(pc)<-c("PC1","PC2")
pc$group <- clinical$group[match(rownames(pc), clinical$sample)]
pc$sample <- rownames(pc)

p <- ggplot(pc, aes(PC1,PC2, color=group, shape=group)) +
  geom_point(size=2.8) +
  stat_ellipse(level=0.95, linetype=2, linewidth=0.8, show.legend=FALSE) +
  geom_text_repel(aes(label=sample), size=2.8, show.legend=FALSE, max.overlaps=20, box.padding=0.4) +
  scale_color_manual(values=c("Healthy"="#4DBBD5","Periimplantitis"="#E64B35")) +
  labs(x=paste0("PC1 (",ve[1],"%)"), y=paste0("PC2 (",ve[2],"%)")) +
  theme_bw() + theme(panel.grid=element_blank())
pdf("PCA.pdf", width=7, height=6)
print(p)
dev.off()
cat("PCA.pdf  PC1=",ve[1],"% PC2=",ve[2],"%\n")
