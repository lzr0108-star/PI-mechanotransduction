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
# 35 核心基因展示
library(Seurat); library(ggplot2); library(patchwork)

setwd(file.path(root,"YM/35.基因展示"))
seu <- readRDS(file.path(root,"YM/12.Seurat分析/seu_annotated.rds"))
Idents(seu) <- "celltype"
hub <- c("ICAM1","HPN","COL19A1")
hub <- hub[hub %in% rownames(seu)]

# FeaturePlot
pdf("Scatter.pdf", width=12, height=5)
print(FeaturePlot(seu, features=hub, ncol=3, order=TRUE))
dev.off()

# VlnPlot
pdf("Violin.pdf", width=12, height=5)
print(VlnPlot(seu, features=hub, pt.size=0, ncol=3))
dev.off()

# DotPlot
pdf("Bubble.pdf", width=8, height=5)
print(DotPlot(seu, features=hub) +
      coord_flip() +
      scale_color_gradient2(low="#3C5488", mid="white", high="#E64B35") +
      theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()
cat("35 基因展示完成\n")
