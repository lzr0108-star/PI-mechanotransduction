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
# 14 scRNAtoolVis 可视化
library(Seurat); library(dplyr); library(ggplot2)
library(ggunchull); library(scRNAtoolVis); library(ggsci)

setwd(file.path(root,"YM/14.scRNAtoolVis"))
seu <- readRDS(file.path(root,"YM/12.Seurat分析/seu_annotated.rds"))
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
Idents(seu) <- "celltype"

# marker（若已保存则直接读取）
if(file.exists("allGroupMarkers.txt")){
  mk <- read.delim("allGroupMarkers.txt", check.names=FALSE)
} else {
  mk <- FindAllMarkers(seu, only.pos=TRUE, min.pct=.25, logfc.threshold=.25)
  write.table(mk, "allGroupMarkers.txt", sep="\t", quote=FALSE, row.names=FALSE)
}
cat("marker 数:", nrow(mk), "\n")
tilecol <- jjAnno::useMyCol("paired", n=12)

# 01 细胞类型图
pdf("01_cellType.pdf", width=9, height=7)
print(clusterCornerAxes(seu, reduction="umap", clusterCol="celltype",
                       noSplit=TRUE, cellLabel=TRUE, addCircle=TRUE,
                       cicDelta=.2, cicAlpha=.15, cicLineSize=1,
                       cellLabelSize=4, show.legend=FALSE))
dev.off()

# 02 比例柱状图
pdf("02_barplot.pdf", width=7, height=6)
print(cellRatioPlot(seu, sample.name="condition", celltype.name="celltype", col.width=.7))
dev.off()

# 03 markerVolcano
pdf("03_markerVolcano.pdf", width=11, height=6)
print(markerVolcano(markers=mk, topn=5, labelCol=pal_npg()(10), log2FC=1))
dev.off()

# 04 jjVolcano
pdf("04_jjVolcano1.pdf", width=12, height=6)
print(jjVolcano(diffData=mk, legend.position=c(.08,.95), tile.col=tilecol,
                log2FC.cutoff=1, adjustP.cutoff=.05) + labs(x="Cell type"))
dev.off()
pdf("04_jjVolcano2.pdf", width=9, height=7)
print(jjVolcano(diffData=mk, legend.position=c(.05,.95), flip=TRUE, tile.col=tilecol,
                log2FC.cutoff=1, adjustP.cutoff=.05) + labs(x="Cell type"))
dev.off()
pdf("04_jjVolcano3.pdf", width=8, height=7)
print(jjVolcano(diffData=mk, legend.position=c(.1,.85), polar=TRUE, size=3, tile.col=tilecol,
                log2FC.cutoff=1, adjustP.cutoff=.05) + labs(x="Cell type"))
dev.off()
cat("14 scRNAtoolVis 完成\n")
