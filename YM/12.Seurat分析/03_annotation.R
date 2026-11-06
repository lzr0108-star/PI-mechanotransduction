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
# 12-③ 细胞类型注释（marker 打分）
library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(cowplot)

setwd(file.path(root,"YM/12.Seurat分析"))
seu <- readRDS("seu_integrated.rds")

# 细胞类型 marker（口腔黏膜/种植体周围组织）
markers <- list(
  Epithelial = c("EPCAM","KRT5","KRT14","KRT13","KRT4","KRT17","TP63","SFN","KRT19"),
  T          = c("CD3D","CD3E","CD3G","CD2","TRAC","IL7R","CD4","CD8A","CD8B","LCK"),
  B          = c("MS4A1","CD79A","CD79B","CD19","CD74","BANK1","CD22"),
  Plasma     = c("MZB1","JCHAIN","SDC1","IGHG1","IGKC","XBP1","DERL3"),
  Myeloid    = c("LYZ","LST1","CD68","CSF1R","TYROBP","AIF1","C1QA","C1QB","C1QC","CD163"),
  Neutrophil = c("FCGR3B","CSF3R","S100A12","CXCR2","FUT4","S100A8","S100A9"),
  Mast       = c("TPSAB1","TPSB2","CPA3","MS4A2","HDC","KIT"),
  Endothelial= c("PECAM1","VWF","CDH5","CLDN5","RAMP2","ENG","KDR","EMCN"),
  Fibroblast = c("COL1A1","COL3A1","COL1A2","DCN","LUM","CXCL12","PDGFRA","COL6A1"),
  SmoothMuscle=c("ACTA2","MYH11","TAGLN","MYL9","RGS5","MCAM")
)

# cluster 平均表达（lognorm）
Idents(seu) <- "seurat_clusters"
avg <- AverageExpression(seu, assays="RNA", group.by="seurat_clusters")$RNA
# 仅用 marker 基因
mg <- unique(unlist(markers)); mg <- intersect(mg, rownames(avg))
avgM <- avg[mg, , drop=FALSE]
# 基因在 cluster 间 z-score
z <- t(scale(t(avgM)))
# 每类型打分
scoreTab <- sapply(names(markers), function(ct){
  g <- intersect(markers[[ct]], rownames(z))
  if(length(g)==0) return(rep(NA,ncol(z)))
  colMeans(z[g,,drop=FALSE])
})
scoreTab <- as.data.frame(scoreTab)
scoreTab$cluster <- rownames(scoreTab)
# 最高分类型
assign <- setNames(colnames(scoreTab)[max.col(as.matrix(scoreTab[,names(markers)]))],
                   scoreTab$cluster)
# 标记低置信（最高分<0.3）为 unknown
maxv <- apply(scoreTab[,names(markers)],1,max)
assign[maxv < 0.3] <- "Unknown"
print(data.frame(cluster=names(assign), celltype=assign, score=round(maxv,2)), row.names=FALSE)

# 写回对象
seu$celltype <- unname(assign[paste0("g", seu$seurat_clusters)])
Idents(seu) <- "celltype"
saveRDS(seu, "seu_annotated.rds")

# 注释 UMAP
p1 <- DimPlot(seu, group.by="celltype", label=TRUE) + NoLegend()
p2 <- DimPlot(seu, group.by="condition", cols=c(Healthy="#4DBBD5",Periimplantitis="#E64B35"))
pdf("UMAP_annotated.pdf", width=15, height=6)
print(plot_grid(p1,p2))
dev.off()

# 各细胞类型比例（按condition）
prop <- prop.table(table(seu$celltype, seu$condition), margin=2)
print(round(prop,3))
write.table(as.data.frame.matrix(round(prop,3)), "celltype_proportion.txt", sep="\t", quote=FALSE)
cat("\n已输出注释结果\n")
