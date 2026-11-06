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
# 12-② Harmony 整合 + 聚类 + UMAP
library(Seurat); library(harmony); library(ggplot2); library(cowplot)

setwd(file.path(root,"YM/12.Seurat分析"))
seuList <- readRDS("seuList_qc.rds")

# 每个样本标准化 + 高可变基因
seuList <- lapply(seuList, function(x){
  x <- NormalizeData(x, normalization.method="LogNormalize", verbose=FALSE)
  x <- FindVariableFeatures(x, nfeatures=2000, verbose=FALSE)
  x
})
hvf <- SelectIntegrationFeatures(object.list=seuList, nfeatures=2000)
cat("整合特征数:", length(hvf), "\n")

# merge（细胞名已唯一）
seu <- merge(seuList[[1]], seuList[-1])
VariableFeatures(seu) <- hvf
seu <- ScaleData(seu, features=hvf, verbose=FALSE)
seu <- RunPCA(seu, npcs=50, verbose=FALSE)
seu <- RunHarmony(seu, group.by.vars="sample", reduction.use="pca", verbose=FALSE)

seu <- FindNeighbors(seu, reduction="harmony", dims=1:30, verbose=FALSE)
seu <- FindClusters(seu, resolution=0.6, verbose=FALSE)
seu <- RunUMAP(seu, reduction="harmony", dims=1:30, verbose=FALSE)

saveRDS(seu, "seu_integrated.rds")
cat("聚类数:", length(levels(seu)), " 细胞:", ncol(seu), "\n")

# 可视化
p1 <- DimPlot(seu, reduction="umap", group.by="condition", cols=c(Healthy="#4DBBD5",Periimplantitis="#E64B35"))
p2 <- DimPlot(seu, reduction="umap", group.by="seurat_clusters", label=TRUE) + NoLegend()
p3 <- DimPlot(seu, reduction="umap", group.by="sample")
pdf("UMAP_integrated.pdf", width=18, height=6)
print(plot_grid(p1,p2,p3, ncol=3))
dev.off()
# Elbow
pdf("PCA_elbow.pdf", width=7, height=5)
print(ElbowPlot(seu, ndims=50))
dev.off()
cat("已输出 UMAP_integrated.pdf\n")
