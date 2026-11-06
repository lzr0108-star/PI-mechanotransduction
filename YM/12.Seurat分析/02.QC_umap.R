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
# 12 (2) QC 诊断 + 复用 CELLxGENE UMAP/注释
library(reticulate); library(Seurat); library(Matrix); library(cowplot); library(ggplot2)

setwd(file.path(base, "12.Seurat分析"))
srt <- readRDS("seurat_ref.rds")

# ---- QC 诊断（基于 count）----
pdf("QC.pdf", width=16, height=6)
p1 <- VlnPlot(srt, features=c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"),
              pt.size=0, group.by="disease") & theme_bw()
print(p1)
p2 <- FeatureScatter(srt, "nCount_RNA","nFeature_RNA", pt.size=0.3)
p3 <- FeatureScatter(srt, "nCount_RNA","percent.mt", pt.size=0.3)
print(plot_grid(p2,p3))
dev.off()
cat("QC summary:\n")
print(summary(srt$nFeature_RNA)); print(summary(srt$percent.mt))

# ---- 从 h5ad 取已算好的 X_umap ----
ad <- import("anndata", convert=FALSE)
a <- ad$read_h5ad(file.path(base,"10.单细胞数据下载/periodontal_atlas.h5ad"))
um <- py_to_r(a$obsm["X_umap"])
rownames(um) <- colnames(srt); colnames(um) <- c("umap_1","umap_2")
srt[["umap_ref"]] <- CreateDimReducObject(embeddings=um, key="umapref_")

# ---- UMAP 可视化 ----
pdf("UMAP.overview.pdf", width=18, height=6)
d1 <- DimPlot(srt, reduction="umap_ref", group.by="cell_type", label=TRUE, label.size=2.5, raster=TRUE) + NoLegend() + ggtitle("Cell type")
d2 <- DimPlot(srt, reduction="umap_ref", group.by="disease", raster=TRUE) + ggtitle("Disease")
d3 <- DimPlot(srt, reduction="umap_ref", group.by="Study", raster=TRUE) + ggtitle("Study")
print(d1+d2+d3)
dev.off()

# 分组比例
cat("\n各 disease 中细胞类型比例（主要类）:\n")
ct_keep <- names(sort(table(srt$cell_type), decreasing=TRUE)[1:12])
kl <- srt$cell_type %in% ct_keep
tab <- prop.table(table(ct=factor(srt$cell_type[kl], levels=ct_keep),
                       disease=srt$disease[kl]), margin=2)
print(round(tab,3))

saveRDS(srt, "seurat_ready.rds")
cat("\n已保存 seurat_ready.rds（分析就绪）\n")
