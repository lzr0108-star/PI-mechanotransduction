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
# 36 细胞轨迹（slingshot）
library(Seurat); library(slingshot); library(ggplot2)
library(dplyr); library(patchwork); library(RColorBrewer)

setwd(file.path(root,"YM/36.细胞轨迹"))
seu <- readRDS(file.path(root,"YM/12.Seurat分析/seu_annotated.rds"))
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
Idents(seu) <- "celltype"

# ---- slingshot（UMAP 坐标 + celltype 标签）----
umap <- Embeddings(seu, "umap")
colnames(umap) <- c("UMAP_1","UMAP_2")
set.seed(1)
sds <- slingshot(umap, clusterLabels=seu$celltype)
pst <- slingPseudotime(sds)
cat("lineage 数:", ncol(pst), "\n")
print(colnames(pst))

# 取每个细胞的最小pseudotime（跨lineage）
pst_min <- apply(pst, 1, function(x) if(all(is.na(x))) NA else min(x, na.rm=TRUE))
seu$pseudotime <- pst_min

# 曲线坐标转 data.frame
curves <- slingCurves(sds)
curvedf <- do.call(rbind, lapply(seq_along(curves), function(i){
  data.frame(curves[[i]]$s[curves[[i]]$ord,], lineage=paste0("L",i))
}))
colnames(curvedf)[1:2] <- c("UMAP_1","UMAP_2")

# 细胞坐标
cdf <- data.frame(umap, celltype=seu$celltype, pseudotime=pst_min)

# 01 轨迹（按celltype + 曲线）
p1 <- ggplot(cdf, aes(UMAP_1,UMAP_2)) +
  geom_point(aes(color=celltype), size=.4, alpha=.7) +
  geom_path(data=curvedf, aes(group=lineage), color="black", linewidth=.8) +
  theme_bw() + labs(title="Trajectory by cell type")
ggsave("Trajectory01_cellType.pdf", p1, width=7.5, height=6.5)

# 02 轨迹（pseudotime 着色）
p2 <- ggplot(cdf, aes(UMAP_1,UMAP_2)) +
  geom_point(aes(color=pseudotime), size=.4, na.rm=TRUE) +
  geom_path(data=curvedf, aes(group=lineage), color="black", linewidth=.8) +
  scale_color_viridis_c(option="magma") + theme_bw() +
  labs(title="Pseudotime", color="Pseudotime")
ggsave("Trajectory02_pseudotime.pdf", p2, width=7.5, height=6.5)

# 03 pseudotime 按 celltype（小提琴）
p3 <- ggplot(cdf, aes(x=celltype, y=pseudotime, fill=celltype)) +
  geom_violin(na.rm=TRUE) +
  geom_boxplot(width=.15, fill="white", outlier.shape=NA, na.rm=TRUE) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none") +
  labs(x="Cell type", y="Pseudotime")
ggsave("Trajectory03_pseudotimeViolin.pdf", p3, width=9, height=5.5)

# 04 核心基因 pseudotime 热图
hub <- c("ICAM1","HPN","COL19A1")
hub <- hub[hub %in% rownames(seu)]
expr <- as.matrix(seu[["RNA"]]$data[hub, ])
# 按 pseudotime 排序细胞
o <- order(pst_min, na.last=NA)
exprO <- expr[, o]
library(pheatmap)
pdf("Trajectory04_geneHeatmap.pdf", width=10, height=3.5)
pheatmap(exprO, cluster_cols=FALSE, cluster_rows=FALSE,
         show_colnames=FALSE,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         annotation_legend=TRUE,
         main="Core genes along pseudotime")
dev.off()

saveRDS(seu, "seu_trajectory.rds")
cat("36 轨迹分析完成\n")
