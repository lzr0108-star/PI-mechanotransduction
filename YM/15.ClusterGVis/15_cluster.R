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
# 15 ClusterGVis：marker 基因聚类热图 + GO 注释
library(tidyverse); library(clusterProfiler); library(org.Hs.eg.db)
library(ComplexHeatmap); library(ClusterGVis); library(ggsci); library(Seurat)
set.seed(123456)

setwd(file.path(root,"YM/15.ClusterGVis"))
seu <- readRDS(file.path(root,"YM/12.Seurat分析/seu_annotated.rds"))
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
Idents(seu) <- "celltype"

# marker（读取已保存）
mk <- read.delim(file.path(root,"YM/14.scRNAtoolVis/allGroupMarkers.txt"), check.names=FALSE)
top_markers <- mk %>% group_by(cluster) %>% top_n(20, avg_log2FC)

# AverageExpression by celltype
avg_exp <- AverageExpression(seu, assays="RNA", layer="data")$RNA
heatmat <- avg_exp[top_markers$gene, , drop=FALSE]
heatmat <- heatmat[!duplicated(rownames(heatmat)), ]
heatmat <- as.matrix(heatmat)
cat("heatmat:", nrow(heatmat), "genes x", ncol(heatmat), "cell types\n")

# clusterData（kmeans，聚类数 = 细胞类型数）
ck <- clusterData(obj=heatmat, clusterMethod="kmeans",
                  clusterNum=length(unique(seu$celltype)))
cat("cluster.list 数量:", length(ck$cluster.list), "\n")

# GO 富集
all_go <- data.frame()
for(clu_name in names(ck$cluster.list)){
  genes <- ck$cluster.list[[clu_name]]
  if(length(genes)<5) next
  idmap <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  if(nrow(idmap)==0) next
  go <- enrichGO(idmap$ENTREZID, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                 ont="BP", pvalueCutoff=1)
  if(nrow(as.data.frame(go))==0) next
  go_df <- as.data.frame(go)[1:5, ]
  go_df$cluster <- clu_name
  all_go <- rbind(all_go, go_df)
}
go_plot <- all_go[, c("cluster","Description")]
write.table(go_plot, "GO_terms.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("GO 条目:", nrow(go_plot), "\n")

# visCluster
top_mk10 <- mk %>% group_by(cluster) %>% top_n(10, avg_log2FC)
pdf("heatmap.pdf", width=16, height=12, onefile=FALSE)
print(visCluster(object=ck, column_names_rot=45, lineSide="left",
                 markGenes=top_mk10$gene, markGenesSide="left",
                 plotType="both", annoTermData=go_plot,
                 annoTermMside="right", annnoblockText=FALSE))
dev.off()
cat("15 ClusterGVis 完成\n")
