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
# 13 CellChat：Healthy vs PI 细胞通讯比较
library(Seurat); library(CellChat); library(ggplot2)
library(patchwork); library(cowplot)

setwd(file.path(root,"YM/13.细胞通讯"))
seu <- readRDS(file.path(root,"YM/12.Seurat分析/seu_annotated.rds"))

# ---- 抽样：每 condition x celltype 最多 800 细胞 ----
cells <- c()
for(cond in c("Healthy","Periimplantitis")){
  for(ct in unique(seu$celltype)){
    id <- colnames(seu)[seu$condition==cond & seu$celltype==ct]
    if(length(id)>800) id <- sample(id,800)
    cells <- c(cells, id)
  }
}
seuSub <- subset(seu, cells=cells)
cat("抽样后细胞:", ncol(seuSub), "\n"); print(table(seuSub$condition))

# ---- 单组 CellChat 流程 ----
runCC <- function(s){
  s[["RNA"]] <- JoinLayers(s[["RNA"]])
  dataInput <- s[["RNA"]]$data
  meta <- data.frame(labels=as.character(s$celltype)); rownames(meta)<-colnames(s)
  cc <- createCellChat(object=dataInput, meta=meta, group.by="labels")
  cc@DB <- subsetDB(CellChatDB.human, search="Secreted Signaling")
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- projectData(cc, PPI.human)
  cc <- computeCommunProb(cc, raw.use=TRUE)
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}
ccH <- runCC(subset(seuSub, subset=condition=="Healthy"))
ccP <- runCC(subset(seuSub, subset=condition=="Periimplantitis"))
object.list <- list(Healthy=ccH, PI=ccP)
cc <- mergeCellChat(object.list, add.names=names(object.list))
saveRDS(cc, "cellchat_compare.rds")

# ---- 可视化 ----
# 1) 各组网络 count / weight
g1 <- netVisual_circle(ccH@net$count, title.name="Healthy: Number")
g2 <- netVisual_circle(ccP@net$count, title.name="PI: Number")
pdf("01_network_count.pdf", width=14, height=7)
print(netVisual_circle(ccH@net$count, title.name="Healthy: interactions"))
print(netVisual_circle(ccP@net$count, title.name="PI: interactions"))
dev.off()

# 2) 比较总交互数与强度
p1 <- compareInteractions(cc, show.legend=FALSE, group=c(1,2))
p2 <- compareInteractions(cc, show.legend=FALSE, group=c(1,2), measure="weight")
pdf("02_compare_overall.pdf", width=10, height=5)
print(p1 + p2)
dev.off()

# 3) 差异交互网络（PI vs Healthy）
pdf("03_diff_interaction.pdf", width=12, height=6)
print(netVisual_diffInteraction(cc, weight=FALSE, title.name="Differential #interactions"))
print(netVisual_diffInteraction(cc, weight=TRUE, title.name="Differential strength"))
dev.off()

# 4) 通路水平 heatmap
pdf("04_pathway_heatmap.pdf", width=11, height=9)
print(netVisual_heatmap(cc, measure="weight", color.heatmap="Reds", title.name=""))
dev.off()

# 5) bubble（两组）
pdf("05_bubble.pdf", width=11, height=8)
print(netVisual_bubble(cc, comparison=c(1,2), angle.x=45, font.size=8, thresh=0.05))
dev.off()

# 汇总
cat("\n各组总交互数: Healthy", sum(ccH@net$count), " PI", sum(ccP@net), "\n")
cat("已输出 CellChat 结果\n")
