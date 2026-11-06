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
# 38 高低表达组通讯（PI 组内按 ICAM1 中位分 High/Low）
library(Seurat); library(CellChat); library(ggplot2)
library(dplyr); library(patchwork); library(ggplotify)

setwd(file.path(root,"YM/38.高低表达通讯"))
seu <- readRDS(file.path(root,"YM/12.Seurat分析/seu_annotated.rds"))
seu <- subset(seu, subset=condition=="Periimplantitis")
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])

# 抽样：每 celltype 最多 600
cells <- c()
for(ct in unique(seu$celltype)){
  id <- colnames(seu)[seu$celltype==ct]
  if(length(id)>600) id <- sample(id,600)
  cells <- c(cells,id)
}
seuS <- subset(seu, cells=cells)
cat("抽样后 PI 细胞:", ncol(seuS), "\n")

# 按 ICAM1 中位分 High/Low
gene <- "ICAM1"
gex <- FetchData(seuS, vars=gene)[[1]]
med <- median(gex)
seuS$gene_group <- ifelse(gex>med, "High", "Low")
seuH <- subset(seuS, gene_group=="High")
seuL <- subset(seuS, gene_group=="Low")
cat("High:", ncol(seuH), " Low:", ncol(seuL), "\n")

runCC <- function(s){
  s[["RNA"]] <- JoinLayers(s[["RNA"]])
  dataInput <- s[["RNA"]]$data
  meta <- data.frame(labels=as.character(s$celltype)); rownames(meta)<-colnames(s)
  cc <- createCellChat(dataInput, meta=meta, group.by="labels")
  cc@DB <- CellChatDB.human
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
ccH <- runCC(seuH); ccL <- runCC(seuL)
mlist <- list(High=ccH, Low=ccL)
cc <- mergeCellChat(mlist, add.names=names(mlist))
saveRDS(cc, "cellchat_HL.rds")

# 01 总体数目/强度
p1 <- compareInteractions(cc, show.legend=FALSE, group=c(1,2))
p2 <- compareInteractions(cc, show.legend=FALSE, group=c(1,2), measure="weight")
ggsave("01_stat.pdf", p1+p2, width=10, height=5)

# 02 circle（count）
gsH <- as.vector(table(seuH$celltype)); gsL <- as.vector(table(seuL$celltype))
mc <- max(c(ccH@net$count, ccL@net$count))
c1 <- as.ggplot(~netVisual_circle(ccH@net$count, vertex.weight=gsH, weight.scale=TRUE,
                                 label.edge=FALSE, title.name="High: number", edge.weight.max=mc))
c2 <- as.ggplot(~netVisual_circle(ccL@net$count, vertex.weight=gsL, weight.scale=TRUE,
                                 label.edge=FALSE, title.name="Low: number", edge.weight.max=mc))
ggsave("02_cellChat.pdf", c1+c2, width=12, height=6)

# 03 rankNet 通路比较
ggsave("04_rankNet.pdf", rankNet(cc, mode="comparison", stacked=TRUE, do.stat=FALSE),
       width=7, height=6)

# 04 signaling role scatter
ccH <- netAnalysis_computeCentrality(ccH); ccL <- netAnalysis_computeCentrality(ccL)
r1 <- netAnalysis_signalingRole_scatter(ccH) + ggtitle("High ICAM1")
r2 <- netAnalysis_signalingRole_scatter(ccL) + ggtitle("Low ICAM1")
ggsave("05_SignalingRole.pdf", r1+r2, width=11, height=5)

# 05 bubble（全部，两组对比）
ggsave("06_bubble.pdf",
       netVisual_bubble(cc, comparison=c(1,2), angle.x=45, font.size=8, thresh=.05),
       width=12, height=9)
cat("38 高低表达通讯完成\n")
