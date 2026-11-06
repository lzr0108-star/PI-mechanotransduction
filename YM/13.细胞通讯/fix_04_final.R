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
# 04 最终：通路 x 细胞 outgoing/incoming 热图（pheatmap，Healthy vs PI）
library(pheatmap)

setwd(file.path(root,"YM/13.细胞通讯"))
cc <- readRDS("cellchat_compare.rds")
npH <- cc@netP$Healthy
npP <- cc@netP$PI
probH <- npH$prob
probP <- npP$prob

outH <- apply(probH, c(3,1), sum)
outP <- apply(probP, c(3,1), sum)
inH  <- apply(probH, c(3,2), sum)
inP  <- apply(probP, c(3,2), sum)

pw <- intersect(rownames(outH), rownames(outP))
cells <- intersect(colnames(outH), colnames(outP))
outH <- outH[pw, cells]; outP <- outP[pw, cells]
inH  <- inH[pw, cells];  inP  <- inP[pw, cells]

reds <- colorRampPalette(c("#fff5f5","#fee0d2","#fc9272","#de2d26"))(60)
blues<- colorRampPalette(c("#f7fbff","#deebf7","#6baed6","#08519c"))(60)
bR <- c(0, max(c(outH,outP)))
bB <- c(0, max(c(inH,inP)))

hm <- function(mat, cols, brks, ttl){
  pheatmap(mat, color=cols, breaks=seq(brks[1],brks[2],length=61),
           cluster_rows=TRUE, cluster_cols=FALSE,
           treeheight_row=14, fontsize=9, fontsize_row=9,
           border_color=NA, main=ttl, angle_col=45)
}
pdf("04_pathway_heatmap.pdf", width=12, height=11)
hm(outH, reds, bR, "Healthy: outgoing (pathway x cell)")
hm(outP, reds, bR, "Periimplantitis: outgoing (pathway x cell)")
hm(inH, blues, bB, "Healthy: incoming (pathway x cell)")
hm(inP, blues, bB, "Periimplantitis: incoming (pathway x cell)")
dev.off()
cat("04 final done, pathways", length(pw), "\n")
