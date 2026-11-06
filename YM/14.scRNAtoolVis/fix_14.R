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
# 修复 14 章 03 / 04 标签阻挡
library(dplyr); library(ggplot2); library(scRNAtoolVis)
library(jjAnno); library(ggsci)

setwd(file.path(root,"YM/14.scRNAtoolVis"))
mk <- read.delim("allGroupMarkers.txt", check.names=FALSE)
tilecol <- jjAnno::useMyCol("paired", n=12)

# ===== 03 markerVolcano：减少标注（topn 5->4），标签下移避免压轴 =====
pdf("03_markerVolcano.pdf", width=12, height=6.4)
print(markerVolcano(markers=mk, topn=4, labelCol=pal_npg()(10), log2FC=1,
                    base_size=12, nudge_x=.12, pnudge_y=.6, nnudge_y=-.6))
dev.off()

# ===== 04 jjVolcano：topGeneN 减少，celltypeSize 加大，底部留白 =====
pdf("04_jjVolcano1.pdf", width=12.5, height=6.8)
print(jjVolcano(diffData=mk, legend.position=c(.07,.96), tile.col=tilecol,
                log2FC.cutoff=1, adjustP.cutoff=.05, topGeneN=4)
      + labs(x="Cell type"))
dev.off()

pdf("04_jjVolcano2.pdf", width=9.5, height=7.6)
print(jjVolcano(diffData=mk, legend.position=c(.05,.96), flip=TRUE,
                tile.col=tilecol, log2FC.cutoff=1, adjustP.cutoff=.05, topGeneN=4) + labs(x="Cell type"))
dev.off()

pdf("04_jjVolcano3.pdf", width=8.4, height=7.6)
print(jjVolcano(diffData=mk, legend.position=c(.1,.86), polar=TRUE, size=3,
                tile.col=tilecol, log2FC.cutoff=1, adjustP.cutoff=.05, topGeneN=4) + labs(x="Cell type"))
dev.off()
cat("14 fixed\n")
