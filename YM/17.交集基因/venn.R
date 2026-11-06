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
# 17 三方交集：MRG ∩ PI-DEG ∩ WGCNA red（H+PI，无牙周炎）
library(ggvenn)

ym <- file.path(root, "YM")
setwd(file.path(ym, "17.交集基因"))

red <- readLines(file.path(ym,"16.WGCNA/module_red.txt"))
MRG <- readLines(file.path(root,"17.venn/MRG.txt"))
degTab <- read.delim(file.path(ym,"09.差异分析/DE.PI_vs_H.sig.txt"), check.names=FALSE)
DEG <- as.character(degTab$gene)

sets <- list(MRG=unique(MRG), DEG=unique(DEG), WGCNA_red=unique(red))
cat("集合大小: MRG=%d DEG=%d red=%d\n", length(sets$MRG), length(sets$DEG), length(sets$WGCNA_red))

inter <- Reduce(intersect, list(sets$MRG, sets$DEG, sets$WGCNA_red))
cat("★ 三方交集基因数:", length(inter), "\n")
cat(paste(sort(inter), collapse=" "), "\n")
write.table(data.frame(gene=sort(inter)), "interGenes.txt", sep="\t", quote=FALSE, row.names=FALSE)

p <- ggvenn(sets, fill_color=c("#E64B35","#4DBBD5","#00A087"),
            stroke_size=.6, set_name_size=4.5, text_size=3.5)
ggsave("venn.pdf", p, width=6.5, height=6)
cat("已输出 interGenes.txt / venn.pdf\n")
