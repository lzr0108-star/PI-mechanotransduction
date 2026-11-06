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
# 27 GSEA + 28 GSVA（bulk PI vs H，本地 KEGG GMT）
library(clusterProfiler); library(enrichplot); library(GSVA); library(GSEABase); library(limma)
library(ggplot2); library(pheatmap); library(cowplot)

ym <- file.path(root,"YM")
setwd(file.path(ym,"27.GSEA"))

clinical <- read.csv(file.path(ym,"05.bulk标准化质控/clinical.csv"))
lcpm <- as.matrix(read.table(file.path(ym,"05.bulk标准化质控/expr.logCPM.txt"),
                             header=TRUE, sep="\t", row.names=1, check.names=FALSE))
lcpm <- lcpm[, clinical$sample]
de <- read.delim(file.path(ym,"09.差异分析/DE.PI_vs_H.all.txt"), check.names=FALSE)
kgmt <- read.gmt(file.path(root,"27.GSEA/c2.cp.kegg.Hs.symbols.gmt"))

# ---------- 27 GSEA（按 t 排序）----------
de <- de[!is.na(de$t) & !duplicated(de$gene), ]
gl <- sort(setNames(de$t, de$gene), decreasing=TRUE)
gl <- gl[is.finite(gl)]
gse <- GSEA(gl, TERM2GENE=kgmt, pvalueCutoff=1, pAdjustMethod="BH")
gres <- as.data.frame(gse)
write.table(gres, "GSEA.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("GSEA 通路数:", nrow(gres), " 显著(NES,p<0.05):", sum(gres$p.adjust<.05),"\n")
# top 富集（按 NES 取正/负 top）
gres2 <- gres[order(-abs(gres$NES)),]
topSet <- head(gres2$ID, 6)
pdf("GSEA_plot.pdf", width=12, height=8)
print(gseaplot2(gse, geneSet=topSet, base_size=9))
dev.off()
print(head(gres2[,c("ID","NES","p.adjust")],10), row.names=FALSE)

# ---------- 28 GSVA（ssGSEA）----------
gsc <- getGmt(file.path(root,"27.GSEA/c2.cp.kegg.Hs.symbols.gmt"))
params <- ssgseaParam(exprData=lcpm, geneSets=gsc)
es <- GSVA::gsva(params, verbose=FALSE)
# 规范化名称
rownames(es) <- gsub("KEGG_|HSA_","",rownames(es))
grp <- factor(clinical$group, levels=c("Healthy","Periimplantitis"))
design <- model.matrix(~0+grp); colnames(design)<-levels(grp)
fit <- lmFit(es, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(PI_vs_H=Periimplantitis-Healthy, levels=design)))
gsvaDE <- topTable(fit2, number=Inf)
gsvaDE$ID <- rownames(gsvaDE)
write.table(gsvaDE, file.path(ym,"28.GSVA/GSVA_DE.txt"), sep="\t", quote=FALSE, row.names=FALSE)
# top 差异通路热图
gs <- head(gsvaDE$ID[order(gsvaDE$adj.P.Val)], 20)
mat <- es[gs, , drop=FALSE]
ann <- data.frame(group=clinical$group); rownames(ann)<-clinical$sample
pdf(file.path(ym,"28.GSVA/GSVA_heatmap.pdf"), width=9, height=8)
pheatmap(mat, annotation_col=ann, color=colorRampPalette(c("navy","white","firebrick3"))(50),
         scale="row", fontsize_row=7, fontsize_col=8, border_color=NA)
dev.off()
cat("\nGSVA 差异通路(adj.P<.05):", sum(gsvaDE$adj.P.Val<.05),"\n")
print(head(gsvaDE[,c("ID","logFC","adj.P.Val")],10), row.names=FALSE)
cat("\n已输出 GSEA / GSVA 结果\n")
