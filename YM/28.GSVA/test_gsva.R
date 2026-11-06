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
library(GSVA); library(GSEABase); library(limma); library(pheatmap)
ym <- file.path(root,"YM")
setwd(file.path(ym,"28.GSVA"))
clinical <- read.csv(file.path(ym,"05.bulk标准化质控/clinical.csv"))
lcpm <- as.matrix(read.table(file.path(ym,"05.bulk标准化质控/expr.logCPM.txt"),
                             header=TRUE, sep="\t", row.names=1, check.names=FALSE))
lcpm <- lcpm[, clinical$sample]
cat("expr:", nrow(lcpm),"x",ncol(lcpm),"\n")
gsc <- getGmt(file.path(root,"27.GSEA/c2.cp.kegg.Hs.symbols.gmt"))
cat("gene sets:", length(gsc),"\n")
params <- ssgseaParam(exprData=lcpm, geneSets=gsc)
es <- GSVA::gsva(params, verbose=FALSE)
cat("ES matrix:", nrow(es),"x",ncol(es),"\n")
rownames(es) <- gsub("KEGG_|HSA_","",rownames(es))
grp <- factor(clinical$group, levels=c("Healthy","Periimplantitis"))
design <- model.matrix(~0+grp); colnames(design)<-levels(grp)
fit <- lmFit(es, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(PI_vs_H=Periimplantitis-Healthy, levels=design)))
gsvaDE <- topTable(fit2, number=Inf); gsvaDE$ID<-rownames(gsvaDE)
write.table(gsvaDE, "GSVA_DE.txt", sep="\t", quote=FALSE, row.names=FALSE)
gs <- head(gsvaDE$ID[order(gsvaDE$adj.P.Val)],20)
ann <- data.frame(group=clinical$group); rownames(ann)<-clinical$sample
pdf("GSVA_heatmap.pdf", width=9, height=8)
pheatmap(es[gs,,drop=FALSE], annotation_col=ann,
         color=colorRampPalette(c("navy","white","firebrick3"))(50),
         scale="row", fontsize_row=7, fontsize_col=8, border_color=NA)
dev.off()
cat("GSVA 显著差异通路:", sum(gsvaDE$adj.P.Val<.05),"\n")
print(head(gsvaDE[,c("ID","logFC","adj.P.Val")],10), row.names=FALSE)
