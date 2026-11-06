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
# 16 WGCNA（H+PI，20 样本，二组）
library(WGCNA); library(limma)
WGCNA::disableWGCNAThreads()
options(stringsAsFactors=FALSE)

setwd(file.path(base, "16.WGCNA"))

logcpm <- as.matrix(read.table(file.path(base,"05.bulk标准化质控/expr.logCPM.txt"),
                               header=TRUE, sep="\t", row.names=1, check.names=FALSE))
sdv <- apply(logcpm, 1, sd)
datExpr <- t(logcpm[sdv > 0.8, ])
cat("保留基因:", ncol(datExpr), " 样本:", nrow(datExpr), "\n")
gsg <- goodSamplesGenes(datExpr, verbose=0)
if(!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

sampleTree <- hclust(dist(datExpr), "average")
pdf("1_sample_cluster.pdf", width=10, height=5)
plot(sampleTree, main="Sample clustering", xlab="", sub="")
dev.off()

clinical <- read.csv(file.path(base,"05.bulk标准化质控/clinical.csv"))
sm <- rownames(datExpr)
grp <- clinical$group[match(sm, clinical$sample)]
trait <- data.frame(Healthy=as.integer(grp=="Healthy"),
                    Periimplantitis=as.integer(grp=="Periimplantitis"))
rownames(trait) <- sm
pdf("2_sample_heatmap.pdf", width=10, height=6)
plotDendroAndColors(sampleTree, numbers2colors(trait, signed=FALSE),
                    groupLabels=names(trait), main="Sample dendrogram and traits")
dev.off()

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=5)
softPower <- sft$powerEstimate
if(is.na(softPower)) softPower <- 10
cat("soft power:", softPower, "\n")
pdf("3_scale_independence.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold", ylab="SFT R^2", type="n", main="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=.85, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold",
     ylab="Mean Connectivity", type="n", main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

adj <- adjacency(datExpr, power=softPower)
TOM <- TOMsimilarity(adj)
dissTOM <- 1-TOM
geneTree <- hclust(as.dist(dissTOM), "average")
pdf("4_gene_clustering.pdf", width=9, height=6)
plot(geneTree, xlab="", sub="", main="Gene clustering on TOM", labels=FALSE, hang=.04)
dev.off()

dynamicMods <- cutreeDynamic(dendro=geneTree, distM=dissTOM, deepSplit=2,
                             pamRespectsDendro=FALSE, minClusterSize=50)
dynamicColors <- labels2colors(dynamicMods)
pdf("5_Dynamic_Tree.pdf", width=9, height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels=FALSE,
                    hang=.03, addGuide=TRUE, guideHang=.05, main="Dynamic Tree")
dev.off()

MEs0 <- moduleEigengenes(datExpr, dynamicColors)$eigengenes
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight=.25, verbose=0)
mergedColors <- merge$colors; mergedMEs <- merge$newMEs
pdf("6_Clustering_module.pdf", width=8, height=6)
plot(hclust(as.dist(1-cor(MEs0)), "average"), main="Module eigengene clustering", xlab="")
abline(h=.25, col="red")
dev.off()
pdf("7_merged_dynamic.pdf", width=9, height=6)
plotDendroAndColors(geneTree, mergedColors, "Merged", dendroLabels=FALSE,
                    hang=.03, addGuide=TRUE, guideHang=.05, main="Merged dynamic")
dev.off()

moduleTraitCor <- cor(mergedMEs, trait, use="p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",signif(moduleTraitP,1),")")
dim(textMatrix) <- dim(moduleTraitCor)
pdf("8_Module_trait.pdf", width=6, height=7)
par(mar=c(6,10,3,3))
labeledHeatmap(Matrix=moduleTraitCor, xLabels=names(trait), yLabels=names(mergedMEs),
               ySymbols=names(mergedMEs), colorLabels=FALSE, colors=blueWhiteRed(50),
               textMatrix=textMatrix, setStdMargins=FALSE, cex.text=.85, zlim=c(-1,1),
               main="Module-trait relationships")
dev.off()

pi.cor <- moduleTraitCor[,"Periimplantitis"]; pi.p <- moduleTraitP[,"Periimplantitis"]
cand <- names(pi.cor)[pi.cor>0.3 & pi.p<.05]; cand <- cand[order(-pi.cor[cand])]
keyModule <- if(length(cand)) cand[1] else names(which.max(pi.cor))
cat("\nPI 模块相关:\n"); print(round(pi.cor[order(-pi.cor)],3))
cat("★ PI 关键模块:", keyModule, " r=", round(pi.cor[keyModule],3),
    " p=", signif(pi.p[keyModule],2), "\n")

GS1 <- abs(as.numeric(cor(trait$Periimplantitis, datExpr, use="p")))
pdf("9_GeneSignificance.pdf", width=12, height=6)
plotModuleSignificance(GS1, mergedColors)
dev.off()

genes <- colnames(datExpr)
for(m in sort(unique(mergedColors))){
  write.table(genes[mergedColors==m], paste0("module_",m,".txt"),
              sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}
writeLines(keyModule, "keyModule.txt")
saveRDS(list(datExpr=datExpr, trait=trait, mergedColors=mergedColors, mergedMEs=mergedMEs,
             moduleTraitCor=moduleTraitCor, moduleTraitP=moduleTraitP, keyModule=keyModule,
             softPower=softPower), "wgcna.rds")
cat("已输出全部 WGCNA 结果\n")
