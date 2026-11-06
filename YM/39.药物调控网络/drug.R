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
# 39 药物富集 + 药物-基因调控网络（核心基因 ICAM1 HPN COL19A1）
library(clusterProfiler); library(igraph); library(tidyr); library(dplyr)
library(ggplot2)

ym <- file.path(root,"YM")
setwd(file.path(ym,"39.药物调控网络"))

core <- readLines(file.path(ym,"23.交集特征基因/core_genes.txt"))
cat("核心基因:", paste(core,collapse=" "),"\n")

# DSigDB（前2列：drug, gene）
drugRT <- read.table(file.path(root,"39.drugEnrich/DSigDB_All_detailed.txt"),
                     header=TRUE, sep="\t", check.names=FALSE,
                     quote="", comment.char="")[,1:2]
colnames(drugRT) <- c("drug","gene")
cat("DSigDB 记录:", nrow(drugRT),"\n")

# 药物富集（输入3基因）
kk <- enricher(core, TERM2GENE=drugRT, pvalueCutoff=1, qvalueCutoff=1,
               minGSSize=1, maxGSSize=5000)
DRUG <- as.data.frame(kk)
DRUG <- DRUG[!is.na(DRUG$pvalue) & DRUG$pvalue < .05, ]
DRUG <- DRUG[order(DRUG$pvalue), ]
write.table(DRUG, "DRUG.enrich.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("显著药物条目:", nrow(DRUG),"\n")

# ---- 构建药物-基因网络（每个基因保留 p 最小的若干药物）----
net <- DRUG[, c("ID","geneID","p.adjust")]
net <- separate_rows(net, geneID, sep="/")
colnames(net) <- c("drug","gene","padj")
net <- net[net$gene %in% core, ]
# 每个基因最多展示 top15 药物
net <- net %>% group_by(gene) %>% arrange(padj) %>% slice_head(n=15) %>% as.data.frame()
write.table(net, "network.txt", sep="\t", quote=FALSE, row.names=FALSE)

edges <- net[, c("drug","gene")]
g <- graph_from_data_frame(edges, directed=FALSE)
vtype <- ifelse(V(g)$name %in% core, "Gene","Drug")
V(g)$type <- vtype
V(g)$shape <- ifelse(vtype=="Gene","circle","square")
V(g)$color <- ifelse(vtype=="Gene","#E64B35","#4DBBD5")
V(g)$size  <- ifelse(vtype=="Gene", 14, 10)
V(g)$label.cex <- .75
E(g)$color <- adjustcolor("grey60",.5)

pdf("Drug_Gene_Network.pdf", width=10, height=8)
set.seed(42)
plot(g, layout=layout_with_fr, vertex.frame.color="white",
     main="Drug-Gene regulatory network")
legend("topright", legend=c("Core gene","Candidate drug"),
       pch=c(19,15), col=c("#E64B35","#4DBBD5"), bty="n")
dev.off()

cat("\nTop 药物:\n")
if(nrow(DRUG)>0) print(head(DRUG[,c("ID","pvalue","Count")],15), row.names=FALSE)
cat("\n已输出药物富集表与网络图\n")
