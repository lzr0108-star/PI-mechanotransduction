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
# 25 独立验证集 GSE33774：核心基因差异表达 + ROC
library(ggplot2); library(ggpubr); library(reshape2)

setwd(file.path(root,"YM/25.验证集差异"))

# ---------- 解析 series matrix ----------
ml <- readLines(gzfile("GSE33774_series_matrix.txt.gz"))
field <- function(key){
  l <- ml[grepl(key, ml)][1]
  gsub('"','', strsplit(l,"\t")[[1]][-1])
}
src <- field("!Sample_source_name_ch1")
acc <- field("!Sample_geo_accession")
t1 <- grep("!series_matrix_table_begin", ml)
t2 <- grep("!series_matrix_table_end", ml)
mat <- read.delim(textConnection(ml[(t1+1):(t2-1)]),
                  check.names=FALSE, quote='"')
rownames(mat) <- mat[,1]
expr <- as.matrix(mat[, -1])
colnames(expr) <- acc
grp <- ifelse(grepl("healthy", src, ignore.case=TRUE), "Healthy",
       ifelse(grepl("periimplantitis", src, ignore.case=TRUE), "PI", NA))
expr <- expr[, !is.na(grp)]; grp <- grp[!is.na(grp)]
cat("Healthy:", sum(grp=="Healthy"), " PI:", sum(grp=="PI"), "\n")

# ---------- 解析 annot：ID -> Gene symbol ----------
al <- readLines("GPL6244.annot")
a1 <- grep("!platform_table_begin", al)
aend <- grep("platform_table_end", al); if(length(aend)==0) aend <- length(al)+1
adf <- read.delim(textConnection(al[(a1+1):(aend-1)]),
                  check.names=FALSE, quote='"')
cat("annot 探针:", nrow(adf), "\n")
ids <- intersect(rownames(expr), as.character(adf$ID))
sym <- adf[match(ids, as.character(adf$ID)), "Gene symbol"]
ann <- data.frame(ID=ids, symbol=as.character(sym), stringsAsFactors=FALSE)
ann <- ann[ann$symbol!="" & !grepl("---", ann$symbol), ]
# 一个 symbol 可能多探针：取均值
sub <- expr[ann$ID, , drop=FALSE]
glev <- strsplit(ann$symbol, " /// ")
geneSym <- rep(ann$symbol, sapply(glev,length))
# 简化：直接用主 symbol（不拆分多基因），按 symbol 聚合探针
agg <- rowsum(sub, group=ann$symbol) / as.vector(table(ann$symbol)[ann$symbol])
# 上式近似；改用更稳妥的 tapply 聚合
geneExpr <- apply(sub, 2, function(v) tapply(v, INDEX=ann$symbol, FUN=mean))
cat("注释到基因:", nrow(geneExpr), "\n")

# ---------- 核心基因验证 ----------
hub <- c("ICAM1","HPN","COL19A1")
found <- hub[hub %in% rownames(geneExpr)]
cat("验证集中找到:", paste(found, collapse=", "), "\n")
df <- data.frame(t(geneExpr[found, , drop=FALSE]), Group=grp, check.names=FALSE)
df <- melt(df, id.vars="Group"); colnames(df) <- c("Group","Gene","Expression")

p1 <- ggboxplot(df, x="Group", y="Expression", fill="Group",
                palette=c(Healthy="#0088FF", PI="#FF5555"), add="jitter") +
  facet_wrap(~Gene, scales="free_y", nrow=1) +
  stat_compare_means(comparisons=list(c("Healthy","PI")), method="wilcox.test",
                     label="p.format", size=3.5) +
  theme(axis.text.x=element_text(angle=30,hjust=1))
ggsave("Validation_boxplot.pdf", p1, width=9, height=5)

# ---------- ROC ----------
if(requireNamespace("pROC", quietly=TRUE)){
  library(pROC); library(patchwork)
  rocs <- lapply(found, function(g){
    x <- as.numeric(geneExpr[g, ])
    r <- roc(grp, x, quiet=TRUE, levels=c("Healthy","PI"))
    data.frame(TPR=r$sensitivities, FPR=1-r$specificities,
               Gene=paste0(g," (AUC=",sprintf("%.3f",as.numeric(auc(r))),")"),
               AUC=as.numeric(auc(r)))
  })
  rdf <- do.call(rbind, rocs)
  p2 <- ggplot(rdf, aes(FPR,TPR,color=Gene)) +
    geom_line(linewidth=.8) + geom_abline(slope=1,intercept=0,lty=2,color="gray60") +
    geom_point(data=subset(rdf, FPR %in% c(0,1)), size=1.5) +
    scale_color_manual(values=c("#d62728","#1f77b4","#2ca02c")) +
    labs(x="1 - Specificity", y="Sensitivity", title="GSE33774 validation ROC") +
    theme_bw() + theme(plot.title=element_text(hjust=.5), legend.position=c(.72,.18))
  ggsave("Validation_ROC.pdf", p2, width=6.5, height=6)
  print(rdf[!duplicated(rdf$Gene), c("Gene","AUC")], row.names=FALSE)
}
cat("25 独立验证完成\n")
