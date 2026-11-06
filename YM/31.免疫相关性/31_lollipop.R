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
# 31 核心基因-免疫细胞相关性棒棒糖图
library(ggplot2); library(dplyr)

setwd(file.path(root,"YM/31.免疫相关性"))

# bulk logCPM
ex <- read.table(file.path(root,"YM/05.bulk标准化质控/expr.logCPM.txt"),
                 header=TRUE, sep="\t", check.names=FALSE, row.names=1)
# CIBERSORT
im <- read.table(file.path(root,"YM/29.CIBERSORT/CIBERSORT-Results.txt"),
                 header=TRUE, sep="\t", check.names=FALSE)
rownames(im) <- im$Mixture
imm <- im[, 2:23]

genes <- c("ICAM1","HPN","COL19A1")
samples <- intersect(colnames(ex), rownames(imm))
cat("对齐样本:", length(samples), "\n")
ex <- ex[, samples]; imm <- imm[samples, ]

res <- data.frame()
for(g in genes){
  x <- as.numeric(ex[g, ])
  for(cell in colnames(imm)){
    y <- imm[[cell]]
    if(sd(y)==0) y[1] <- 1e-5
    ct <- suppressWarnings(cor.test(x, y, method="spearman"))
    res <- rbind(res, data.frame(Gene=g, Cell=cell,
                                 cor=unname(ct$estimate), pvalue=ct$p.value))
  }
}
write.table(res, "core_immune_corP.txt", sep="\t", quote=FALSE, row.names=FALSE)

# 棒棒糖图（每基因一张）
pcol <- c('darkgreen','LimeGreen','orange','pink','gold')
for(g in genes){
  d <- res[res$Gene==g, ]
  d$Cell <- factor(d$Cell, levels=d$Cell[order(d$cor)])
  d$pvlab <- ifelse(d$pvalue<.001,"<.001",sprintf("%.3f",d$pvalue))
  p <- ggplot(d, aes(x=cor, y=Cell)) +
    geom_vline(xintercept=0, color="gray60") +
    geom_segment(aes(xend=0, yend=Cell), linewidth=1, color="gray50") +
    geom_point(aes(size=abs(cor), color=pvalue)) +
    geom_text(aes(x=cor+ifelse(d$cor>=0,.09,-.09), label=pvlab,
                  hjust=ifelse(d$cor>=0,0,1)), size=2.8,
              color=ifelse(d$pvalue<.05,"red","black")) +
    scale_color_gradientn(colors=pcol, limits=c(0,1),
                          breaks=c(.2,.4,.6,.8,1)) +
    scale_size_continuous(limits=c(0,1), range=c(1,5)) +
    scale_x_continuous(limits=c(-1.15,1.15)) +
    labs(title=paste0(g," vs immune cells"), x="Spearman r", y="") +
    theme_bw() +
    theme(plot.title=element_text(hjust=.5),
          axis.text.y=element_text(size=8))
  ggsave(paste0("Lollipop.",g,".pdf"), p, width=9, height=7)
}
cat("31 免疫相关性棒棒糖图完成\n")
