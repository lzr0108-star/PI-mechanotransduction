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
# 07 箱线图（H+PI，TMM 前后，二组着色）
library(edgeR); library(ggplot2); library(reshape2); library(cowplot)

setwd(file.path(base, "07.箱线图"))
y <- readRDS(file.path(base,"05.bulk标准化质控/dge.rds"))
clinical <- read.csv(file.path(base,"05.bulk标准化质控/clinical.csv"))

ord <- clinical$sample[order(clinical$group, clinical$sample)]
group <- clinical$group[match(ord, clinical$sample)]
before <- log2(y$count[, ord] + 1)
after  <- cpm(y, log=TRUE, prior.count=2)[, ord]

mk <- function(mat, title){
  df <- melt(mat); colnames(df) <- c("gene","sample","expression")
  df$group <- group[match(df$sample, ord)]
  ggplot(df, aes(x=sample, y=expression, fill=group)) +
    geom_boxplot(outlier.size=0.1, outlier.alpha=0.3, width=0.7) +
    scale_fill_manual(values=c("Healthy"="#4DBBD5","Periimplantitis"="#E64B35")) +
    labs(title=title, x="", y="Expression") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=7),
          plot.title=element_text(hjust=0.5), panel.grid=element_blank())
}
pdf("boxplot.pdf", width=12, height=10)
plot_grid(mk(before,"Before TMM"), mk(after,"After TMM"), ncol=1)
dev.off()
cat("已输出 boxplot.pdf\n")
