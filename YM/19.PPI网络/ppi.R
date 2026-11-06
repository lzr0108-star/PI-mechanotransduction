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
# 19 PPI（817 基因；STRING REST 分批查询 + igraph，top30 hub）
library(httr); library(igraph); library(ggplot2)

setwd(file.path(base, "19.PPI网络"))
gdf <- read.delim(file.path(base,"17.交集基因/interGenes.txt"), check.names=FALSE)
genes <- unique(gdf$gene)

# 随机分批（固定种子，可复现），10 批
set.seed(2026)
bg <- sample(genes)
k <- 10
chunks <- split(bg, cut(seq_along(bg), breaks=k, labels=FALSE))
res <- list()
for(i in seq_along(chunks)){
  idstr <- paste(chunks[[i]], collapse="\r")
  r <- httr::GET("https://string-db.org/api/tsv/network",
                 query=list(species=9606, identifiers=idstr),
                 httr::timeout(90))
  txt <- httr::content(r, as="text", encoding="UTF-8")
  d <- read.delim(text=txt, check.names=FALSE)
  if(nrow(d)>0) res[[i]] <- d
  cat("批次", i, "返回互作", nrow(d), "\n")
  Sys.sleep(0.4)
}
all <- unique(do.call(rbind, res))
cat("\n汇总互作(去重前):", nrow(all), "\n")

# score>=0.4，且两端均在 817 内
all <- all[as.numeric(all$score) >= 0.4, ]
A <- all$preferredName_A; B <- all$preferredName_B
keep <- A %in% genes & B %in% genes & A != B
edges0 <- data.frame(A=A[keep], B=B[keep], score=as.numeric(all$score[keep]))
# 无向去重
key <- apply(edges0[,1:2], 1, function(x) paste(sort(x), collapse="|"))
edges0 <- edges0[!duplicated(key), ]
colnames(edges0)[1:2] <- c("from","to")
write.table(edges0, "PPI_interactions.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("最终互作边数:", nrow(edges0), "\n")

# igraph 节点度
g <- graph_from_data_frame(edges0[,1:2], directed=FALSE)
E(g)$weight <- edges0$score
deg <- sort(degree(g), decreasing=TRUE)
dd <- data.frame(gene=names(deg), degree=as.integer(deg), stringsAsFactors=FALSE)
write.table(dd, "PPI_degree.txt", sep="\t", quote=FALSE, row.names=FALSE)

# top30 barplot
pd <- head(dd,30); pd$gene <- factor(pd$gene, levels=rev(pd$gene))
p <- ggplot(pd, aes(degree, gene, fill=degree)) +
  geom_bar(stat="identity") +
  scale_fill_gradient(low="#F39B7F", high="#E64B35") +
  labs(x="Node degree", y="") +
  theme_bw() + theme(panel.grid=element_blank(), legend.position="none")
ggsave("PPI_degree_barplot.pdf", p, width=7, height=9)

# top30 网络
top <- names(deg)[1:30]
g2 <- induced_subgraph(g, top)
pdf("PPI_network_top30.pdf", width=11, height=10)
set.seed(42); wc <- cluster_walktrap(g2)
plot(g2, layout=layout_with_fr, vertex.size= 5+sqrt(degree(g2))*2.2,
     vertex.color=membership(wc), vertex.label.cex=.62,
     vertex.label.color="black", edge.color=adjustcolor("grey60",.4))
title("PPI network (top30 hub)")
dev.off()

cat("\nTop30 hub:\n"); cat(paste(top, collapse=" "), "\n")
cat("已输出互作表/度/barplot/network\n")
