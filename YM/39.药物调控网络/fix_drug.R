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
# 重做 39 药物-基因网络：双中心星型布局，标签外排
library(igraph); library(dplyr)

setwd(file.path(root,"YM/39.药物调控网络"))
core <- readLines(file.path(root,"YM/23.交集特征基因/core_genes.txt"))
net <- read.table("network.txt", header=TRUE, sep="\t", check.names=FALSE)

# 每基因 top8；过滤超长/低质量药名（IUPAC 全名），过长截断
net <- net %>% group_by(gene) %>% arrange(padj) %>% slice_head(n=8) %>% ungroup()
shortName <- function(x){
  x <- gsub("\\(.*","", x)
  x <- trimws(x)
  ifelse(nchar(x)>22, paste0(substr(x,1,21),"~"), x)
}
net$drug <- shortName(net$drug)
net <- net[nchar(net$drug)>=2, ]
net <- net[!duplicated(net$drug), ]

genes <- unique(net$gene)
k <- length(genes)
gx <- seq(-1,1,length=k) * 2.6
gy <- rep(0, k)
gcenter <- cbind(gx, gy)
rownames(gcenter) <- genes

# 计算每个药物坐标（围绕其基因，半径 1.55，全圆周）
pos <- list()
for(i in seq_len(k)){
  sub <- net[net$gene==genes[i], ]
  n <- nrow(sub)
  ang <- seq(0, 2*pi, length=n+1)[seq_len(n)] + pi/2
  c0 <- gcenter[i,]
  xy <- cbind(c0[1] + 1.55*cos(ang), c0[2] + 1.45*sin(ang))
  rownames(xy) <- sub$drug
  pos[[i]] <- xy
}
drugPos <- do.call(rbind, pos)
allPos <- rbind(gcenter, drugPos)

# 构建图（仅用于画点/边）
g <- graph_from_data_frame(net[,c("drug","gene")], directed=FALSE,
                          vertices=data.frame(name=c(genes, rownames(drugPos))))
vtype <- ifelse(V(g)$name %in% genes, "Gene","Drug")

pdf("Drug_Gene_Network.pdf", width=13, height=8.5)
par(mar=c(2,2,3,2))
plot(g, layout=allPos[V(g)$name, ], rescale=FALSE,
     vertex.shape=ifelse(vtype=="Gene","circle","square"),
     vertex.color=ifelse(vtype=="Gene","#E64B35","#4DBBD5"),
     vertex.size=ifelse(vtype=="Gene",36,12),
     vertex.frame.color="white", vertex.label=NA,
     edge.color=adjustcolor("grey60",.55), edge.width=.9,
     xlim=c(-4.9,4.9), ylim=c(-2.9,2.9),
     main="Drug-Gene regulatory network")
# 基因标签
for(i in seq_len(k)) text(gx[i], gy[i], genes[i], cex=.78, font=2, col="white")
# 药物标签沿径向向外
for(i in seq_len(k)){
  sub <- net[net$gene==genes[i], ]
  n <- nrow(sub)
  ang <- seq(0,2*pi,length=n+1)[seq_len(n)] + pi/2
  lx <- gx[i] + 1.82*cos(ang); ly <- gy[i] + 1.72*sin(ang)
  ha <- ifelse(cos(ang)>0, 0, 1)
  for(j in seq_len(n)) text(lx[j], ly[j], sub$drug[j], cex=.72, adj=c(ha[j],.5))
}
legend("topright", legend=c("Core gene","Candidate drug"),
       pch=c(19,15), col=c("#E64B35","#4DBBD5"), bty="n", cex=.9)
dev.off()
cat("Drug network 重做完成，药物数", nrow(net), "\n")
