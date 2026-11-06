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
# 37 虚拟敲除 scTenifoldKnk（敲除 ICAM1）
library(scTenifoldKnk); library(Seurat); library(ggplot2)
library(dplyr); library(ggrepel)
set.seed(12345)

setwd(file.path(root,"YM/37.虚拟敲除"))
seu <- readRDS(file.path(root,"YM/12.Seurat分析/seu_annotated.rds"))
seu <- subset(seu, subset=condition=="Periimplantitis")
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])

target <- "ICAM1"
countMat <- GetAssayData(seu, layer="counts")
seu <- FindVariableFeatures(seu, nfeatures=5000)
hvgs <- VariableFeatures(seu)
keep <- unique(c(target, hvgs))
keep <- intersect(keep, rownames(countMat))
data <- as.matrix(countMat[keep, ])
cat("输入基因:", length(keep), " 细胞:", ncol(data), "\n")

result <- scTenifoldKnk(countMatrix=data, gKO=target,
                       qc=TRUE, qc_minLibSize=500, qc_maxMTratio=.1,
                       nc_nNet=10, nc_nCells=500, nCores=4)
df <- result$diffRegulation
df <- df[df$gene != target, ]
outTab <- df[!is.na(df$p.adj) & df$p.adj<.05, ]
write.table(outTab, "sigDiff.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("显著调控基因:", nrow(outTab), "\n")

# 柱状图 top20
top <- head(df[order(-df$FC),], 20)
p1 <- ggplot(top, aes(x=reorder(gene,FC), y=FC)) +
  geom_bar(stat="identity", fill="#5A9BD4") + coord_flip() +
  labs(x="Gene", y="FC", title="Top 20 regulated (ICAM1 KO)") +
  theme_minimal() + theme(plot.title=element_text(hjust=.5))
ggsave("barplot.pdf", p1, width=7, height=5.5)

# 火山图
df$log_p.adj <- -log10(df$p.adj)
df$significant <- ifelse(!is.na(df$p.adj)&df$p.adj<.05,"Significant","NS")
lab <- subset(df, significant=="Significant")
yu <- quantile(df$log_p.adj,.999,na.rm=TRUE)
p2 <- ggplot(df, aes(x=Z,y=log_p.adj,color=significant)) +
  geom_point(alpha=.7,size=1) +
  scale_color_manual(values=c(Significant="red",NS="gray70")) +
  geom_hline(yintercept=-log10(.05),linetype="dashed",color="red") +
  geom_text_repel(data=lab, aes(label=gene),size=3,max.overlaps=40) +
  coord_cartesian(ylim=c(0,yu)) +
  labs(x="Z-score",y="-log10(p.adj)") +
  theme_classic() + theme(legend.position="none")
ggsave("vol.pdf", p2, width=7, height=5.5)

# 饼图
sc <- table(df$significant); sdf <- as.data.frame(sc)
colnames(sdf) <- c("category","count")
sdf$percentage <- paste0(round(sdf$count/sum(sdf$count)*100,1),"%")
p3 <- ggplot(sdf,aes(x="",y=count,fill=category)) +
  geom_bar(stat="identity",width=1) +
  geom_text(aes(label=percentage),position=position_stack(vjust=.5),size=4) +
  coord_polar("y",start=0) +
  scale_fill_manual(values=c(NS="lightgray",Significant="red")) +
  labs(title="Significant proportion",fill="") +
  theme_minimal()+theme(axis.title=element_blank(),axis.text=element_blank(),
                       panel.grid=element_blank())
ggsave("pie.pdf", p3, width=6, height=5.5)
cat("37 虚拟敲除完成\n")
