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
# 18 GO/KEGG 富集（817 交集基因；GO本地org库，KEGG本地GMT，全程离线）
library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2); library(circlize)

base <- file.path(root, "YM")
setwd(file.path(base, "18.富集分析"))
gdf <- read.delim(file.path(base,"17.交集基因/interGenes.txt"), check.names=FALSE)
genes <- unique(gdf$gene)

# ---- GO（本地）----
eg <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
ego <- enrichGO(eg$ENTREZID, OrgDb=org.Hs.eg.db, ont="all", readable=TRUE,
                pvalueCutoff=.05, qvalueCutoff=.5)
gores <- as.data.frame(ego)
write.table(gores, "GO.txt", sep="\t", quote=FALSE, row.names=FALSE)

# ---- KEGG（本地 GMT + enricher）----
kgmt <- read.gmt(file.path(root,"27.GSEA/c2.cp.kegg.Hs.symbols.gmt"))
# 规范 term 名（去掉 KEGG_ 前缀，保留可读名）
ekk <- enricher(genes, TERM2GENE=kgmt, pvalueCutoff=.1, qvalueCutoff=.5)
kres <- as.data.frame(ekk)
if(nrow(kres)>0) kres$Description <- gsub("KEGG_|HSA_","",kres$ID)
write.table(kres, "KEGG.txt", sep="\t", quote=FALSE, row.names=FALSE)

# ---- 分组 barplot ----
mkdf <- function(d, ont, n=7){
  if(is.null(d)||nrow(d)==0) return(NULL)
  d <- head(d[order(d$p.adjust),], n)
  data.frame(Description=d$Description, p.adjust=d$p.adjust, ont=ont,
             Count=as.integer(d$Count), stringsAsFactors=FALSE)
}
bp <- mkdf(gores[gores$ONTOLOGY=="BP",],"BP")
cc <- mkdf(gores[gores$ONTOLOGY=="CC",],"CC")
mf <- mkdf(gores[gores$ONTOLOGY=="MF",],"MF")
kp <- if(nrow(kres)) mkdf(kres,"KEGG") else NULL
alld <- do.call(rbind,list(bp,cc,mf,kp))
alld$ont <- factor(alld$ont, levels=c("BP","CC","MF","KEGG"))
alld$Description <- factor(alld$Description, levels=rev(unique(alld$Description)))
cols <- c(BP="#E64B35",CC="#4DBBD5",MF="#00A087",KEGG="#3C5488")
p <- ggplot(alld, aes(-log10(p.adjust),Description,fill=ont)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=cols) + labs(x="-log10(p.adjust)",y="") +
  theme_bw() + theme(panel.grid=element_blank())
ggsave("barplot.pdf", p, width=10, height=9)

# ---- circlize 环形图（BP top10）----
cir <- head(gores[gores$ONTOLOGY=="BP",][order(gores[gores$ONTOLOGY=="BP",]$p.adjust),],10)
glist <- strsplit(as.character(cir$geneID),"/")
ug <- unique(unlist(glist)); nm<-length(ug)
lab <- paste0("T",1:nrow(cir))
tcol <- colorRampPalette(c("#E64B35","#F39B7F","#4DBBD5","#00A087","#3C5488"))(nrow(cir))
gcol <- rep("grey88",nm)
sectors <- c(lab, paste0("G",1:nm))
pdf("circlize.pdf", width=9, height=9)
par(mar=c(0,0,0,0)); circos.clear()
circos.par(gap.degree=c(1), start.degree=90)
circos.initialize(factors=sectors, xlim=cbind(rep(0,length(sectors)),rep(1,length(sectors))))
circos.track(ylim=c(0,1), track.height=.12, bg.col=c(tcol,gcol), bg.border=NA)
circos.track(ylim=c(0,1), track.height=.06, bg.col="white", bg.border=NA,
             panel.fun=function(x,y){
               i<-get.cell.meta.data("sector.numeric.index")
               if(i<=nrow(cir)) circos.text(.5,.5,paste0(strwrap(cir$Description[i],14),collapse="\n"),
                 cex=.62, niceFacing=TRUE, adj=c(.5,.5))
             })
for(i in 1:nrow(cir)){
  gs <- match(glist[[i]],ug)
  for(gg in gs) circos.link(lab[i],.5,paste0("G",gg),.5,
                           col=adjustcolor(tcol[i],.28), h.ratio=.4)
}
dev.off(); circos.clear()

cat("BP top:\n"); print(head(gores[gores$ONTOLOGY=="BP",c("Description","p.adjust")],10),row.names=FALSE)
cat("\nKEGG:\n"); if(nrow(kres)) print(head(kres[,c("Description","p.adjust")],10),row.names=FALSE)
cat("\n已输出 GO/KEGG 表、barplot.pdf、circlize.pdf\n")
