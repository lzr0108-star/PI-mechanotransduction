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
# 重新生成 18 章 circlize.pdf（读 GO.txt，BP top10）
library(circlize)

setwd(file.path(root,"YM/18.富集分析"))
gores <- read.delim("GO.txt", check.names=FALSE)
cat("GO 行数:", nrow(gores), " 列:", paste(colnames(gores),collapse=","), "\n")

bp <- gores[gores$ONTOLOGY=="BP",]
cir <- head(bp[order(bp$p.adjust),], 10)
cat("BP term 数:", nrow(cir), "\n")

glist0 <- strsplit(as.character(cir$geneID), "/")
# 只保留跨多 term 出现的关键基因，限制总数，避免 sector 过密
allg <- unlist(glist0)
freq <- sort(table(allg), decreasing=TRUE)
ug <- names(freq)[seq_len(min(45, length(freq)))]
glist <- lapply(glist0, function(x) intersect(x, ug))
nm <- length(ug)
lab <- paste0("T", seq_len(nrow(cir)))
tcol <- colorRampPalette(c("#E64B35","#F39B7F","#4DBBD5","#00A087","#3C5488"))(nrow(cir))
gcol <- rep("grey90", nm)
sectors <- c(lab, paste0("G", seq_len(nm)))

circos.clear()
pdf("circlize.pdf", width=9.5, height=9.5)
par(mar=c(0,0,0,0))
gdeg <- min(0.8, 240/length(sectors))
circos.par(gap.degree=gdeg, start.degree=90)
circos.initialize(factors=sectors, xlim=cbind(rep(0,length(sectors)), rep(1,length(sectors))))
circos.track(ylim=c(0,1), track.height=.12, bg.col=c(tcol,gcol), bg.border=NA)
circos.track(ylim=c(0,1), track.height=.07, bg.col="white", bg.border=NA,
             panel.fun=function(x,y){
               i <- get.cell.meta.data("sector.numeric.index")
               if(i <= nrow(cir))
                 circos.text(.5,.5, lab[i], cex=.95, font=2,
                             niceFacing=TRUE, adj=c(.5,.5))
             })
for(i in seq_len(nrow(cir))){
  gs <- match(glist[[i]], ug)
  for(gg in gs)
    circos.link(lab[i], .5, paste0("G",gg), .5,
                col=adjustcolor(tcol[i],.26), h.ratio=.4)
}
dev.off()
circos.clear()
writeLines(paste0(lab,": ",cir$Description), "term_legend.txt")
cat("circlize 重新生成完成\n")
