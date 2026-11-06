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
# 30 免疫浸润柱状图 / 差异箱线 / 相关性热图
library(reshape2); library(ggpubr); library(corrplot)

setwd(file.path(root,"YM/30.免疫差异"))
rt <- read.table(file.path(root,"YM/29.CIBERSORT/CIBERSORT-Results.txt"),
                 header=TRUE, sep="\t", check.names=FALSE)
rownames(rt) <- rt$Mixture
imm <- rt[, 2:23]                      # 22 种免疫细胞
grp <- ifelse(grepl("^PI", rownames(imm)), "Periimplantitis", "Healthy")
# 按分组排序：Healthy 在前，PI 在后
o <- order(grp)
imm <- imm[o, ]; grp <- grp[o]
conNum <- sum(grp=="Healthy"); treatNum <- sum(grp=="Periimplantitis")
data <- t(imm)

# 01 堆叠柱状图
pdf("barplot.pdf", width=13, height=7)
col <- rainbow(nrow(data), s=.7, v=.7)
par(las=1, mar=c(8,5,4,16), mgp=c(3,.1,0), cex.axis=1.3)
a1 <- barplot(data, col=col, xaxt="n", yaxt="n", ylab="Relative Percent", cex.lab=1.4)
a2 <- axis(2, tick=F, labels=F); axis(2, a2, paste0(a2*100,"%"))
par(srt=0, xpd=T)
rect(a1[1]-.5,-.01,a1[conNum]+.5,-.08,col="#0088FF"); text(a1[conNum]/2,-.045,"Healthy",cex=1.6)
rect(a1[conNum]+.5,-.01,a1[length(a1)]+.5,-.08,col="#FF5555"); text((a1[length(a1)]+a1[conNum])/2,-.045,"PI",cex=1.6)
legend(par('usr')[2]*.98, par('usr')[4], legend=rownames(data), col=col, pch=15, bty="n", cex=.9)
dev.off()

# 02 差异箱线图
df <- cbind(as.data.frame(imm), Type=grp)
df <- melt(df, id.vars="Type"); colnames(df) <- c("Type","Immune","Expression")
bx <- ggboxplot(df, x="Immune", y="Expression", fill="Type", ylab="Fraction",
                palette=c(Healthy="#0088FF", Periimplantitis="#FF5555")) +
  rotate_x_text(50) +
  stat_compare_means(aes(group=Type),
                     symnum.args=list(cutpoints=c(0,.001,.01,.05,1), symbols=c("***","**","*","")),
                     label="p.signif")
ggsave("immune.diff.pdf", bx, width=10, height=6)

# 03 PI 组免疫相关性热图
pid <- imm[grp=="Periimplantitis", ]
pid <- pid[, apply(pid,2,sd)>0]
pdf("corHeatmap.pdf", width=11, height=11)
corrplot(cor(pid, method="spearman"), method="color", order="hclust",
         tl.col="black", number.cex=.7, addCoef.col="black",
         col=colorRampPalette(c("blue","white","red"))(50), type="lower")
dev.off()
cat("30 免疫差异三图完成\n")
