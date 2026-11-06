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
# 20-23 机器学习 —— 稳定性选择（stability selection，B=100 次重抽样）
library(glmnet); library(e1071); library(randomForest); library(ggplot2)
library(ggvenn)

ym <- file.path(root,"YM")
setwd(file.path(ym,"20.LASSO"))

lcpm0 <- as.matrix(read.table(file.path(ym,"05.bulk标准化质控/expr.logCPM.txt"),
                              header=TRUE, sep="\t", row.names=1, check.names=FALSE))
cand0 <- readLines(file.path(ym,"17.交集基因/interGenes.txt"))[-1]
deAll <- read.delim(file.path(ym,"09.差异分析/DE.PI_vs_H.all.txt"), check.names=FALSE)
deAll <- deAll[deAll$gene %in% cand0, ]
cand <- head(deAll$gene[order(-abs(deAll$t))], 50)
clinical <- read.csv(file.path(ym,"05.bulk标准化质控/clinical.csv"))
X <- t(lcpm0[cand, clinical$sample])
y <- clinical$group; y <- factor(ifelse(y=="Periimplantitis","PI","H"), levels=c("H","PI"))
idxH <- which(y=="H"); idxP <- which(y=="PI")
cat("矩阵:", nrow(X),"x",ncol(X),"\n")

# SVM-RFE 排序函数
rfe_rank <- function(X,y){
  feats<-colnames(X); ro<-c()
  while(length(feats)>1){
    sv<-svm(X[,feats,drop=FALSE],y,kernel="linear",cost=1,scale=FALSE)
    w<-t(sv$coefs)%*%sv$SV; w2<-colSums(w^2)
    worst<-which.min(w2); ro<-c(feats[worst],ro); feats<-feats[-worst]
  }
  c(feats,ro)
}

B <- 100
freqL <- setNames(numeric(ncol(X)), colnames(X))
freqS <- setNames(numeric(ncol(X)), colnames(X))
freqR <- setNames(numeric(ncol(X)), colnames(X))
set.seed(2026)
for(b in 1:B){
  tr <- c(sample(idxH,8), sample(idxP,8))     # 分层抽 16（80%）
  # LASSO
  cv <- cv.glmnet(X[tr,], y[tr], family="binomial", alpha=1, nfolds=4)
  cf <- coef(cv, s="lambda.min")
  sel <- setdiff(rownames(cf)[cf[,1]!=0],"(Intercept)")
  freqL[sel] <- freqL[sel]+1
  # SVM-RFE（记 top5）
  r <- rfe_rank(X[tr,,drop=FALSE], y[tr])
  ss <- r[1:min(5,length(r))]
  freqS[ss] <- freqS[ss]+1
  # RF（记 top10 by MeanDecreaseAccuracy）
  rf <- randomForest(X[tr,], y[tr], ntree=500, importance=TRUE)
  im_raw <- importance(rf, type=1)
  if(is.matrix(im_raw)){
    im <- if("MeanDecreaseAccuracy" %in% colnames(im_raw)) im_raw[,"MeanDecreaseAccuracy"] else im_raw[,ncol(im_raw)]
  } else im <- setNames(as.numeric(im_raw), names(im_raw))
  tt <- names(sort(im, decreasing=TRUE))[1:10]
  freqR[tt] <- freqR[tt]+1
  if(b %% 20 == 0) cat("完成", b, "/", B, "\n")
}
fL <- sort(freqL, decreasing=TRUE)
fS <- sort(freqS, decreasing=TRUE)
fR <- sort(freqR, decreasing=TRUE)
write.table(data.frame(gene=names(fL),freq=fL/B), "LASSO_freq.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(gene=names(fS),freq=fS/B), file.path(ym,"21.SVM/SVM_freq.txt"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(gene=names(fR),freq=fR/B), file.path(ym,"22.随机森林/RF_freq.txt"), sep="\t", quote=FALSE, row.names=FALSE)

thr <- 0.5
gL <- names(fL)[fL/B >= thr]
gS <- names(fS)[fS/B >= thr]
gR <- names(fR)[fR/B >= thr]
writeLines(gL, "LASSO_genes.txt")
writeLines(gS, file.path(ym,"21.SVM/SVM_genes.txt"))
writeLines(gR, file.path(ym,"22.随机森林/RF_genes.txt"))
cat("\n稳定基因(频率>=50%): LASSO", length(gL)," SVM",length(gS)," RF",length(gR),"\n")
cat("LASSO:", paste(gL,collapse=" "),"\n")
cat("SVM:", paste(gS,collapse=" "),"\n")
cat("RF:", paste(gR,collapse=" "),"\n")

core <- Reduce(intersect, list(gL,gS,gR))
cat("\n★ 核心基因（三方法稳定交集）:", length(core), "->", paste(core,collapse=" "), "\n")
writeLines(core, file.path(ym,"23.交集特征基因/core_genes.txt"))
sets <- list(LASSO=gL, SVM_RFE=gS, RF=gR)
pv <- ggvenn(sets, fill_color=c("#E64B35","#4DBBD5","#00A087"),
             stroke_size=.6, set_name_size=4, text_size=3.5)
ggsave(file.path(ym,"23.交集特征基因/feature_venn.pdf"), pv, width=7, height=6)
saveRDS(list(X=X,y=y,freqL=fL,freqS=fS,freqR=fR,gL=gL,gS=gS,gR=gR,core=core),
        file.path(ym,"23.交集特征基因/ml.rds"))
cat("已输出稳定性选择结果\n")
