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
library(randomForest)
lcpm0 <- as.matrix(read.table(file.path(ym,"05.bulk标准化质控/expr.logCPM.txt"), header=TRUE, sep="\t", row.names=1, check.names=FALSE))
cand0 <- readLines(file.path(ym,"17.交集基因/interGenes.txt"))[-1]
deAll <- read.delim(file.path(ym,"09.差异分析/DE.PI_vs_H.all.txt"), check.names=FALSE)
deAll <- deAll[deAll$gene %in% cand0,]
cand <- head(deAll$gene[order(-abs(deAll$t))],50)
clinical <- read.csv(file.path(ym,"05.bulk标准化质控/clinical.csv"))
X <- t(lcpm0[cand, clinical$sample])
y <- factor(ifelse(clinical$group=="Periimplantitis","PI","H"), levels=c("H","PI"))
idxH<-which(y=="H"); idxP<-which(y=="PI")
set.seed(2026)
tr <- c(sample(idxH,8), sample(idxP,8))
rf <- randomForest(X[tr,], y[tr], ntree=500)
im_raw <- importance(rf, type=1)
cat("class:", class(im_raw), " is.matrix:", is.matrix(im_raw),"\n")
if(is.matrix(im_raw)){ cat("dim:",dim(im_raw),"\n"); cat("colnames:",paste(colnames(im_raw),collapse=" | "),"\n") }
cat("head names:", paste(head(rownames(im_raw)),collapse=" "),"\n")
if(is.matrix(im_raw)){
  im <- if("MeanDecreaseAccuracy" %in% colnames(im_raw)) im_raw[,"MeanDecreaseAccuracy"] else im_raw[,ncol(im_raw)]
} else im <- im_raw
cat("length im:", length(im), " names head:", paste(head(names(im)),collapse=" "),"\n")
tt <- names(sort(im, decreasing=TRUE))[1:10]
cat("tt:", paste(tt, collapse=" "),"\n")
cat("match freq names:", sum(tt %in% colnames(X)), "/10\n")
