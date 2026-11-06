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
library(ConsensusClusterPlus)
clinical <- read.csv(file.path(ym,"05.bulk标准化质控/clinical.csv"))
lcpm <- as.matrix(read.table(file.path(ym,"05.bulk标准化质控/expr.logCPM.txt"),
                             header=TRUE, sep="\t", row.names=1, check.names=FALSE))
lcpm <- lcpm[, clinical$sample]
madg <- apply(lcpm,1,mad)
topG <- names(sort(madg,decreasing=TRUE))[1:5000]
dat <- t(lcpm[topG,])                 # 20样本 x 5000基因
cat("dat:", nrow(dat),"x",ncol(dat),"  rownames head:",paste(head(rownames(dat)),collapse=" "),"\n")
dat2 <- apply(dat,2,function(z)(z-mean(z))/sd(z))
rownames(dat2) <- rownames(dat)
cat("dat2:", nrow(dat2),"x",ncol(dat2),"\n")
res <- ConsensusClusterPlus(dat2, maxK=6, reps=50, pItem=.8, pFeature=1,
                            clusterAlg="km", distance="euclidean",
                            title=file.path(getwd(),"CCP_debug"), plot="pdf", seed=2026)
cat("\n--- str res[[2]] ---\n")
str(res[[2]], max.level=1)
cat("\nconsensusClass:\n"); print(res[[2]]$consensusClass)
cat("\ntable:\n"); print(table(res[[2]]$consensusClass))
