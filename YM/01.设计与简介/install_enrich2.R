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
# 补装富集相关包（脚本文件，避免命令行截断）
options(repos=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

need <- c("clusterProfiler","org.Hs.eg.db","circlize","enrichplot","DOSE")
cat("安装前:\n"); print(sapply(need, function(p) requireNamespace(p, quietly=TRUE)))
miss <- need[!sapply(need, function(p) requireNamespace(p, quietly=TRUE))]
if(length(miss) > 0){
  if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", quiet=TRUE)
  BiocManager::install(miss, ask=FALSE, update=FALSE, quiet=TRUE)
}
cat("\n安装后:\n"); print(sapply(need, function(p) requireNamespace(p, quietly=TRUE)))
if(all(sapply(need, function(p) requireNamespace(p, quietly=TRUE)))) cat("ENRICH_OK\n")
