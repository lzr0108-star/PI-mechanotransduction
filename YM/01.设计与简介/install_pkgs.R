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
# YM项目 R包安装（bulk部分）—— 用户级库

options(repos=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(install.packages.check.source="no")

# BiocManager
if(!requireNamespace("BiocManager", quietly=TRUE)){
  install.packages("BiocManager", quiet=TRUE)
}
# Bioconductor 包
BiocManager::install(c("limma","edgeR"), ask=FALSE, update=FALSE, quiet=TRUE)
# CRAN 包
cran.pkgs <- c("pheatmap","ggplot2","ggrepel","reshape2","ggpubr","RColorBrewer","dplyr","cowplot")
to.install <- cran.pkgs[!sapply(cran.pkgs, function(p) requireNamespace(p, quietly=TRUE))]
if(length(to.install)>0) install.packages(to.install, quiet=TRUE)

# 验证加载
need <- c("limma","edgeR","pheatmap","ggplot2","ggrepel","reshape2","ggpubr","RColorBrewer","dplyr")
ok <- sapply(need, function(p) requireNamespace(p, quietly=TRUE))
print(ok)
if(all(ok)) cat("ALL_PACKAGES_OK\n") else cat("MISSING:", paste(names(ok)[!ok], collapse=","), "\n")
