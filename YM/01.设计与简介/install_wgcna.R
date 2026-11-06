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
# 安装 WGCNA
options(repos=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# Bioconductor 依赖
bioc <- c("impute","GO.db","preprocessCore","AnnotationDbi")
BiocManager::install(bioc, ask=FALSE, update=FALSE, quiet=TRUE)

# WGCNA（CRAN，自动装 CRAN 依赖）
if(!requireNamespace("WGCNA", quietly=TRUE)){
  install.packages("WGCNA", quiet=TRUE)
}
cat("WGCNA:", requireNamespace("WGCNA", quietly=TRUE), "\n")
if(requireNamespace("WGCNA", quietly=TRUE)) cat("WGCNA_OK\n")
