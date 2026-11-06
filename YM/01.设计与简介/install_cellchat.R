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
options(repos=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN"),
        BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# CRAN 依赖
cran.dep <- c("NMF","ggalluvial","remotes","patchwork","ggrepel","RColorBrewer")
to.in <- cran.dep[!sapply(cran.dep, function(p) requireNamespace(p, quietly=TRUE))]
if(length(to.in)>0) install.packages(to.in, quiet=TRUE)
# Bioc 依赖
if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager", quiet=TRUE)
bioc.dep <- c("ComplexHeatmap","BiocNeighbors")
to.bioc <- bioc.dep[!sapply(bioc.dep, function(p) requireNamespace(p, quietly=TRUE))]
if(length(to.bioc)>0) BiocManager::install(to.bioc, ask=FALSE, update=FALSE, quiet=TRUE)
# CellChat via GitHub proxy
if(!requireNamespace("CellChat", quietly=TRUE)){
  remotes::install_git("https://ghfast.top/https://github.com/sqjin/CellChat",
                       upgrade="never", quiet=TRUE)
}
cat("CellChat:", requireNamespace("CellChat", quietly=TRUE), "\n")
