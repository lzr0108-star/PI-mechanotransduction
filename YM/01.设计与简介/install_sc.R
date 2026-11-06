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
# 安装单细胞核心 R 包（CRAN + Bioconductor；不含 github 包）
options(repos=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(install.packages.check.source="no")

# CRAN
cran <- c("Seurat","SeuratObject","harmony","remotes","BiocManager","patchwork","cowplot","Matrix")
ic <- cran[!sapply(cran, function(p) requireNamespace(p, quietly=TRUE))]
if(length(ic)>0) install.packages(ic, quiet=TRUE)

# Bioconductor
bioc <- c("SingleR","celldex","BiocParallel","zellkonverter","scater","scran","SummarizedExperiment","SingleCellExperiment")
BiocManager::install(bioc, ask=FALSE, update=FALSE, quiet=TRUE)

# 验证
need <- unique(c(cran, bioc))
ok <- sapply(need, function(p) requireNamespace(p, quietly=TRUE))
print(ok)
if(all(ok)) cat("ALL_SC_OK\n") else cat("MISSING:", paste(names(ok)[!ok], collapse=","), "\n")
