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
if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",quiet=TRUE)
if(!requireNamespace("monocle3",quietly=TRUE))
  BiocManager::install("monocle3", ask=FALSE, update=FALSE, quiet=TRUE)
cat("monocle3:", requireNamespace("monocle3",quietly=TRUE), "\n")
