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
options(repos=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
# 依赖
cran.dep <- c("ggforce","scatterpie","ggraph","ggalluvial","ggrepel","ggdendro","cowplot")
to.in <- cran.dep[!sapply(cran.dep, function(p) requireNamespace(p, quietly=TRUE))]
if(length(to.in)>0) install.packages(to.in, quiet=TRUE)
# ggunchull
if(!requireNamespace("ggunchull", quietly=TRUE))
  remotes::install_git("https://ghfast.top/https://github.com/sajuukLyu/ggunchull", upgrade="never", quiet=TRUE)
# scRNAtoolVis
if(!requireNamespace("scRNAtoolVis", quietly=TRUE))
  remotes::install_git("https://ghfast.top/https://github.com/junjunlab/scRNAtoolVis", upgrade="never", quiet=TRUE)
# ClusterGVis
if(!requireNamespace("ClusterGVis", quietly=TRUE))
  remotes::install_git("https://ghfast.top/https://github.com/junjunlab/ClusterGVis", upgrade="never", quiet=TRUE)
cat("ggunchull:", requireNamespace("ggunchull",quietly=TRUE),
    "scRNAtoolVis:", requireNamespace("scRNAtoolVis",quietly=TRUE),
    "ClusterGVis:", requireNamespace("ClusterGVis",quietly=TRUE), "\n")
