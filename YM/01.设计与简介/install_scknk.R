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
if(!requireNamespace("remotes",quietly=TRUE)) install.packages("remotes",quiet=TRUE)
if(!requireNamespace("scTenifoldKnk",quietly=TRUE))
  remotes::install_git("https://ghfast.top/https://github.com/cailab-tamu/scTenifoldKnk",
                       upgrade="never", quiet=TRUE)
cat("scTenifoldKnk:", requireNamespace("scTenifoldKnk",quietly=TRUE), "\n")
