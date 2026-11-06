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
# 机器学习相关 R 包
options(repos=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN"))

pkgs <- c("glmnet","e1071","randomForest","pROC","caret","Matrix","foreach","doParallel")
miss <- pkgs[!sapply(pkgs, function(p) requireNamespace(p, quietly=TRUE))]
if(length(miss)>0) install.packages(miss, quiet=TRUE)
ok <- sapply(pkgs, function(p) requireNamespace(p, quietly=TRUE))
print(ok)
if(all(ok)) cat("ML_OK\n") else cat("MISSING:", paste(names(ok)[!ok], collapse=","), "\n")
