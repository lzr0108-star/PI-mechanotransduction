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
# 03 解析核对 GSE223924 raw count

cnt <- read.csv("GSE223924_raw_count.csv", row.names=1, check.names=FALSE)
cnt <- as.matrix(cnt)

cat("==== 1. 维度 ====\n")
cat("基因数:", nrow(cnt), " 样本数:", ncol(cnt), "\n")

cat("\n==== 2. 样本分组（按前缀 H/PI/PT）====\n")
prefix <- sub("[0-9]+$", "", colnames(cnt))
num    <- sub("^[A-Z]+", "", colnames(cnt))
grp <- c(H="Healthy", PI="Peri-implantitis", PT="Periodontitis")[prefix]
print(table(grp))

cat("\n==== 3. 配对关系（PT 与 PI 同患者，编号相同）====\n")
pi.n <- sort(as.integer(num[prefix=="PI"]))
pt.n <- sort(as.integer(num[prefix=="PT"]))
cat("PI 编号:", paste(pi.n, collapse=","), "\n")
cat("PT 编号:", paste(pt.n, collapse=","), "\n")
cat("PT/PI 编号完全对应(可配对):", identical(pi.n, pt.n), "\n")

cat("\n==== 4. 数据类型/缺失/重复 ====\n")
cat("是否全为非负整数:", all(cnt >= 0) && all(cnt == floor(cnt)), "\n")
cat("缺失值 NA 数:", sum(is.na(cnt)), "\n")
cat("重复基因名:", sum(duplicated(rownames(cnt))), "\n")
cat("重复样本名:", sum(duplicated(colnames(cnt))), "\n")

cat("\n==== 5. 文库大小（每样本总count, 百万）====\n")
lib <- colSums(cnt)
cat("整体文库(百万):\n"); print(round(summary(lib)/1e6, 2))
cat("各组平均文库(百万):\n"); print(round(tapply(lib, grp, mean)/1e6, 2))

cat("\n==== 6. 零计数基因 ====\n")
allzero <- rowSums(cnt) == 0
cat("全样本零表达基因数:", sum(allzero), "\n")
detected <- colSums(cnt > 0)
cat("每样本检出基因数: 最小", min(detected), " 最大", max(detected), "\n")
cat("各组平均检出基因数:\n"); print(round(tapply(detected, grp, mean), 0))

cat("\n==== 7. 样本列名（核对顺序）====\n")
cat(paste(colnames(cnt), collapse=" "), "\n")
