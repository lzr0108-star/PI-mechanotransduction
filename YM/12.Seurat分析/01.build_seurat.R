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
# 12 (1) 从 h5ad 构建 Seurat（reticulate，raw count）
library(reticulate); library(Seurat); library(Matrix)

h5 <- file.path(base, "10.单细胞数据下载/periodontal_atlas.h5ad")
outdir <- file.path(base, "12.Seurat分析"); setwd(outdir)

np <- import("numpy", convert=TRUE)
ad <- import("anndata", convert=FALSE)
a <- ad$read_h5ad(h5)

# ---- raw count: 转成 genes x cells 的 csc ----
csc <- a$raw$X$T$tocsc()
csc$sort_indices()
idx <- as.integer(py_to_r(csc$indices))
ptr <- as.integer(py_to_r(csc$indptr))
xdt <- as.numeric(py_to_r(csc$data))
shp <- as.integer(py_to_r(csc$shape))
dg <- new("dgCMatrix", i=idx, p=ptr, x=xdt, Dim=shp)
cat("count 矩阵 (genes x cells):", nrow(dg), "x", ncol(dg), "\n")

# ---- 基因名（feature_name=symbol；NaN 用 ENSG；去重）----
sym <- as.character(py_to_r(a$raw$var["feature_name"]$values$tolist()))
ensg <- as.character(py_to_r(a$raw$var$index$tolist()))
sym[is.na(sym) | sym==""] <- ensg[is.na(sym) | sym==""]
if(any(duplicated(sym))) sym <- make.unique(sym)
rownames(dg) <- sym

# ---- 细胞名 ----
bcs <- as.character(py_to_r(a$obs_names$tolist()))
colnames(dg) <- bcs

# ---- metadata（obs；python端选列并转字符串，避免nullable-int转换bug）----
keepcols <- c("cell_type","disease","Sample ID","donor_id","Study","GSE","sex","tissue",
              "Annotated Cell Sets","Tooth Type","Anterior vs Posterior")
obs_py <- a$obs[r_to_py(as.list(keepcols))]$astype("str")
meta <- as.data.frame(py_to_r(obs_py))
rownames(meta) <- bcs
for(cc in colnames(meta)) meta[[cc]] <- as.character(meta[[cc]])

# ---- 建 Seurat ----
srt <- CreateSeuratObject(dg, min.cells=0, min.features=0, meta.data=meta)
srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern="^MT-")
srt[["percent.ribo"]] <- PercentageFeatureSet(srt, pattern="^RP[SL]")

cat("\n各组(disease)细胞数:\n"); print(table(srt$disease))
cat("\n各细胞类型(cell_type)细胞数:\n"); print(table(srt$cell_type))
cat("\n各Study细胞数:\n"); print(table(srt$Study))

saveRDS(srt, "seurat_ref.rds")
cat("\n已保存 seurat_ref.rds\n")
