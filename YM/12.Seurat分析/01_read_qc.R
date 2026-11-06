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
# 12-① 读取 8 样本 + QC 过滤
library(Seurat); library(ggplot2); library(cowplot)

dataDir <- file.path(root,"YM/10.单细胞数据下载/GSE272774")
setwd(file.path(root,"YM/12.Seurat分析"))

mergeList <- function(lst){
  m <- lst[[1]]
  for(i in 2:length(lst)) m <- merge(m, lst[[i]])
  m
}

sampleDirs <- list.dirs(dataDir, recursive=FALSE)
sampleDirs <- sampleDirs[grepl("GSM", sampleDirs)]
seuList <- list()
for(sd in sampleDirs){
  sname <- basename(sd)
  cond <- if(grepl("_N_", sname)) "Healthy" else "Periimplantitis"
  mat <- Read10X(sd)
  genes <- read.delim(file.path(sd,"genes.tsv"), header=FALSE)
  sym <- genes$V2[match(rownames(mat), genes$V1)]
  sym[is.na(sym)] <- rownames(mat)[is.na(sym)]
  rownames(mat) <- make.unique(sym)
  seu <- CreateSeuratObject(mat, project=sname, min.cells=3, min.features=200)
  seu <- RenameCells(seu, add.cell.id=sname)
  seu$sample <- sname
  seu$condition <- cond
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
  seuList[[sname]] <- seu
  cat(sname, cond, "细胞:", ncol(seu), "\n")
}

before <- mergeList(seuList)
p1 <- VlnPlot(before, c("nFeature_RNA","nCount_RNA","percent.mt"), group.by="sample", pt.size=.1)
ggsave("QC_before.pdf", p1, width=13, height=5)

seuList <- lapply(seuList, function(seu){
  subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 &
                   nCount_RNA > 500 & percent.mt < 20)
})
after <- mergeList(seuList)
cat("QC 后总细胞数:", ncol(after), "\n")
print(table(after$condition))
p2 <- VlnPlot(after, c("nFeature_RNA","nCount_RNA","percent.mt"), group.by="sample", pt.size=.1)
ggsave("QC_after.pdf", p2, width=13, height=5)
saveRDS(seuList, "seuList_qc.rds")
cat("已保存 seuList_qc.rds\n")
