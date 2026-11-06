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
library(Seurat)
mk <- function(sub){
  sd <- file.path(dataDir, sub)
  sname <- sub
  mat <- Read10X(sd)
  genes <- read.delim(file.path(sd,"genes.tsv"), header=FALSE)
  sym <- genes$V2[match(rownames(mat),genes$V1)]
  sym[is.na(sym)] <- rownames(mat)[is.na(sym)]
  rownames(mat) <- make.unique(sym)
  seu <- CreateSeuratObject(mat, min.cells=3, min.features=200)
  seu <- RenameCells(seu, add.cell.id=sname)
  seu$sample <- sname
  seu$condition <- if(grepl("_N_",sname)) "Healthy" else "Periimplantitis"
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
  seu
}
a <- mk("GSM8410650_N_1")
p <- mk("GSM8410653_PI_1")
cat("a:", ncol(a), " p:", ncol(p), "\n")
cat("a gene NA:", sum(is.na(rownames(a))), " p gene NA:", sum(is.na(rownames(p))), "\n")
m <- merge(a, p)
cat("merge OK:", ncol(m), "\n")
