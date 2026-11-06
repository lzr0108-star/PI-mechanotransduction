# ===========================================================================
# Project configuration
# Reproducible code for the peri-implantitis mechanotransduction study.
# This file is located automatically and sourced by every analysis script,
# so no hard-coded absolute paths are needed.
# ===========================================================================

# ---- Locate the repository root (it contains 'YM/' and this file) --------
.find_repo_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  for (i in seq_len(12)) {
    if (dir.exists(file.path(d, "YM")) &&
        file.exists(file.path(d, "config.R"))) {
      return(d)
    }
    up <- dirname(d)
    if (identical(up, d)) break
    d <- up
  }
  stop("Repository root not found. Please run from within the repository.")
}

# Standard directory handles used across the analysis scripts.
root <- .find_repo_root()
ym   <- file.path(root, "YM")
base <- ym                                   # alias used by bulk/QC scripts

# ---- External curated / third-party inputs (see README for download) -----
mrg_file    <- file.path(root, "17.venn", "MRG.txt")
kegg_gmt    <- file.path(root, "27.GSEA", "c2.cp.kegg.Hs.symbols.gmt")
dsigdb_file <- file.path(root, "39.drugEnrich", "DSigDB_All_detailed.txt")

# ---- Convenience helper: ensure a directory exists -----------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

cat("Configuration loaded. Repository root:", root, "\n")
