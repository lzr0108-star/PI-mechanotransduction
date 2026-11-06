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
# 修复 13 章 01 / 03（as.ggplot + patchwork 并排）
library(CellChat); library(patchwork); library(cowplot); library(ggplotify)

setwd(file.path(root,"YM/13.细胞通讯"))
cc <- readRDS("cellchat_compare.rds")
countH <- cc@net$Healthy$count
countP <- cc@net$PI$count
gsH <- as.vector(table(cc@idents$Healthy))
gsP <- as.vector(table(cc@idents$PI))
maxCount <- max(c(countH, countP))

# ===== 01：两组并排 =====
g1 <- as.ggplot(~netVisual_circle(countH, vertex.weight=gsH, weight.scale=TRUE,
                                 edge.weight.max=maxCount, edge.width.max=7,
                                 alpha.edge=.42, vertex.label.cex=1.05,
                                 title.name="Healthy: interactions"))
g2 <- as.ggplot(~netVisual_circle(countP, vertex.weight=gsP, weight.scale=TRUE,
                                 edge.weight.max=maxCount, edge.width.max=7,
                                 alpha.edge=.42, vertex.label.cex=1.05,
                                 title.name="Periimplantitis: interactions"))
ggsave("01_network_count.pdf", g1+g2, width=15, height=7.4)

# ===== 03：手动差异，增强/减弱并排 =====
types <- intersect(rownames(countH), rownames(countP))
diff <- countP[types, types] - countH[types, types]
up <- diff; up[up < 0] <- 0
down <- -diff; down[down < 0] <- 0
maxD <- max(c(up, down))
u1 <- as.ggplot(~netVisual_circle(up, vertex.weight=gsP, weight.scale=TRUE,
                                 edge.weight.max=maxD, edge.width.max=7,
                                 alpha.edge=.5, vertex.label.cex=1.05,
                                 title.name="Enhanced in PI"))
u2 <- as.ggplot(~netVisual_circle(down, vertex.weight=gsH, weight.scale=TRUE,
                                 edge.weight.max=maxD, edge.width.max=7,
                                 alpha.edge=.5, vertex.label.cex=1.05,
                                 title.name="Reduced in PI"))
ggsave("03_diff_interaction.pdf", u1+u2, width=15, height=7.4)
cat("01 / 03 已重新生成\n")
