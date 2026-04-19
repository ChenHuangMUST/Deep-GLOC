setwd("D:/phd cases/metabolic analysis/code/Deep_GLOC/getG_liver")
{library(data.table)
  library(edgeR)
  library(dplyr)
}
infile <- "D:/phd cases/metabolic analysis/data/GTEx/gene_reads_v11_liver.gct"
outdir <- "D:/phd cases/metabolic analysis/code/Deep_GLOC/getG_liver"
gtex <- fread(infile, skip = 2, data.table = FALSE)
cat("原始数据维度: ", dim(gtex), "\n")
print(colnames(gtex)[1:6])

# 前两列
gene_id <- gtex[, 1]
gene_symbol <- gtex[, 2]

# 表达矩阵
expr <- gtex[, -(1:2), drop = FALSE]
expr <- as.matrix(expr)
mode(expr) <- "numeric"
rownames(expr) <- gene_id

cat("表达矩阵维度: ", dim(expr), "\n")
cat("样本数: ", ncol(expr), "\n")

# -----------------------------
# 4. 去除重复基因和全零基因
# -----------------------------
# 对于重复 gene_id，保留总表达量最大的那一条
if(any(duplicated(rownames(expr)))) {
  rs <- rowSums(expr, na.rm = TRUE)
  
  ord <- order(rownames(expr), -rs)
  expr <- expr[ord, , drop = FALSE]
  gene_id <- gene_id[ord]
  gene_symbol <- gene_symbol[ord]
  
  keep <- !duplicated(rownames(expr))
  expr <- expr[keep, , drop = FALSE]
  gene_id <- gene_id[keep]
  gene_symbol <- gene_symbol[keep]
}

# 去掉全0基因
keep_nonzero <- rowSums(expr) > 0
expr <- expr[keep_nonzero, , drop = FALSE]
gene_id <- gene_id[keep_nonzero]
gene_symbol <- gene_symbol[keep_nonzero]

# -----------------------------
# 5. 计算 CPM
# -----------------------------
dge <- DGEList(counts = expr)
cpm_mat <- cpm(dge)

# -----------------------------
# 6. 定义 G_liver
# CPM > 1 in at least 20% of samples
# -----------------------------
min_prop <- 0.20
min_samples <- ceiling(ncol(cpm_mat) * min_prop)

keep_expr <- rowSums(cpm_mat > 1) >= min_samples
G_liver <- rownames(cpm_mat)[keep_expr]

cat("CPM阈值: > 1\n")
cat("最少样本数: ", min_samples, " / ", ncol(cpm_mat), "\n")
cat("G_liver基因数: ", length(G_liver), "\n")

# 汇总表
res_df <- data.frame(
  gene_id = rownames(cpm_mat),
  gene_symbol = gene_symbol[match(rownames(cpm_mat), gene_id)],
  mean_count = rowMeans(expr),
  median_count = apply(expr, 1, median),
  mean_cpm = rowMeans(cpm_mat),
  median_cpm = apply(cpm_mat, 1, median),
  prop_cpm_gt1 = rowMeans(cpm_mat > 1),
  in_G_liver = keep_expr,
  stringsAsFactors = FALSE
)

res_df <- res_df %>%
  arrange(desc(in_G_liver), desc(mean_cpm))

# -----------------------------
# 7. 加入蛋白编码注释
# -----------------------------
{library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)}

edb <- EnsDb.Hsapiens.v86
gene_anno <- genes(edb, return.type = "data.frame")

gene_anno2 <- unique(gene_anno[, c("gene_id", "gene_biotype")])
colnames(gene_anno2) <- c("gene_id", "GENETYPE")

# 去版本号
res_df$gene_id <- sub("\\..*$", "", res_df$gene_id)

# 保留顺序
res_df$gene_order <- seq_len(nrow(res_df))

# 合并注释
res_df1 <- merge(res_df, gene_anno2, by = "gene_id", all.x = TRUE)
res_df1 <- res_df1[order(res_df1$gene_order), ]
res_df1$gene_order <- NULL

# 全部 G_liver
G_liver_all_df <- res_df1[res_df1$in_G_liver, ]
G_liver_all <- G_liver_all_df$gene_id
length(G_liver_all)
library(biomaRt)
# 1. 选择ensembl数据库
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# 2. NA基因列表（示例）
na_genes <- G_liver_all_df$gene_symbol[which(is.na(G_liver_all_df$GENETYPE))] # 这里替换成你全部NA基因
# 3. 查询gene_biotype
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  filters = c( "hgnc_symbol"),
  values = na_genes,
  mart = ensembl
)
protein_coding_genes <- gene_info[gene_info$gene_biotype == "protein_coding", ]
head(protein_coding_genes)
res_df1$GENETYPE[which(res_df1$gene_symbol%in%protein_coding_genes$hgnc_symbol)] <- "protein_coding"
# protein-coding G_liver
G_liver_pc_df <- res_df1[res_df1$in_G_liver & res_df1$GENETYPE == "protein_coding", ]
G_liver_pc_df <- na.omit(G_liver_pc_df)
G_liver_pc <- G_liver_pc_df$gene_id

cat("全部 G_liver 基因数: ", length(G_liver_all), "\n")
cat("protein-coding G_liver 基因数: ", length(G_liver_pc), "\n")

# -----------------------------
# 8. 保存结果
# -----------------------------
# 所有基因汇总
write.table(
  res_df1,
  file = file.path(outdir, "GTEx_liver_gene_summary_with_biotype.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 全部 G_liver
write.table(
  data.frame(gene_id = G_liver_all),
  file = file.path(outdir, "G_liver_gene_ids.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  G_liver_all_df[, c("gene_id", "gene_symbol", "mean_cpm", "prop_cpm_gt1", "GENETYPE")],
  file = file.path(outdir, "G_liver_genes_with_symbol.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# protein-coding G_liver
write.table(
  data.frame(gene_id = G_liver_pc),
  file = file.path(outdir, "G_liver_protein_coding_gene_ids.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  G_liver_pc_df[, c("gene_id", "gene_symbol", "mean_cpm", "prop_cpm_gt1", "GENETYPE")],
  file = file.path(outdir, "G_liver_protein_coding_with_symbol.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------------
# 9. 保存 CPM 表达矩阵
# -----------------------------

# 保证 cpm_mat 行名为无版本号 ENSG
rownames(cpm_mat) <- sub("\\..*$", "", rownames(cpm_mat))
anyDuplicated(rownames(cpm_mat))
#cpm_mat <- aggregate(cpm_mat,by = list(gene_id = rownames(cpm_mat)),FUN = max)

# 1) 全部基因的 CPM 表达矩阵
all_cpm_df <- data.frame(
  gene_id = rownames(cpm_mat),
  gene_symbol = res_df1$gene_symbol[match(rownames(cpm_mat), res_df1$gene_id)],
  cpm_mat,
  check.names = FALSE
)

write.table(
  all_cpm_df,
  file = file.path(outdir, "GTEx_liver_all_genes_CPM.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 2) G_liver 的 CPM 表达矩阵
G_liver_all_ids <- G_liver_all_df$gene_id
G_liver_cpm_mat <- cpm_mat[rownames(cpm_mat) %in% G_liver_all_ids, , drop = FALSE]

G_liver_cpm_df <- data.frame(
  gene_id = rownames(G_liver_cpm_mat),
  gene_symbol = G_liver_all_df$gene_symbol[match(rownames(G_liver_cpm_mat), G_liver_all_df$gene_id)],
  G_liver_cpm_mat,
  check.names = FALSE
)

write.table(
  G_liver_cpm_df,
  file = file.path(outdir, "GTEx_liver_G_liver_CPM.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 3) 全部 protein-coding 基因的 CPM 表达矩阵（可选）
all_pc_ids <- res_df1$gene_id[res_df1$GENETYPE == "protein_coding"]
all_pc_cpm_mat <- cpm_mat[rownames(cpm_mat) %in% all_pc_ids, , drop = FALSE]

all_pc_cpm_df <- data.frame(
  gene_id = rownames(all_pc_cpm_mat),
  gene_symbol = res_df1$gene_symbol[match(rownames(all_pc_cpm_mat), res_df1$gene_id)],
  all_pc_cpm_mat,
  check.names = FALSE
)

write.table(
  all_pc_cpm_df,
  file = file.path(outdir, "GTEx_liver_all_protein_coding_CPM.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 4) G_liver 的 protein-coding CPM 表达矩阵
G_liver_pc_ids <- G_liver_pc_df$gene_id
G_liver_pc_cpm_mat <- cpm_mat[rownames(cpm_mat) %in% G_liver_pc_ids, , drop = FALSE]

G_liver_pc_cpm_df <- data.frame(
  gene_id = rownames(G_liver_pc_cpm_mat),
  gene_symbol = G_liver_pc_df$gene_symbol[match(rownames(G_liver_pc_cpm_mat), G_liver_pc_df$gene_id)],
  G_liver_pc_cpm_mat,
  check.names = FALSE
)

write.table(
  G_liver_pc_cpm_df,
  file = file.path(outdir, "GTEx_liver_G_liver_protein_coding_CPM.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

save.image("code_Gliver.RData")
