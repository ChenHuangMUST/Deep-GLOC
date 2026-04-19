setwd("D:/phd cases/metabolic analysis/code/Deep_GLOC/get_meta_graph_new")
G_liver_file <- "D:/phd cases/metabolic analysis/code/Deep_GLOC/getG_liver/G_liver_protein_coding_with_symbol.txt"
G_liver <- read.table(G_liver_file, header = TRUE, sep = "\t", check.names = FALSE)
G_liver <- na.omit(G_liver)
G_liver_genes <- unique(G_liver$gene_id)
length(G_liver_genes)
##读取数据集的表达数据
read_logCPM <- function(gse){
  file <- paste0(
    "D:/phd cases/metabolic analysis/data/GEO/fatty_liver/human/exp_human_CPM/",
    gse,
    "_logCPM.txt"
  )
  
  dat <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE)
  rownames(dat) <- dat$gene_id
  dat <- dat[, -(1:2), drop = FALSE]   # 去掉 gene_id 和 gene_symbol
  
  expr <- as.matrix(dat)
  mode(expr) <- "numeric"
  
  return(expr)
}
#G_liver 过滤
filter_G_liver <- function(expr, G_liver_genes){
  common_genes <- intersect(rownames(expr), G_liver_genes)
  expr <- expr[common_genes, , drop = FALSE]
  
  cat("Genes after G_liver filter:", nrow(expr), "\n")
  return(expr)
}
#variance filter 保留 top 50% variance genes
filter_variance <- function(expr, prop = 0.5){
  vars <- apply(expr, 1, var, na.rm = TRUE)
  cutoff <- quantile(vars, probs = 1 - prop, na.rm = TRUE)
  expr <- expr[vars > cutoff, , drop = FALSE]
  
  cat("Genes after variance filter:", nrow(expr), "\n")
  return(expr)
}

#构建单个数据集网络输入矩阵
prepare_dataset_expr <- function(gse, G_liver_genes, var_prop = 0.5){
  expr <- read_logCPM(gse)
  cat("Dataset:", gse, "\n")
  cat("Original genes:", nrow(expr), " Samples:", ncol(expr), "\n")
  
  expr <- filter_G_liver(expr, G_liver_genes)
  expr <- filter_variance(expr, prop = var_prop)
  
  return(expr)
}

datasets <- c("GSE126848","GSE167523","GSE174478",
              "GSE193066","GSE193080","GSE225740")
###批量读取并过滤
expr_list <- lapply(datasets, prepare_dataset_expr, G_liver_genes = G_liver_genes, var_prop = 0.5)
names(expr_list) <- datasets
###看每个数据集保留下来的基因数和样本数
sapply(expr_list, nrow)
sapply(expr_list, ncol)
###取交集基因
genes_meta <- Reduce(intersect, lapply(expr_list, rownames))
length(genes_meta)
expr_list2 <- lapply(expr_list, function(x){
  x[genes_meta, , drop = FALSE]
})
##确实是否全为TRUE
all(sapply(expr_list2, function(x) identical(rownames(x), genes_meta)))

###构建每个队列的相关矩阵
build_cor_mat <- function(expr){
  cor_mat <- cor(
    t(expr),
    method = "pearson",
    use = "pairwise.complete.obs"
  )
  return(cor_mat)
}
cor_list <- lapply(expr_list2, build_cor_mat)

##Fisher z 转换
fisher_z <- function(cor_mat, eps = 1e-7){
  cor_mat[cor_mat >=  1] <-  1 - eps
  cor_mat[cor_mat <= -1] <- -1 + eps
  z_mat <- atanh(cor_mat)
  diag(z_mat) <- 0
  return(z_mat)
}
z_list <- lapply(cor_list, fisher_z)

###按样本量加权整合 network
sample_sizes <- sapply(expr_list2, ncol)
sample_sizes
meta_network <- function(z_list, sample_sizes){
  weights <- sample_sizes - 3
  
  # 加权平均 z
  z_meta <- Reduce("+", Map(function(z, w) z * w, z_list, weights)) / sum(weights)
  
  # 转回相关系数
  r_meta <- tanh(z_meta)
  
  diag(r_meta) <- 1
  return(list(
    z_meta = z_meta,
    r_meta = r_meta,
    weights = weights
  ))
}
meta_res <- meta_network(z_list, sample_sizes)

meta_z <- meta_res$z_meta
meta_r <- meta_res$r_meta

##增加“方向一致性”指标
calc_sign_consistency <- function(cor_list){
  arr <- simplify2array(cor_list)   # gene × gene × dataset
  
  sign_arr <- sign(arr)
  
  pos_prop <- apply(sign_arr == 1, c(1,2), mean, na.rm = TRUE)
  neg_prop <- apply(sign_arr == -1, c(1,2), mean, na.rm = TRUE)
  
  consistency <- pmax(pos_prop, neg_prop, na.rm = TRUE)
  diag(consistency) <- 1
  
  return(consistency)
}

sign_consistency <- calc_sign_consistency(cor_list)

##定义 consensus edges
edge_threshold_r <- 0.3
edge_threshold_consistency <- 0.75

consensus_adj <- (abs(meta_r) >= edge_threshold_r) & (sign_consistency >= edge_threshold_consistency)
diag(consensus_adj) <- FALSE

write.table(
  data.frame(gene_id = rownames(meta_r), meta_r, check.names = FALSE),
  file = "meta_network_r.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  data.frame(gene_id = rownames(sign_consistency), sign_consistency, check.names = FALSE),
  file = "meta_network_sign_consistency.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  data.frame(gene_id = rownames(consensus_adj), consensus_adj, check.names = FALSE),
  file = "consensus_network_adj.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

get_edge_list <- function(meta_r, consistency, adj){
  idx <- which(adj, arr.ind = TRUE)
  idx <- idx[idx[,1] < idx[,2], , drop = FALSE]   # 只保留上三角
  
  edge_df <- data.frame(
    gene1 = rownames(meta_r)[idx[,1]],
    gene2 = colnames(meta_r)[idx[,2]],
    meta_r = meta_r[idx],
    consistency = consistency[idx],
    stringsAsFactors = FALSE
  )
  
  edge_df <- edge_df[order(-abs(edge_df$meta_r), -edge_df$consistency), ]
  return(edge_df)
}
edge_df <- get_edge_list(meta_r, sign_consistency, consensus_adj)

write.table(
  edge_df,
  file = "consensus_network_edge_list.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

dim(edge_df)
length(unique(c(edge_df$gene1,edge_df$gene2)))
bb <- unique(c(edge_df$gene1,edge_df$gene2))
core_genes <- read.csv("D:/phd cases/metabolic analysis/code/Deep_GLOC/run_GAT/core_gene/core_gene_all.csv")
intersect(bb,core_genes$gene_id)
setdiff(core_genes$gene_id,bb)
core_genes[core_genes$gene_id%in%bb,]


save.image("code_meta_graph_new.RData")
