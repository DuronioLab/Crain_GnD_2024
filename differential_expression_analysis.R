##FIGURE 3 - differential expression analysis Set8 and H4K20 mutants
#import all whole larvae RNA-seq data
samples <- read_tsv("sample_sheet_RNA-seq_3wl.txt")
rownames(samples) <- samples$useName

#set file names to match Salmon output
files <- file.path(samples$useName, "quant.sf")

#check if file labels are correct and exist in current directory
all(file.exists(files))
which(!file.exists(files))

#import gtf file
gtf <- import("dmel-all-r6.55.gtf") %>% 
  as.data.frame()

#filter for protein-coding genes
gtf_protein_genes <- gtf %>% 
  filter(type == "gene") %>% 
  filter(!grepl("RNA",gene_symbol)) %>% 
  filter(!grepl("mir-",gene_symbol)) %>% 
  filter(!grepl("CR[0-9]",gene_symbol)) %>% 
  unique()

#transcript name to gene name conversion
tx2gene <- gtf_protein_genes %>% 
  dplyr::select(transcript_id,gene_symbol) %>% 
  na.omit()

#functions to combine counts for replication-dependent histone gene transcripts from Salmon
compile_H3H4_genes <- function(gene_names) {
  combined_histone_genes <- gsub("^(His[0-9]+):.*$", "\\1", gene_names)
  return(combined_histone_genes)
}
compile_H2AB_genes <- function(gene_names) {
  combined_histone_genes <- gsub("^(His[0-9]+[AB]):.*$", "\\1", gene_names)
  return(combined_histone_genes)
}

#combine counts for RD histone gene transcripts from Salmon
tx2gene$gene_symbol <- sapply(tx2gene$gene_symbol,compile_H3H4_genes)
tx2gene$gene_symbol <- sapply(tx2gene$gene_symbol,compile_H2AB_genes)
tx2gene <- unique(tx2gene)

#import Salmon count files
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

#make design matrix for DESeq2
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ sample)

#Run DESeq2
ddsTxi.test <- DESeq(ddsTxi)

# #Data exploration (optional)
# #Normalize data and make PCA plot
# #filter low counts
# filt <- rowSums(assay(ddsTxi.test)) > 1
# 
# #Normalize count data using variance stabilizing transformation
# vsd <- vst(ddsTxi.test[filt,],blind = FALSE)
# 
# #plot PCA
# plotPCA(vsd,intgroup="sample")

#Differential expression analysis of pairwise comparisons
#histone mutants pairwise comparisons
results_K20A_HWT <- results(ddsTxi.test_all,contrast=c("sample","K20A","HWT"),alpha=0.05)
results_K20R_HWT <- results(ddsTxi.test_all,contrast=c("sample","K20R","HWT"),alpha=0.05)

#Set8 mutants pairwise comparisons
results_Set8null_OregonR <- results(ddsTxi.test_all,contrast=c("sample","Set8null","OregonR"),alpha=0.05) 
results_Set8rg_Set8wt <- results(ddsTxi.test_all,contrast=c("sample","Set8rg","Set8wt"),alpha=0.05)  

#log2FC shrinakge
#Set8 mutants
#Set8null
resLFC_Set8null_ashr2 <- lfcShrink(ddsTxi.test_all,
                                   contrast = c("sample","Set8null","OregonR"),
                                   res = results_Set8null_OregonR,
                                   type = "ashr")
resLFC_Set8null_ashr.df <- as.data.frame(resLFC_Set8null_ashr) %>% rownames_to_column("gene_symbol")

#Set8rg
resLFC_Set8rg_ashr <- lfcShrink(ddsTxi.test_all,
                                contrast = c("sample","Set8rg","Set8wt"),
                                res = results_Set8rg_Set8wt,
                                type = "ashr")
resLFC_Set8rg_ashr.df <- as.data.frame(resLFC_Set8rg_ashr) %>% rownames_to_column("gene_symbol")

#K20 mutants
resLFC_K20A_ashr <- lfcShrink(ddsTxi.test_all,
                              contrast = c("sample","K20A","HWT"),
                              res = results_K20A_HWT,
                              type = "ashr")
resLFC_K20A_ashr.df <- as.data.frame(resLFC_K20A_ashr) %>% rownames_to_column("gene_symbol")

resLFC_K20R_ashr <- lfcShrink(ddsTxi.test_all,
                              contrast = c("sample","K20R","HWT"),
                              res = results_K20R_HWT,
                              type = "ashr")
resLFC_K20R_ashr.df <- as.data.frame(resLFC_K20R_ashr) %>% rownames_to_column("gene_symbol")

# #MA plots (optional)
# plotMA(resLFC_Set8null_ashr,alpha=0.05,ylim=c(-10,10))
# 
# plotMA(resLFC_Set8null_Set8wt_ashr,alpha=0.05,ylim=c(-10,10))
# 
# plotMA(resLFC_Set8rg_ashr,alpha=0.05,ylim=c(-10,10))
# 
# plotMA(resLFC_K20A_ashr,alpha=0.05,ylim=c(-10,10))
# 
# plotMA(resLFC_K20R_ashr,alpha=0.05,ylim=c(-10,10))

#filter for significant genes
significant_genes <- function(results.df,l2fc,signif_thresh){
  sig_genes <- results.df %>% 
    filter(abs(log2FoldChange) >= l2fc) %>% 
    filter(padj <= signif_thresh) %>% 
    filter(gene_symbol %in% gtf__protein_genes$gene_symbol) %>% 
    na.omit()
  return(sig_genes)
}

#TABLE S1
#set8null
set8null_vs_OR_sig_shrunk <- significant_genes(resLFC_Set8null_ashr.df,1,0.01)
set8null_sig_down <- set8null_vs_OR_sig_shrunk %>% 
  filter(log2FoldChange < 0)
set8null_sig_up <- set8null_vs_OR_sig_shrunk %>% 
  filter(log2FoldChange > 0)
set8null_sig_unchanged <- resLFC_Set8null_ashr.df %>% 
  filter(!gene_symbol %in% set8null_vs_OR_sig_shrunk$gene_symbol)

#set8rg
set8rg_vs_set8wt_sig_shrunk <- significant_genes(resLFC_Set8rg_ashr.df,1,0.01)
set8rg_sig_down <- set8rg_vs_set8wt_sig_shrunk %>% 
  filter(log2FoldChange < 0)
set8rg_sig_up <- set8rg_vs_set8wt_sig_shrunk %>% 
  filter(log2FoldChange > 0)
set8rg_sig_unchanged <- resLFC_Set8rg_ashr.df %>% 
  filter(!gene_symbol %in% set8rg_vs_set8wt_sig_shrunk$gene_symbol)

#k20a
k20a_vs_hwt_sig_shrunk <- significant_genes(resLFC_K20A_ashr.df,1,0.01)
k20a_sig_down <- k20a_vs_hwt_sig_shrunk %>% 
  filter(log2FoldChange < 0)
k20a_sig_up <- k20a_vs_hwt_sig_shrunk %>% 
  filter(log2FoldChange > 0)
k20a_sig_unchanged <- resLFC_K20A_ashr.df %>% 
  filter(!gene_symbol %in% k20a_vs_hwt_sig_shrunk$gene_symbol)

#k20r
k20r_vs_hwt_sig_shrunk <- significant_genes(resLFC_K20R_ashr.df,1,0.01)
k20r_sig_down <- k20r_vs_hwt_sig_shrunk %>% 
  filter(log2FoldChange < 0)
k20r_sig_up <- k20r_vs_hwt_sig_shrunk %>% 
  filter(log2FoldChange > 0)
k20r_sig_unchanged <- resLFC_K20R_ashr.df %>% 
  filter(!gene_symbol %in% k20r_vs_hwt_sig_shrunk$gene_symbol)

##Export significant gene tables
make_bed_files <- function(gene_list_list) {
  for (gene_list_name in names(gene_list_list)) {
    gene_list <- as.data.frame(gene_list_list[[gene_list_name]])
    colnames(gene_list) <- sub("^[^.]*\\.", "", colnames(gene_list))
    gene_list.bed <- gtf_genes %>% 
      dplyr::filter(gene_symbol %in% gene_list$SYMBOL) %>%
      dplyr::select(seqnames, start, end, gene_symbol, width, strand) %>% 
      dplyr::filter(!seqnames == "mitochondrion_genome")
    
    gene_list.bed$seqnames <- paste("chr", gene_list.bed$seqnames, sep = "")
    write_tsv(gene_list.bed, file = paste0(gene_list_name, ".bed"), col_names = FALSE)
  }
}

list_for_bed <- list("Set8null.vs.Oregon-R.txt" = set8null_vs_OR_sig_shrunk,
                     "Set8RG.vs.Set8WT.txt" = set8rg_vs_set8wt_sig_shrunk,
                     "H4K20A.vs.HWT.txt" = k20a_vs_hwt_sig_shrunk,
                     "H4K20R.vs.HWT.txt" = k20r_vs_hwt_sig_shrunk
)

make_bed_files(list_for_bed)

#FIGURE 3A-D volcano plots
volcano <- function(results.df,l2fc,sig_thresh,filter_genes,filter_name,dot_cols) {
  for_volcano <- results.df %>% 
    mutate(Expression = case_when(log2FoldChange >= l2fc & padj <= sig_thresh ~ "Up-regulated",
                                  log2FoldChange <= -l2fc & padj <= sig_thresh ~ "Down-regulated",
                                  TRUE ~ "Unchanged")) %>% 
    mutate(status = if_else(gene_symbol %in% filter_genes, 
                            paste0(filter_name),
                            paste0("no",filter_name))) %>% 
    mutate(sig_status = case_when(gene_symbol %in% filter_genes & log2FoldChange >= l2fc & padj <= sig_thresh ~ paste0(filter_name,"_","up"),
                                  gene_symbol %in% filter_genes & log2FoldChange <= -l2fc & padj <= sig_thresh ~ paste0(filter_name,"_","down"),
                                  !gene_symbol %in% filter_genes & log2FoldChange >= l2fc & padj <= sig_thresh ~ paste0("no","_",filter_name,"_","up"),
                                  !gene_symbol %in% filter_genes & log2FoldChange <= -l2fc & padj <= sig_thresh ~ paste0("no","_",filter_name,"_","down"),
                                  TRUE ~ "unchanged")) %>% 
    mutate(outliers = case_when(log2FoldChange > 10 ~ "outlier",
                                log2FoldChange < -10 ~ "outlier",
                                padj < 1e-50 ~ "outlier",
                                TRUE ~ "not_outlier")) %>% 
    na.omit() %>% 
    arrange(desc(sig_status))
  
  for_volcano$log2FoldChange[for_volcano$log2FoldChange > 10] <- 10
  for_volcano$log2FoldChange[for_volcano$log2FoldChange < -10] <- -10
  for_volcano$padj[for_volcano$padj < 1e-50] <- 1e-50
  
  volcano_plot <- ggplot(for_volcano, aes(x=log2FoldChange,y=-log(padj,10),shape=outliers)) +
    geom_point(aes(color=sig_status), size=2.5) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"padj")) +
    scale_color_manual(values = dot_cols) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    scale_x_continuous(limits = c(-10,10)) +
    scale_y_continuous(limits = c(0,50)) +
    geom_vline(xintercept=1, 
               linetype="dashed", color = "black",linewidth=1) +
    geom_vline(xintercept=-1, 
               linetype="dashed", color = "black",linewidth=1) +
    geom_hline(yintercept=2, 
               linetype="dashed", color = "black",linewidth=1) +
    # coord_cartesian(xlim = c(x_lower, x_upper), 
    #                 ylim = c(y_lower, y_upper), 
    #                 clip = "off") +
    theme_classic2()
  
  return(volcano_plot)
}

dot_cols <- c("dodgerblue3", "firebrick3", "#99CCFF", "#EEAAAA", "grey")
dot_cols2 <- c("#99CCFF", "#EEAAAA", "grey")

#set8null
volcano(resLFC_Set8null_ashr.df,1,0.01,k20me1_genes$gene_symbol,"h4k20me1",dot_cols)

#set8rg
volcano(resLFC_Set8rg_ashr.df,1,0.01,k20me1_genes$gene_symbol,"h4k20me1",dot_cols)

#K20A
volcano(resLFC_K20A_ashr.df,1,0.01,k20me1_genes$gene_symbol,"h4k20me1",dot_cols2)

#K20R
volcano(resLFC_K20R_ashr.df,1,0.01,k20me1_genes$gene_symbol,"h4k20me1",dot_cols)

#TABLE 2
k20me1_genes <- read.tsv("k20me1_genes_0.5.bed")

nok20me1_genes <- gtf_protein_genes %>% 
  filter(!gene_name %in% k20me1_genes$gene_symbol)
colnames(nok20me1_genes)[11] <- "gene_symbol"

sig_k20me1_genes <- function(results.df,l2fc,signif_thresh){
  sig_genes <- results.df %>% 
    filter(abs(log2FoldChange) >= l2fc) %>% 
    filter(padj <= signif_thresh) %>% 
    filter(gene_symbol %in% gtf_protein_genes$gene_symbol) %>% 
    mutate(k20me1_status = case_when(
      gene_symbol %in% k20me1_genes$gene_symbol ~ "k20me1",
      TRUE ~ "no_k20me1")) %>% 
    na.omit()
  return(sig_genes)
}

#set8null
set8null_sig_k20me1 <- sig_k20me1_genes(resLFC_Set8null_ashr.df,1,0.01) %>% 
  filter(k20me1_status == "k20me1")
set8null_sig_nok20me1 <- sig_k20me1_genes(resLFC_Set8null_ashr.df,1,0.01) %>% 
  filter(k20me1_status == "no_k20me1")

set8null_sig_k20me1_up <- set8null_sig_k20me1 %>% 
  filter(log2FoldChange > 0)
set8null_sig_k20me1_down <- set8null_sig_k20me1 %>% 
  filter(log2FoldChange < 0)
set8null_sig_nok20me1_up <- set8null_sig_nok20me1 %>% 
  filter(log2FoldChange > 0)
set8null_sig_nok20me1_down <- set8null_sig_nok20me1 %>% 
  filter(log2FoldChange < 0)

#set8rg
set8rg_sig_k20me1 <- sig_k20me1_genes(resLFC_Set8rg_ashr.df,1,0.01) %>% 
  filter(k20me1_status == "k20me1")
set8rg_sig_nok20me1 <- sig_k20me1_genes(resLFC_Set8rg_ashr.df,1,0.01) %>% 
  filter(k20me1_status == "no_k20me1")

set8rg_sig_k20me1_up <- set8rg_sig_k20me1 %>% 
  filter(log2FoldChange > 0)
set8rg_sig_k20me1_down <- set8rg_sig_k20me1 %>% 
  filter(log2FoldChange < 0)
set8rg_sig_nok20me1_up <- set8rg_sig_nok20me1 %>% 
  filter(log2FoldChange > 0)
set8rg_sig_nok20me1_down <- set8rg_sig_nok20me1 %>% 
  filter(log2FoldChange < 0)

#K20A
K20A_sig_k20me1 <- sig_k20me1_genes(resLFC_K20A_ashr.df,1,0.01) %>% 
  filter(k20me1_status == "k20me1")
K20A_sig_nok20me1 <- sig_k20me1_genes(resLFC_K20A_ashr.df,1,0.01) %>% 
  filter(k20me1_status == "no_k20me1")

k20a_sig_k20me1_up <- K20A_sig_k20me1 %>% 
  filter(log2FoldChange > 0)
k20a_sig_k20me1_down <- K20A_sig_k20me1 %>% 
  filter(log2FoldChange < 0)
k20a_sig_nok20me1_up <- K20A_sig_nok20me1 %>% 
  filter(log2FoldChange > 0)
k20a_sig_nok20me1_down <- K20A_sig_nok20me1 %>% 
  filter(log2FoldChange < 0)

#K20R
K20R_sig_k20me1 <- sig_k20me1_genes(resLFC_K20R_ashr.df,1,0.01) %>% 
  filter(k20me1_status == "k20me1")
K20R_sig_nok20me1 <- sig_k20me1_genes(resLFC_K20R_ashr.df,1,0.01) %>% 
  filter(k20me1_status == "no_k20me1")

k20r_sig_k20me1_up <- K20R_sig_k20me1 %>% 
  filter(log2FoldChange > 0)
k20r_sig_k20me1_down <- K20R_sig_k20me1 %>% 
  filter(log2FoldChange < 0)
k20r_sig_nok20me1_up <- K20R_sig_nok20me1 %>% 
  filter(log2FoldChange > 0)
k20r_sig_nok20me1_down <- K20R_sig_nok20me1 %>% 
  filter(log2FoldChange < 0)

##Figure 3G - k-means clustered normalized centered count heatmaps - Set8 and H4K20 mutants
#make centered count heatmaps

#select normalized counts from appropriate genotypes
vsd_some.df <- vsd.df %>% 
  column_to_rownames("gene_symbol") %>% 
  select(contains(c("OregonR","Set8null","Set8wt","Set8rg","HWT","K20A","K20R")))

#calculate mean counts for controls across replicates
vsd_some_means.df <- vsd_some.df %>% 
  rowwise() %>% 
  mutate(OregonR_mean = mean(c_across(OregonR_rep1:OregonR_rep4))) %>% 
  mutate(HWT_mean = mean(c_across(HWT_rep1:HWT_rep4))) %>% 
  mutate(Set8wt_mean = mean(c_across(Set8wt_rep1:Set8wt_rep4)))

#subtract control mean from replicates to obtain centered counts relative to control
vsd_some_means_center.df <- vsd_some_means.df %>%
  mutate_at(vars(contains("Set8null") | contains("OregonR_rep")), 
            ~ . - OregonR_mean) %>%
  mutate_at(vars(contains("HWT_rep") | contains("K20A") | contains("K20R")), 
            ~ . - HWT_mean) %>%
  mutate_at(vars(contains("Set8wt_rep") | contains("Set8rg")), 
            ~ . - Set8wt_mean)

#filter out control mean columns and clean up df
vsd_some_ctrl_center <- vsd_some_means_center.df %>% 
  select(!contains("mean"))
rownames(vsd_some_ctrl_center) <- rownames(vsd_some.df)
vsd_some_ctrl_center.mat <- as.matrix(vsd_some_ctrl_center)

#filter df for k20me1 genes and nok20me1 genes
vsd_some_center_k20me1.mat <- vsd_some_ctrl_center %>% 
  rownames_to_column("SYMBOL") %>% 
  filter(SYMBOL %in% k20me1_genes$gene_symbol) %>%
  column_to_rownames("SYMBOL") %>% 
  as.matrix()

vsd_some_center_nok20me1.mat <- vsd_some_ctrl_center %>% 
  rownames_to_column("SYMBOL") %>% 
  filter(SYMBOL %in% nok20me1_genes$gene_symbol) %>% 
  column_to_rownames("SYMBOL") %>% 
  as.matrix()

##Determine appropriate number of clusters using WCSS method
# Function to calculate WCSS for a given number of clusters (k)
calculate_wcss <- function(data, k) {
  k_values <- 
    kmeans(data, centers = k)$tot.withinss
}

# Determine a range of k values (e.g., from 2 to 10)
k_values <- 2:10

# Calculate WCSS for each k
wcss_values <- map_dbl(k_values, ~calculate_wcss(vsd_some_center_k20me1.mat, .))

# Create a data frame for the scree plot
scree_data <- data.frame(k = k_values, WCSS = wcss_values)

# Plot the scree plot
ggplot(scree_data, aes(x = k, y = WCSS)) +
  geom_line() +
  geom_point() +
  labs(title = "Scree Plot for K-means Clustering",
       x = "Number of Clusters (k)",
       y = "Within-Cluster Sum of Squares (WCSS)")

#FIGURE 3E - make heatmaps
library(ComplexHeatmap)

col_fun <- colorRamp2(c(-1,0,1), c("dodgerblue3","white","firebrick"))
norm_counts_k20me1_htmap <- Heatmap(vsd_some_center_k20me1.mat,
                                    col = col_fun,
                                    cluster_rows = TRUE,
                                    row_km = 6,
                                    row_km_repeats = 20,
                                    cluster_columns = FALSE,
                                    show_row_names = FALSE,
                                    name = "K20me1 centered counts"
)
set.seed(122)
norm_counts_k20me1_htmap

#extract genes from clusters in heatmaps
heatmap_extract_cluster <- function(heatmap_obj, matrix_obj, which = "row",seed_num){
  set.seed(seed_num)
  if(which == "column") c <- ComplexHeatmap::column_order(heatmap_obj)
  if(which == "row") c <- ComplexHeatmap::row_order(heatmap_obj)
  n <- names(c)
  l <- seq(1, length(n), 1)
  clu <- purrr::map2(c, l, function(x, y){
    data.frame(x, y)
  })
  clu <- dplyr::bind_rows(clu)
  clu <- dplyr::arrange(clu, x)
  clu <- dplyr::mutate(clu, y = as.character(y))
  if(which == "column") r <- tibble::tibble(SYMBOL = colnames(matrix_obj), cluster = clu$y)
  if(which == "row") r <- tibble::tibble(SYMBOL = rownames(matrix_obj), cluster = clu$y)
  return(r)
}

#k20me1 genes
k20me1_gene_clusters <- heatmap_extract_cluster(norm_counts_k20me1_htmap,
                                                vsd_some_center_k20me1.mat,
                                                "row",
                                                122)
k20me1_normCounts_cluster1 <- k20me1_gene_clusters %>% 
  filter(cluster == 1)
k20me1_normCounts_cluster2 <- k20me1_gene_clusters %>% 
  filter(cluster == 2)
k20me1_normCounts_cluster3 <- k20me1_gene_clusters %>% 
  filter(cluster == 5)
k20me1_normCounts_cluster4 <- k20me1_gene_clusters %>% 
  filter(cluster == 6)
k20me1_normCounts_cluster5 <- k20me1_gene_clusters %>% 
  filter(cluster == 4)
k20me1_normCounts_cluster6 <- k20me1_gene_clusters %>% 
  filter(cluster == 3)

#GO terms
k20me1_normCounts_cluster1_go <- go_terms(k20me1_normCounts_cluster1$SYMBOL)
k20me1_normCounts_cluster2_go <- go_terms(k20me1_normCounts_cluster2$SYMBOL)
k20me1_normCounts_cluster3_go <- go_terms(k20me1_normCounts_cluster3$SYMBOL)
k20me1_normCounts_cluster4_go <- go_terms(k20me1_normCounts_cluster4$SYMBOL)
k20me1_normCounts_cluster5_go <- go_terms(k20me1_normCounts_cluster5$SYMBOL)
k20me1_normCounts_cluster6_go <- go_terms(k20me1_normCounts_cluster6$SYMBOL)

#Supplemental table S3 - Gene ontology terms for clusters in FIGURE 3E HIGH H4K20me1
write_tsv(k20me1_normCounts_cluster1_go[,1:11],file = "k20me1_normCounts_cluster1_go.txt", col_names = TRUE)
write_tsv(k20me1_normCounts_cluster2_go[,1:11],file = "k20me1_normCounts_cluster2_go.txt", col_names = TRUE)
write_tsv(k20me1_normCounts_cluster3_go[,1:11],file = "k20me1_normCounts_cluster3_go.txt", col_names = TRUE)
write_tsv(k20me1_normCounts_cluster4_go[,1:11],file = "k20me1_normCounts_cluster4_go.txt", col_names = TRUE)
write_tsv(k20me1_normCounts_cluster5_go[,1:11],file = "k20me1_normCounts_cluster5_go.txt", col_names = TRUE)
write_tsv(k20me1_normCounts_cluster6_go[,1:11],file = "k20me1_normCounts_cluster6_go.txt", col_names = TRUE)

norm_counts_htmap_nok20me1 <- Heatmap(vsd_some_center_nok20me1.mat,
                                      col = col_fun,
                                      cluster_rows = TRUE,
                                      row_km = 6,
                                      row_km_repeats = 20,
                                      cluster_columns = FALSE,
                                      show_row_names = FALSE,
                                      name = "nok20me1 centered counts"
)
set.seed(123)
norm_counts_htmap_nok20me1

#nok20me1 genes
nok20me1_gene_clusters <- heatmap_extract_cluster(norm_counts_htmap_nok20me1,
                                                  vsd_some_center_nok20me1.mat,
                                                  "row",
                                                  123)

nok20me1_normCounts_cluster1 <- nok20me1_gene_clusters %>% 
  filter(cluster == 2)
nok20me1_normCounts_cluster2 <- nok20me1_gene_clusters %>% 
  filter(cluster == 4)
nok20me1_normCounts_cluster3 <- nok20me1_gene_clusters %>% 
  filter(cluster == 5)
nok20me1_normCounts_cluster4 <- nok20me1_gene_clusters %>% 
  filter(cluster == 6)
nok20me1_normCounts_cluster5 <- nok20me1_gene_clusters %>% 
  filter(cluster == 3)
nok20me1_normCounts_cluster6 <- nok20me1_gene_clusters %>% 
  filter(cluster == 1)

#GO terms
nok20me1_normCounts_cluster1_go <- go_terms(nok20me1_normCounts_cluster1$SYMBOL)
nok20me1_normCounts_cluster2_go <- go_terms(nok20me1_normCounts_cluster2$SYMBOL)
nok20me1_normCounts_cluster3_go <- go_terms(nok20me1_normCounts_cluster3$SYMBOL)
nok20me1_normCounts_cluster4_go <- go_terms(nok20me1_normCounts_cluster4$SYMBOL)
nok20me1_normCounts_cluster5_go <- go_terms(nok20me1_normCounts_cluster5$SYMBOL)
nok20me1_normCounts_cluster6_go <- go_terms(nok20me1_normCounts_cluster6$SYMBOL)

#Supplemental table S4 - Gene ontology terms for clusters in FIGURE 3E LOW H4K20me1
write_tsv(nok20me1_normCounts_cluster1_go[,1:11],file = "nok20me1_normCounts_cluster1_go.txt", col_names = TRUE)
write_tsv(nok20me1_normCounts_cluster2_go[,1:11],file = "nok20me1_normCounts_cluster2_go.txt", col_names = TRUE)
write_tsv(nok20me1_normCounts_cluster3_go[,1:11],file = "nok20me1_normCounts_cluster3_go.txt", col_names = TRUE)
write_tsv(nok20me1_normCounts_cluster4_go[,1:11],file = "nok20me1_normCounts_cluster4_go.txt", col_names = TRUE)
write_tsv(nok20me1_normCounts_cluster5_go[,1:11],file = "nok20me1_normCounts_cluster5_go.txt", col_names = TRUE)
write_tsv(nok20me1_normCounts_cluster6_go[,1:11],file = "nok20me1_normCounts_cluster6_go.txt", col_names = TRUE)

norm_counts_htmap_all <- Heatmap(vsd_some_ctrl_center.mat,
                                 col = col_fun,
                                 cluster_rows = TRUE,
                                 row_km = 6,
                                 row_km_repeats = 20,
                                 cluster_columns = FALSE,
                                 show_row_names = FALSE,
                                 use_raster = FALSE,
                                 name = "Counts"
)
set.seed(121)
norm_counts_htmap_all

#all genes
all_gene_clusters <- heatmap_extract_cluster(norm_counts_htmap_all,
                                             vsd_some_ctrl_center.mat,
                                             "row",
                                             121)

all_normCounts_cluster1 <- all_gene_clusters %>% 
  filter(cluster == 2)
all_normCounts_cluster2 <- all_gene_clusters %>% 
  filter(cluster == 5)
all_normCounts_cluster3 <- all_gene_clusters %>% 
  filter(cluster == 6)
all_normCounts_cluster4 <- all_gene_clusters %>% 
  filter(cluster == 4)
all_normCounts_cluster5 <- all_gene_clusters %>% 
  filter(cluster == 3)
all_normCounts_cluster6 <- all_gene_clusters %>% 
  filter(cluster == 1)

#GO terms
all_normCounts_cluster1_go <- go_terms(all_normCounts_cluster1$SYMBOL)
all_normCounts_cluster2_go <- go_terms(all_normCounts_cluster2$SYMBOL)
all_normCounts_cluster3_go <- go_terms(all_normCounts_cluster3$SYMBOL)
all_normCounts_cluster4_go <- go_terms(all_normCounts_cluster4$SYMBOL)
all_normCounts_cluster5_go <- go_terms(all_normCounts_cluster5$SYMBOL)
all_normCounts_cluster6_go <- go_terms(all_normCounts_cluster6$SYMBOL)

#Supplemental table S5 - Gene ontology terms for clusters in FIGURE 3E ALL genes
write_tsv(all_normCounts_cluster1_go[,1:11],file = "all_normCounts_cluster1_go.txt", col_names = TRUE)
write_tsv(all_normCounts_cluster2_go[,1:11],file = "all_normCounts_cluster2_go.txt", col_names = TRUE)
write_tsv(all_normCounts_cluster3_go[,1:11],file = "all_normCounts_cluster3_go.txt", col_names = TRUE)
write_tsv(all_normCounts_cluster4_go[,1:11],file = "all_normCounts_cluster4_go.txt", col_names = TRUE)
write_tsv(all_normCounts_cluster5_go[,1:11],file = "all_normCounts_cluster5_go.txt", col_names = TRUE)
write_tsv(all_normCounts_cluster6_go[,1:11],file = "all_normCounts_cluster6_go.txt", col_names = TRUE)

#FIGURE 5 - Differential expression analysis of L(3)mbt mutants
#FIGURE 5A-C - volcano plots

#L(3)mbt mutants
results_L3WT_25_OregonR <- results(ddsTxi.test_all,contrast=c("sample","L3mbtWT_Df25","OregonR"),alpha=0.05)
results_L3GM76_25_OregonR <- results(ddsTxi.test_all,contrast=c("sample","L3mbtGM76_Df25","OregonR"),alpha=0.05)
results_L3PBac_25_OregonR <- results(ddsTxi.test_all,contrast=c("sample","L3mbtPBac_Df25","OregonR"),alpha=0.05)
results_L3GM76_29_L3GM76_25 <- results(ddsTxi.test_all,contrast=c("sample","L3mbtGM76_Df29","L3mbtGM76_Df25"),alpha=0.05)
results_L3PBac_29_L3PBac_25 <- results(ddsTxi.test_all,contrast=c("sample","L3mbtPBac_Df29","L3mbtPBac_Df25"),alpha=0.05)

#L2FC shrinkage
resLFC_L3WT_OR_ashr <- lfcShrink(ddsTxi.test_all,
                                 contrast = c("sample","L3mbtWT_Df25","OregonR"),
                                 res = results_L3WT_25_OregonR,
                                 type = "ashr")
resLFC_L3WT_OR_ashr.df <- as.data.frame(resLFC_L3WT_OR_ashr) %>% rownames_to_column("gene_symbol")

resLFC_L3GM76_OR_ashr <- lfcShrink(ddsTxi.test_all,
                                   contrast = c("sample","L3mbtGM76_Df25","OregonR"),
                                   res = results_L3GM76_25_OregonR,
                                   type = "ashr")
resLFC_L3GM76_OR_ashr.df <- as.data.frame(resLFC_L3GM76_OR_ashr) %>% rownames_to_column("gene_symbol")

resLFC_L3PBac_25_OR_ashr <- lfcShrink(ddsTxi.test_all,
                                      contrast = c("sample","L3mbtPBac_Df25","OregonR"),
                                      res = results_L3PBac_25_OregonR,
                                      type = "ashr")
resLFC_L3PBac_25_OR_ashr.df <- as.data.frame(resLFC_L3PBac_25_OR_ashr) %>% rownames_to_column("gene_symbol")

#29C vs. 25C
resLFC_L3WT_29_L3WT_25_ashr <- lfcShrink(ddsTxi.test_all,
                                         contrast = c("sample","L3mbtWT_Df29","L3mbtWT_Df25"),
                                         res = results_L3WT_29_L3WT_25,
                                         type = "ashr")
resLFC_L3WT_29_L3WT_25_ashr.df <- as.data.frame(resLFC_L3WT_29_L3WT_25_ashr) %>% rownames_to_column("gene_symbol")

resLFC_L3GM76_29_L3GM76_25_ashr <- lfcShrink(ddsTxi.test_all,
                                             contrast = c("sample","L3mbtGM76_Df29","L3mbtGM76_Df25"),
                                             res = results_L3GM76_29_L3GM76_25,
                                             type = "ashr")
resLFC_L3GM76_29_L3GM76_25_ashr.df <- as.data.frame(resLFC_L3GM76_29_L3GM76_25_ashr) %>% rownames_to_column("gene_symbol")

resLFC_L3PBac_29_L3PBac_25_ashr <- lfcShrink(ddsTxi.test_all,
                                             contrast = c("sample","L3mbtPBac_Df29","L3mbtPBac_Df25"),
                                             res = results_L3PBac_29_L3PBac_25,
                                             type = "ashr")
resLFC_L3PBac_29_L3PBac_25_ashr.df <- as.data.frame(resLFC_L3PBac_29_L3PBac_25_ashr) %>% rownames_to_column("gene_symbol")

# #MA plots
# plotMA(resLFC_L3WT_OR_ashr,alpha=0.05,ylim=c(-10,10))
# plotMA(resLFC_L3GM76_OR_ashr,alpha=0.05,ylim=c(-10,10))
# plotMA(resLFC_L3PBac_25_OR_ashr,alpha=0.05,ylim=c(-10,10))
# 
# plotMA(resLFC_L3WT_29_L3WT_25_ashr,alpha=0.05,ylim=c(-10,10))
# plotMA(resLFC_L3GM76_29_L3GM76_25_ashr,alpha=0.05,ylim=c(-10,10))
# plotMA(resLFC_L3PBac_29_L3PBac_25_ashr,alpha=0.05,ylim=c(-10,10))

#FIGURE 5A-C - L(3)mbt volcano plots
volcano(resLFC_L3WT_OR_ashr.df,1,0.01,k20me1_genes$gene_symbol,"k20me1", dot_cols)

volcano(resLFC_L3GM76_OR_ashr.df,1,0.01, k20me1_genes$gene_symbol,"k20me1", dot_cols)

volcano(resLFC_L3PBac_25_OR_ashr.df,1,0.01,k20me1_genes$gene_symbol,"k20me1",dot_cols)

##significant L(3)mbt genes
#TABLE S6
#L(3)mbtWT_25 vs OR
L3WT_25_vs_OR_sig_shrunk <- significant_genes(resLFC_L3WT_OR_ashr.df,1,0.01)
L3WT_25_vs_OR_sig_down <- L3WT_25_vs_OR_sig_shrunk %>% 
  filter(log2FoldChange < 0)
L3WT_25_vs_OR_sig_up <- L3WT_25_vs_OR_sig_shrunk %>% 
  filter(log2FoldChange > 0)

#L(3)mbtGM76_25 vs OR
L3GM76_25_vs_OR_sig_shrunk <- significant_genes(resLFC_L3GM76_OR_ashr.df,1,0.01)
L3GM76_25_vs_OR_sig_down <- L3GM76_25_vs_OR_sig_shrunk %>% 
  filter(log2FoldChange < 0)
L3GM76_25_vs_OR_sig_up <- L3GM76_25_vs_OR_sig_shrunk %>% 
  filter(log2FoldChange > 0)

#L(3)mbtWT_29 vs WT_25
L3WT_29_vs_L3WT_25_sig_shrunk <- significant_genes(resLFC_L3WT_29_L3WT_25_ashr.df,1,0.01)
L3WT_29_vs_L3WT_25_sig_down <- L3WT_29_vs_L3WT_25_sig_shrunk %>% 
  filter(log2FoldChange < 0)
L3WT_29_vs_L3WT_25_sig_up <- L3WT_29_vs_L3WT_25_sig_shrunk %>% 
  filter(log2FoldChange > 0)

#L(3)mbtGM76_29 vs GM76_25
L3GM76_29_vs_L3GM76_25_sig_shrunk <- significant_genes(resLFC_L3GM76_29_L3GM76_25_ashr.df,1,0.01)
L3GM76_29_vs_L3GM76_25_sig_down <- L3GM76_29_vs_L3GM76_25_sig_shrunk %>% 
  filter(log2FoldChange < 0)
L3GM76_29_vs_L3GM76_25_sig_up <- L3GM76_29_vs_L3GM76_25_sig_shrunk %>% 
  filter(log2FoldChange > 0)

#L(3)mbtPBac_25 vs OR
L3PBac_25_vs_OR_sig_shrunk <- significant_genes(resLFC_L3PBac_25_OR_ashr.df,1,0.01)
L3PBac_25_vs_OR_sig_down <- L3PBac_25_vs_OR_sig_shrunk %>% 
  filter(log2FoldChange < 0)
L3PBac_25_vs_OR_sig_up <- L3PBac_25_vs_OR_sig_shrunk %>% 
  filter(log2FoldChange > 0)

#L(3)mbtPBac_29 vs L(3)mbtPBac_25
L3PBac_29_vs_L3PBac_25_sig_shrunk <- significant_genes(resLFC_L3PBac_29_L3PBac_25_ashr.df,1,0.01)
L3PBac_29_vs_L3PBac_25_sig_down <- L3PBac_29_vs_L3PBac_25_sig_shrunk %>% 
  filter(log2FoldChange < 0)
L3PBac_29_vs_L3PBac_25_sig_up <- L3PBac_29_vs_L3PBac_25_sig_shrunk %>% 
  filter(log2FoldChange > 0)

#FIGURE 5D - k-means clustered normalized centered count heatmaps - L(3)mbt mutants
vsd_more.df <- vsd_all.df %>% 
  column_to_rownames("gene_symbol") %>% 
  select(contains(c("OregonR","Set8null","HWT","K20R","Df25")))

vsd_more_means.df <- vsd_more.df %>% 
  rowwise() %>% 
  mutate(OregonR_mean = mean(c_across(OregonR_rep1:OregonR_rep4))) %>% 
  mutate(HWT_mean = mean(c_across(HWT_rep1:HWT_rep4)))

vsd_more_means_center.df <- vsd_more_means.df %>%
  mutate_at(vars(contains("Set8null") | contains("OregonR_rep")), 
            ~ . - OregonR_mean) %>%
  mutate_at(vars(contains("HWT_rep") | contains("K20A") | contains("K20R")), 
            ~ . - HWT_mean) %>%
  mutate_at(vars(contains("Df25")), 
            ~ . - OregonR_mean)

vsd_more_ctrl_center <- vsd_more_means_center.df %>% 
  select(!contains("mean"))
vsd_more_ctrl_center <- vsd_more_ctrl_center[,c("OregonR_rep1","OregonR_rep2","OregonR_rep3","OregonR_rep4",
                                                "Set8null_rep1","Set8null_rep2","Set8null_rep3","Set8null_rep4",
                                                "HWT_rep1","HWT_rep2","HWT_rep3","HWT_rep4",
                                                "K20R_rep1","K20R_rep2","K20R_rep3","K20R_rep4",
                                                "L3mbtWT_Df25_rep1","L3mbtWT_Df25_rep2","L3mbtWT_Df25_rep3","L3mbtGM76_Df25_rep1",
                                                "L3mbtGM76_Df25_rep2","L3mbtGM76_Df25_rep3","L3mbtGM76_Df25_rep4",
                                                "L3mbtPBac_Df25_rep1","L3mbtPBac_Df25_rep2","L3mbtPBac_Df25_rep3")]
rownames(vsd_more_ctrl_center) <- rownames(vsd_more.df)
vsd_more_ctrl_center.mat <- as.matrix(vsd_more_ctrl_center)

vsd_more_ctrl_center_k20me1.mat <- vsd_more_ctrl_center %>% 
  rownames_to_column("gene_symbol") %>% 
  filter(gene_symbol %in% k20me1_genes$gene_symbol) %>% 
  column_to_rownames("gene_symbol") %>% 
  as.matrix()

vsd_more_ctrl_center_nok20me1.mat <- vsd_more_ctrl_center %>% 
  rownames_to_column("gene_symbol") %>% 
  filter(!gene_symbol %in% k20me1_genes$gene_symbol) %>% 
  column_to_rownames("gene_symbol") %>% 
  as.matrix()

norm_counts_k20me1_htmap_wL3mbt <- Heatmap(vsd_more_ctrl_center_k20me1.mat,
                                           col = col_fun,
                                           cluster_rows = TRUE,
                                           row_km = 6,
                                           row_km_repeats = 10,
                                           cluster_columns = FALSE,
                                           show_row_names = FALSE,
                                           name = "centered counts"
)
set.seed(122)
norm_counts_k20me1_htmap_wL3mbt

norm_counts_nok20me1_htmap_wL3mbt <- Heatmap(vsd_more_ctrl_center_nok20me1.mat,
                                             col = col_fun,
                                             cluster_rows = TRUE,
                                             row_km = 6,
                                             row_km_repeats = 10,
                                             cluster_columns = FALSE,
                                             show_row_names = FALSE,
                                             name = "centered counts"
)
set.seed(122)
norm_counts_nok20me1_htmap_wL3mbt

norm_counts_htmap_wL3mbt <- Heatmap(vsd_more_ctrl_center.mat,
                                    col = col_fun,
                                    cluster_rows = TRUE,
                                    row_km = 6,
                                    row_km_repeats = 10,
                                    cluster_columns = FALSE,
                                    show_row_names = FALSE,
                                    name = "K20me1 centered counts"
)
set.seed(122)
norm_counts_htmap_wL3mbt

#FIGURE 5E - MBT tumor gene analysis
##combine pariwise comparison results for l2fc analysis
set8null_l2fc_ashr_for_merge <- resLFC_Set8null_ashr.df %>% 
  select(gene_symbol,log2FoldChange,padj)
colnames(set8null_l2fc_ashr_for_merge) <- c("SYMBOL","set8null_l2fc","set8null_padj")

k20r_l2fc_ashr_for_merge <- resLFC_K20R_ashr.df %>% 
  select(gene_symbol,log2FoldChange,padj)
colnames(k20r_l2fc_ashr_for_merge) <- c("SYMBOL","k20r_l2fc","k20r_padj")

L3WT_OR_l2fc_ashr_for_merge <- resLFC_L3WT_OR_ashr.df %>% 
  select(gene_symbol,log2FoldChange,padj)
colnames(L3WT_OR_l2fc_ashr_for_merge) <- c("SYMBOL","L3WT_OR_l2fc","L3WT_OR_padj")

L3GM76_OR_l2fc_ashr_for_merge <- resLFC_L3GM76_OR_ashr.df %>% 
  select(gene_symbol,log2FoldChange,padj)
colnames(L3GM76_OR_l2fc_ashr_for_merge) <- c("SYMBOL","L3GM76_OR_l2fc","L3GM76_OR_padj")

L3PBac_OR_l2fc_ashr_for_merge <- resLFC_L3PBac_25_OR_ashr.df %>% 
  select(gene_symbol,log2FoldChange,padj)
colnames(L3PBac_OR_l2fc_ashr_for_merge) <- c("SYMBOL","L3PBac_OR_l2fc","L3PBac_OR_padj")

all_l2fc_ashr <- set8null_l2fc_ashr_for_merge %>% 
  left_join(k20r_l2fc_ashr_for_merge) %>% 
  left_join(L3WT_OR_l2fc_ashr_for_merge) %>% 
  left_join(L3GM76_OR_l2fc_ashr_for_merge) %>% 
  left_join(L3PBac_OR_l2fc_ashr_for_merge)

#import tumor genes list from Janic et al.
mbt_tumor_genes <- read_tsv("mbt_tumor_genes.txt") %>% 
  dplyr::select(gene_symbol) %>% 
  na.omit()

all_l2fc_ashr.mat <- all_l2fc_ashr %>%
  filter(SYMBOL %in% mbt_tumor_genes$gene_symbol) %>% 
  column_to_rownames("SYMBOL") %>%
  dplyr::select(set8null_l2fc,
                k20a_l2fc,
                k20r_l2fc,
                L3WT_OR_l2fc,
                L3GM76_OR_l2fc,
                L3PBac_OR_l2fc
  ) %>%
  as.matrix() %>% 
  na.omit()

col_fun <- colorRamp2(c(-2,0,6), c("dodgerblue3","white","firebrick"))
shared_DEGs_htmap_mbt <- Heatmap(all_l2fc_ashr.mat,
                                 col = col_fun,
                                 cluster_rows = TRUE,
                                 row_km = 6,
                                 row_km_repeats = 10,
                                 cluster_columns = FALSE,
                                 show_row_names = FALSE,
                                 name = "Log2 Fold Change"
)
set.seed(123)
shared_DEGs_htmap_mbt