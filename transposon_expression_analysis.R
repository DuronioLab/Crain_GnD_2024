#Supplemental FIGURE S1B-E - Transposon and piRNA cluster expression analysis
samples_transposon <- read.tsv("sample_sheet_RNA-seq_transposon_analysis.txt")

#generate transposon tx2gene (from https://ftp.flybase.net/releases/current/precomputed_files/transposons/transposon_sequence.gff)
transposon_gff <- import("transposon_sequence_set.gff") %>% 
  as.data.frame()

tx2gene_te <- transposon_gff %>% 
  mutate(te_name = gsub(".*\\\\", "", name)) %>% 
  dplyr::select(ID,te_name) %>% 
  na.omit()

#generate piRNA cluster tx2gene (modified from https://www.smallrnagroup.uni-mainz.de/piRNAclusterDB/data/FASTA/Drosophila_melanogaster.piRNAclusters.gtf)
piRNA_clusters_gff <- read.tsv("Dm_piRNA_clusters.gff.txt")

tx2gene_piRNAs <- piRNA_clusters_gff %>% 
  dplyr::select(id,name) %>% 
  na.omit()

#make combined tx2gene
combined_tx2gene <- bind_rows(tx2gene,tx2gene_te,tx2gene_piRNAs)

txi_all <- tximport(files_wK9R, type="salmon", tx2gene=combined_tx2gene)
txi_all.df <- as.data.frame(txi_all)

ddsTxi_all_genes <- DESeqDataSetFromTximport(txi_all,
                                             colData = samples_wK9R,
                                             design = ~ sample)

summary(ddsTxi_all_genes)
ddsTxi_all_genes.df <- as.data.frame(assay(ddsTxi_all_genes)) %>% 
  rownames_to_column("gene_symbol")

ddsTxi_all_genes.test <- DESeq(ddsTxi_all_genes)
ddsTxi_all_genes.test.df <- as.data.frame(assay(ddsTxi_all_genes.test)) %>% 
  rownames_to_column("gene_symbol")

#differential expression analysis
#histone mutants
results_all_genes_K20A_HWT <- results(ddsTxi_all_genes.test,contrast=c("sample","K20A","HWT"),alpha=0.05)
results_all_genes_K20R_HWT <- results(ddsTxi_all_genes.test,contrast=c("sample","K20R","HWT"),alpha=0.05)  
results_all_genes_K9R_HWT <- results(ddsTxi_all_genes.test,contrast=c("sample","RNA-K9R","RNA-WT"),alpha=0.05)

#Set8 mutants
results_all_genes_Set8null_OregonR <- results(ddsTxi_all_genes.test,contrast=c("sample","Set8null","OregonR"),alpha=0.05) 

#log2FC shrinakge
#Set8 mutants
resLFC_all_genes_Set8null_ashr <- lfcShrink(ddsTxi_all_genes.test,
                                            contrast = c("sample","Set8null","OregonR"),
                                            res = results_all_genes_Set8null_OregonR,
                                            type = "ashr")

#K20 mutants
resLFC_all_genes_K20A_ashr <- lfcShrink(ddsTxi_all_genes.test,
                                        contrast = c("sample","K20A","HWT"),
                                        res = results_all_genes_K20A_HWT,
                                        type = "ashr")
resLFC_all_genes_K20A_ashr.df <- as.data.frame(resLFC_all_genes_K20A_ashr) %>% rownames_to_column("gene_symbol")

resLFC_all_genes_K20R_ashr <- lfcShrink(ddsTxi_all_genes.test,
                                        contrast = c("sample","K20R","HWT"),
                                        res = results_all_genes_K20R_HWT,
                                        type = "ashr")
resLFC_all_genes_K20R_ashr.df <- as.data.frame(resLFC_all_genes_K20R_ashr) %>% rownames_to_column("gene_symbol")

#K9R
resLFC_all_genes_K9R_ashr <- lfcShrink(ddsTxi_all_genes.test,
                                       contrast = c("sample","RNA-K9R","RNA-WT"),
                                       res = results_all_genes_K9R_HWT,
                                       type = "ashr")
resLFC_all_genes_K9R_ashr.df <- as.data.frame(resLFC_all_genes_K9R_ashr) %>% rownames_to_column("gene_symbol")

# #MA plots
# plotMA(resLFC_all_genes_Set8null_ashr,alpha=0.05,ylim=c(-10,10))
# 
# plotMA(resLFC_all_genes_K20A_ashr,alpha=0.05,ylim=c(-10,10))
# 
# plotMA(resLFC_all_genes_K20R_ashr,alpha=0.05,ylim=c(-10,10))
# 
# plotMA(resLFC_all_genes_K9R_ashr,alpha=0.05,ylim=c(-10,10))

#volcano plots - transposons
transposon_cols <- c("dodgerblue3", "firebrick3", "grey")
volcano_transposon <- function(results.df,l2fc,sig_thresh,dot_cols) {
  for_volcano <- results.df %>% 
    left_join(tx2gene_te_expanded) %>% 
    filter(gene_symbol %in% tx2gene_te$gene_symbol |
             gene_symbol %in% tx2gene_piRNAs$gene_symbol) %>%
    mutate(Expression = case_when(log2FoldChange >= l2fc & padj <= sig_thresh ~ "Up-regulated",
                                  log2FoldChange <= -l2fc & padj <= sig_thresh ~ "Down-regulated",
                                  TRUE ~ "Unchanged")) %>% 
    mutate(outliers = case_when(log2FoldChange > 10 ~ "outlier",
                                log2FoldChange < -10 ~ "outlier",
                                padj < 1e-50 ~ "outlier",
                                TRUE ~ "not_outlier")) %>% 
    mutate(sig_subtype = case_when(abs(log2FoldChange) < l2fc ~ "Unchanged",
                                   padj > sig_thresh | is.na(padj)  ~ "Unchanged",
                                   gene_symbol %in% tx2gene_te$gene_symbol ~ "transposon",
                                   gene_symbol %in% tx2gene_piRNAs$gene_symbol ~ "piRNA")) %>% 
    arrange(desc(sig_subtype))
  
  for_volcano$log2FoldChange[for_volcano$log2FoldChange > 5] <- 5
  for_volcano$log2FoldChange[for_volcano$log2FoldChange < -5] <- -5
  for_volcano$padj[for_volcano$padj < 1e-30] <- 1e-30
  for_volcano$padj[is.na(for_volcano$padj)] <- 1
  
  volcano_plot <- ggplot(for_volcano, aes(x=log2FoldChange,y=-log(padj,10),shape=outliers)) +
    geom_point(aes(color=sig_subtype), size=2.5) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"padj")) +
    scale_color_manual(values = dot_cols) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    scale_x_continuous(limits = c(-5,5)) +
    scale_y_continuous(limits = c(0,30)) +
    geom_vline(xintercept=1, 
               linetype="dashed", color = "black",size=1) +
    geom_vline(xintercept=-1, 
               linetype="dashed", color = "black",size=1) +
    geom_hline(yintercept=2, 
               linetype="dashed", color = "black",size=1) +
    theme_classic2()
  
  return(volcano_plot)
}

volcano_transposon(resLFC_all_genes_Set8null_ashr.df,1,0.01,transposon_cols)
volcano_transposon(resLFC_all_genes_K20A_ashr.df,1,0.01,transposon_cols)
volcano_transposon(resLFC_all_genes_K20R_ashr.df,1,0.01,transposon_cols)
volcano_transposon(resLFC_all_genes_K9R_ashr.df,1,0.01,transposon_cols)