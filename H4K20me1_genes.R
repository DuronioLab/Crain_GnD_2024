library(tidyverse)
library(rtracklayer)

#import gtf file (downloaded from FlyBase)
gtf <- import("dmel-all-r6.55.gtf") %>% 
  as.data.frame()

#filter for protein-coding genes
gtf_protein_genes <- gtf %>% 
  filter(type == "gene") %>% 
  filter(!grepl("RNA",gene_symbol)) %>% 
  filter(!grepl("mir-",gene_symbol)) %>% 
  filter(!grepl("CR[0-9]",gene_symbol)) %>% 
  unique()

#export as bed file
make_bed_files <- function(gene_list_list) {
  for (gene_list_name in names(gene_list_list)) {
    gene_list <- as.data.frame(gene_list_list[[gene_list_name]])
    colnames(gene_list) <- sub("^[^.]*\\.", "", colnames(gene_list))
    gene_list.bed <- gtf_protein_genes %>% 
      dplyr::filter(gene_symbol %in% gene_list$gene_symbol) %>%
      dplyr::select(seqnames, start, end, gene_symbol, width, strand) %>% 
      dplyr::filter(!seqnames == "mitochondrion_genome")
    
    gene_list.bed$seqnames <- paste("chr", gene_list.bed$seqnames, sep = "")
    write_tsv(gene_list.bed, file = paste0(gene_list_name, ".bed"), col_names = FALSE)
  }
}

list_for_bed <- list("protein_genes_r6.55" = gtf_protein_genes
)
make_bed_files(list_for_bed)

#import h4k20me1 coverage gene overlap lists
k20me1_0.1_genes <- read_tsv("k20me1_genes_0.1.bed",col_names = FALSE) %>% 
  dplyr::select(X4) %>% 
  unique() %>% 
  na.omit()
colnames(k20me1_0.1_genes) <- "gene_symbol"

k20me1_0.25_genes <- read_tsv("k20me1_genes_0.25.bed",col_names = FALSE) %>% 
  dplyr::select(X4) %>% 
  unique() %>% 
  na.omit()
colnames(k20me1_0.25_genes) <- "gene_symbol"

k20me1_0.5_genes <- read_tsv("k20me1_genes_0.5.bed",col_names = FALSE) %>% 
  dplyr::select(X4) %>% 
  unique() %>% 
  na.omit()
colnames(k20me1_0.5_genes) <- "gene_symbol"

k20me1_0.75_genes <- read_tsv("k20me1_genes_0.75.bed",col_names = FALSE) %>% 
  dplyr::select(X4) %>% 
  unique() %>% 
  na.omit()
colnames(k20me1_0.75_genes) <- "gene_symbol"

k20me1_any_genes <- read_tsv("k20me1_genes_anyOverlap.bed",col_names = FALSE) %>% 
  dplyr::select(X4) %>% 
  unique() %>% 
  na.omit()
colnames(k20me1_any_genes) <- "gene_symbol"

##FIGURE 1C Nested venn diagram indicating number of genes with H4K20me1 peak overlap
k20me1_0.1_only <- k20me1_0.1_genes %>% 
  filter(!gene_symbol %in% k20me1_0.25_genes$gene_symbol)

k20me1_0.25_only <- k20me1_0.25_genes %>% 
  filter(!gene_symbol %in% k20me1_0.5_genes$gene_symbol)

k20me1_0.5_only <- k20me1_0.5_genes %>%
  filter(!gene_symbol %in% k20me1_0.75_genes$gene_symbol)

low_k20me1_genes <- gtf_protein_genes %>% 
  filter(!gene_symbol %in% k20me1_0.5_genes$gene_symbol)

no_k20me1_genes <- gtf_protein_genes %>% 
  filter(!gene_symbol %in% k20me1_0.1_genes$gene_symbol)

list_for_bed <- list("k20me1_genes_0.75_only" = k20me1_0.75_genes,
                     "k20me1_genes_0.5_only" = k20me1_0.5_only,
                     "k20me1_genes_0.25_only" = k20me1_0.25_only,
                     "k20me1_genes_0.1_only" = k20me1_0.1_only,
                     "low_k20me1_genes" = low_k20me1_genes,
                     "no_k20me1_genes" = no_k20me1_genes
)
make_bed_files(list_for_bed)
