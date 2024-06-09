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
