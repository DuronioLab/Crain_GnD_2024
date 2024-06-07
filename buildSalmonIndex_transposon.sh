#!/bin/bash 

module load salmon

cat dmel-all-transcript-r6.55.fasta.gz dmel-all-miRNA-r6.55.fasta.gz dmel-all-miscRNA-r6.55.fasta.gz dmel-all-ncRNA-r6.55.fasta.gz dmel-all-tRNA-r6.55.fasta.gz dmel-all-pseudogene-r6.55.fasta.gz transposon_sequence_set.fa.gz piRNA_clusters.fa.gz dmel-all-chromosome-r6.55.fasta.gz > gentrome_all_genes_plus_transposons_piRNAs.fa.gz

salmon index -t gentrome_all_genes_plus_transposons_piRNAs.fa.gz -d decoys.txt -p 12 -i dm6_r6.55_all_genes_plus_transposons_piRNAs_salmon_index -k 31