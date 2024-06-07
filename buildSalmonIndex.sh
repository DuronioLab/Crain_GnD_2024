#!/bin/bash 

module load salmon

cat dmel-all-transcript-r6.55.fasta.gz dmel-all-chromosome-r6.55.fasta.gz > gentrome_protein_genes_plus.fa.gz

salmon index -t gentrome_protein_genes.fa.gz -d decoys.txt -p 12 -i dm6_r6.55_protein_genes_salmon_index -k 31