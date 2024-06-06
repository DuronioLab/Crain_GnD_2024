# Crain_GnD_2024
## Author: Aaron Crain

**Table of Contents**
- [Introduction](#introduction)
- [Quick Start](#quick-start)
- [Call H4K20me1 peaks](#call-H4K20me1-peaks)
    - [Required files](#required-files)
        - H4K20me1 CUT&RUN fastqs (GSE268819)
              - GSM8299933
              - GSM8299934
              - GSM8299935
              - GSM8299936
              - GSM8299937
              - GSM8299938
        - Run (Link to snakemake CUT&RUN pipeline)
        - [call_peaks_H4K20me1_CnR.R](#call_peaks_H4K20me1_CnR.R)
- [Annotate H4K20me1 peaks](#annotate_H4K20me1_peaks)
    - [Make upsetplot](#make_upset_plot)
        - [Required files](#required-files)
            - H4K20me1.vs.no_primary.peaks.txt (GSE268819)
            - [make_upsetplot.R](#make_upset_plot.R)
    - [Calcluate H4K20me1 peak gene overlap](#calculate_H4K20me1_peak_gene_overlap)
        - [Required files](#required-files)
            - H4K20me1.vs.no_primary.peaks.txt (GSE268819)
            - bedtools_intersect.sh
- [Process H4K20me1 wing disc CUT&RUN and whole larvae ChIP-seq](#Process_H4K20me1_wing_disc_CUT&RUN_and_whole_larvae_ChIP-seq)
    - [Required files](#required-files)
        - H4K20me1 ChIP-seq fastqs (GSE47254)
            - GSM1147213
            - GSM1147214
            - GSM1147215
            - GSM1147216
- [Make H4K20me1 gene overlap heatmap]
    - [Required files](#required-files)
        - gene percentage overlap bed files generated in [Calcluate H4K20me1 peak gene overlap](#calculate_H4K20me1_peak_gene_overlap)
        - processed bigwigs from [Process H4K20me1 wing disc CUT&RUN and whole larvae ChIP-seq](#Process_H4K20me1_wing_disc_CUT&RUN_and_whole_larvae_ChIP-seq)
        - deeptools computeMatrix and plotHeatmap.sh
- [H4K20me1 gene expression correlation](#H4K20me1_gene_expression_correlation)
    - [Required files](#required-files)
        - Salmon_protein_coding_index
        - make_Salmon_scripts
        - Oregon-R whole larvae RNA-seq fastqs (GSE268821)
        - yw wing disc RNA-seq fastqs (GSE141632)
            - GSM4210275
            - GSM4210276
            - GSM4210277
        - sample_sheet_wt_RNA-seq.txt
        - gene percentage overlap bed files generated in [Calcluate H4K20me1 peak gene overlap](#calculate_H4K20me1_peak_gene_overlap)
        - [wt_H4K20me1_gene_expression_correlation.R](#wt_H4K20me1_gene_expression_correlation.R)
  - [Call GFP-L(3)mbt peaks](#call-H4K20me1-peaks)
    - [Required files](#required-files)
        - GFP CUT&RUN fastqs (GSE268820)
              - GFP-L3mbt_rep1_L001_R1.fastq.gz GFP-L3mbt_rep1_L001_R2.fastq.gz GFP-L3mbt_rep1_L002_R1.fastq.gz GFP-L3mbt_rep1_L002_R2.fastq.gz
              - GFP-L3mbt_rep2_L001_R1.fastq.gz GFP-L3mbt_rep2_L001_R2.fastq.gz GFP-L3mbt_rep2_L002_R1.fastq.gz GFP-L3mbt_rep2_L002_R2.fastq.gz
              - OregonR_rep1_L001_R1.fastq.gz OregonR_rep1_L001_R2.fastq.gz OregonR_rep1_L002_R1.fastq.gz OregonR_rep1_L002_R2.fastq.gz
              - OregonR_rep2_L001_R1.fastq.gz OregonR_rep2_L001_R2.fastq.gz OregonR_rep2_L002_R1.fastq.gz OregonR_rep2_L002_R2.fastq.gz
        - Link to snakemake CUT&RUN pipeline
        - [call_peaks_GFP_CnR.R](#call_peaks_GFP_CnR.R)
- [Optional files](#optional-files)
- [Required directory structure](#required-directory-structure)
- [Expected Output](#expected-output)
- [Acknowledgements](#acknowledgements)

## Introduction:
This repository is provides code to reproduce results in Crain et al. Genes & Development 2024. 
Each Rscript (*.R) or shell script (*.sh) can be run independently, but some are dependent on outputs from previous scripts. Thes dependencies are noted throughout.

## Quick Start:
Download data from GEO
- GSE268819 - H4K20me1 CUT&RUN
- GSE268820 - GFP-L(3)mbt CUT&RUN
- GSE268821 - Set8, H4K20, and l(3)mbt mutant RNA-seq

## Call H4K20me1 peaks



