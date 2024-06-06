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
        - Run (link for snakemake CUT&RUN pipeline)
        - all reps of RPGC normalized Oregon-R_H4K20me1, Oregon-R_no_primary, input_H4K20me1, ChiP_H4K20me1
        - process_rpgc_bw.sh
- [Make H4K20me1 gene overlap heatmap]
    - [Required files](#required-files)
        - gene percentage overlap bed files generated in [Calcluate H4K20me1 peak gene overlap](#calculate_H4K20me1_peak_gene_overlap)
        - processed bigwigs from [Process H4K20me1 wing disc CUT&RUN and whole larvae ChIP-seq](#Process_H4K20me1_wing_disc_CUT&RUN_and_whole_larvae_ChIP-seq)
        - deeptools computeMatrix and plotHeatmap.sh
- [H4K20me1 gene expression correlation](#H4K20me1_gene_expression_correlation)
    - [Required files](#required-files)
        - Oregon-R whole larvae RNA-seq fastqs (GSE268821)
        - yw wing disc RNA-seq fastqs (GSE141632)
            - GSM4210275
            - GSM4210276
            - GSM4210277
        - trim and fastq check
        - Salmon_protein_coding_index
        - make_Salmon_scripts
        - sample_sheet_wt_RNA-seq.txt
        - gene percentage overlap bed files generated in [Calcluate H4K20me1 peak gene overlap](#calculate_H4K20me1_peak_gene_overlap)
        - [wt_H4K20me1_gene_expression_correlation.R](#wt_H4K20me1_gene_expression_correlation.R)
- [Spike-in normalization for CUT&RUN](#Spike-in_normalization_for_CUT&RUN)
- [Make H4K20me1 spike normalized heatmaps](Make_H4K20me1_spike_normalized_heatmaps)
    - [Required files](#required-files)
        - H4K20me1 spike-normalized bigwigs
            - OregonR_spikeNorm_K20me1_allReps_avg.bw
            - Set8null_spikeNorm_K20me1_allReps_avg.bw
            - Set8wt_spikeNorm_K20me1_allReps_avg.bw
            - Set8rg_spikeNorm_K20me1_allReps_avg.bw
            - HWT_spikeNorm_K20me1_allReps_avg.bw
            - K20A_spikeNorm_K20me1_allReps_avg.bw
            - K20R_spikeNorm_K20me1_allReps_avg.bw
        - computeMatrix and plotHeatmap.sh
- [Differential expression analysis Set8, H4K20, and l(3)mbt mutants](#differential_expression_analysis)
    - [Required files](#required-files)
        - Whole larvae RNA-seq fastqs (GSE268821)
            - GSM8299968-8300014
        - trim and fastq check
        - Salmon_protein_coding_index
        - make_Salmon_scripts
        - sample_sheet_RNA-seq_3wl.txt
        - [differential_expression_analysis.R](#differential_expression_analysis.R)
- [Call GFP-L(3)mbt peaks](#call-H4K20me1-peaks)
    - [Required files](#required-files)
        - GFP CUT&RUN fastqs (GSE268820)
              - GSM8299960
              - GSM8299961
              - GSM8299962
              - GSM8299963
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



