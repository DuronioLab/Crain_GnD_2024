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
        - [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline)
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
            - bedtools intersect -a genes_r6.55.bed -b H4K20me1.vs.no_primary.peaks.bed | sort -u -k4,4 > k20me1_genes_anyOverlap_unique.bed
            - bedtools intersect -a genes_r6.55.bed -b H4K20me1.vs.no_primary.peaks.bed -f 0.1 | sort -u -k4,4 > k20me1_genes_0.1overlap_unique.bed
            - bedtools intersect -a genes_r6.55.bed -b H4K20me1.vs.no_primary.peaks.bed -f 0.25 | sort -u -k4,4 > k20me1_genes_0.25overlap_unique.bed
            - bedtools intersect -a genes_r6.55.bed -b H4K20me1.vs.no_primary.peaks.bed -f 0.50 | sort -u -k4,4 > k20me1_genes_0.50overlap_unique.bed
            - bedtools intersect -a genes_r6.55.bed -b H4K20me1.vs.no_primary.peaks.bed -f 0.75 | sort -u -k4,4 > k20me1_genes_0.75overlap_unique.bed
- [Process H4K20me1 wing disc CUT&RUN and whole larvae ChIP-seq](#Process_H4K20me1_wing_disc_CUT&RUN_and_whole_larvae_ChIP-seq)
    - [Required files](#required-files)
        - H4K20me1 ChIP-seq fastqs (GSE47254)
            - GSM1147213
            - GSM1147214
            - GSM1147215
            - GSM1147216
        - [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline)
        - all reps of RPGC normalized Oregon-R_H4K20me1, Oregon-R_no_primary, input_H4K20me1, ChIP_H4K20me1
        - In cutNrun-pipeline/BigWig/
            - deeptools bigwigAverage -b OregonR_H4K20me1_rep1_all_Frags_rpgcNorm.bw OregonR_H4K20me1_rep2_all_Frags_rpgcNorm.bw OregonR_H4K20me1_rep3_all_Frags_rpgcNorm.bw -bs 1 -o OregonR_H4K20me1_all_Frags_rpgcNorm_allReps_avg.bw
            - deeptools bigwigAverage -b OregonR_no_primary_rep1_all_Frags_rpgcNorm.bw OregonR_no_primary_rep2_all_Frags_rpgcNorm.bw OregonR_no_primary_rep3_all_Frags_rpgcNorm.bw -bs 1 -o OregonR_no_primary_all_Frags_rpgcNorm_allReps_avg.bw
            - deeptools bigwigAverage -b ChIP_H4K20me1_rep1_all_Frags_rpgcNorm.bw ChIP_H4K20me1_rep2_all_Frags_rpgcNorm.bw -bs 1 -o ChIP_H4K20me1_all_Frags_rpgcNorm_allReps_avg.bw
            - deeptools bigwigAverage -b input_H4K20me1_rep1_all_Frags_rpgcNorm.bw input_H4K20me1_rep2_all_Frags_rpgcNorm.bw -bs 1 -o input_H4K20me1_all_Frags_rpgcNorm_allReps_avg.bw
        - deeptools bigwigCompare -b1 OregonR_H4K20me1_all_Frags_rpgcNorm_allReps_avg.bw -b2 OregonR_no_primary_all_Frags_rpgcNorm_allReps_avg.bw -bs 1 -o OregonR_H4K20me1_all_Frags_rpgcNorm_allReps_avg_ratioCtrl.bw
        - deeptools bigwigCompare -b1 ChIP_H4K20me1_all_Frags_rpgcNorm_allReps_avg.bw -b2 input_H4K20me1_all_Frags_rpgcNorm_allReps_avg.bw -bs 1 -o ChIP_H4K20me1_all_Frags_rpgcNorm_allReps_avg_ratioCtrl.bw
        - zNorm.R
- [Make H4K20me1 gene overlap heatmap]
    - [Required files](#required-files)
        - gene percentage overlap bed files generated in [Calcluate H4K20me1 peak gene overlap](#calculate_H4K20me1_peak_gene_overlap)
        - processed bigwigs from [Process H4K20me1 wing disc CUT&RUN and whole larvae ChIP-seq](#Process_H4K20me1_wing_disc_CUT&RUN_and_whole_larvae_ChIP-seq)
        - deeptools computeMatrix scale-regions --regionsFileName 'k20me1_genes_0.75.bed' 'k20me1_genes_0.5.bed' 'k20me1_genes_0.25.bed' 'k20me1_genes_0.1.bed' 'nok20me1_genes.bed'  --scoreFileName 'OregonR_k20me1_rpgc_avg_allReps_ratioCtrl_zNorm.bw' 'h4k20me1_modEncode_rpgc_avg_allReps_ratioCtrl_zNorm.bw'  --samplesLabel 'Oregon-R wing disc CUT&RUN' 'Oregon-R whole larvae ChIP-seq'  --regionBodyLength 1000 --beforeRegionStartLength 1000 --afterRegionStartLength 1000  --unscaled5prime 0 --unscaled3prime 0 --sortRegions 'keep' --sortUsing 'mean' --averageTypeBins 'mean'  --missingDataAsZero --binSize 50 --transcriptID transcript --exonID exon --transcript_id_designator transcript_id
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
    - [Required files](#required-files)
        - H4K20me1 CUT&RUN fastqs (GSE268819)
            - GSM8299933-GSM8299959
        - H4K20me1_SRPMC_scaling_factors.txt
        - [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline)
        - allFrags.bed output from snakemake pipeline
        - spikeNorm_SRPMC.sh
- [Make H4K20me1 spike normalized heatmaps](Make_H4K20me1_spike_normalized_heatmaps)
    - [Required files](#required-files)
        - H4K20me1 spike-normalized bigwigs (GSE268819)
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
        - [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline)
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



