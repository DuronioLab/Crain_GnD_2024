# Crain_GnD_2024
## Author: Aaron Crain

**Table of Contents**
- [Introduction](#introduction)
- [Quick Start](#quick-start)
- [Call H4K20me1 peaks](#call-H4K20me1-peaks)
    - [Required files](#required-files)
        - [H4K20me1 CUT&RUN fastq files (GSE268819)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE268819)
        - [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline)
        - [H4K20me1_CnR_sample_sheet.txt](#H4K20me1_CnR_sample_sheet.txt)
        - [call_peaks_H4K20me1_CnR.R](#call_peaks_H4K20me1_CnR.R)
    - [Run code](#run_code)
    - [Expected output](#expected_output)
- [Annotate H4K20me1 peaks](#Annotate_H4K20me1_peaks)
    - [Make upsetplot](#make_upset_plot)
        - [Required files](#required-files)
            - [H4K20me1.vs.no_primary.peaks.bed (GSE268819)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE268819)
            - [make_upsetplot.R](#make_upset_plot.R)
        - [Run code](#run_code)
        - [Expected output](#expected_output)
    - [Calcluate H4K20me1 peak gene overlap](#calculate_H4K20me1_peak_gene_overlap)
        - [Required files](#required-files)
            - [H4K20me1.vs.no_primary.peaks.bed (GSE268819)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE268819)
            - [make_gene_bed_from_gtf.R](#make_gene_bed_from_gtf.R)
            - [H4K20me1_genes.R](#H4K20me1_genes.R)
        - [Run code](#run_code)
        - [Expected outputs](#expected_outputs)
- [Process H4K20me1 wing disc CUT&RUN and whole larvae ChIP-seq](#Process_H4K20me1_wing_disc_CUT&RUN_and_whole_larvae_ChIP-seq)
    - [Required files](#required-files)
        - [H4K20me1 ChIP-seq fastqs (GSE47254)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47254)
        - [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline)
        - [zNorm.R](https://github.com/snystrom/cutNrun-pipeline/blob/master/scripts/zNorm.r)
    - [Run code](#run_code)
    - [Expected outputs](#expected_outputs)
- [Make H4K20me1 gene overlap heatmap](#Make_H4K20me1_gene_overlap_heatmap)
    - [Required files](#required-files)
        - gene percentage overlap bed files generated in [Calcluate H4K20me1 peak gene overlap](#calculate_H4K20me1_peak_gene_overlap)
        - processed bigwigs from [Process H4K20me1 wing disc CUT&RUN and whole larvae ChIP-seq](#Process_H4K20me1_wing_disc_CUT&RUN_and_whole_larvae_ChIP-seq)
    - [Run code](#run_code)
    - [Expected outputs](#expected_outputs)
        - Plot in FIGURE 1D
- [H4K20me1 gene expression correlation](#H4K20me1_gene_expression_correlation)
    - [Required files](#required-files)
        - [Oregon-R whole larvae RNA-seq fastqs (GSE268821)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE268821)
        - [yw wing disc RNA-seq fastqs (GSE141632)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141632)
        - [RNA-seq pipeline](https://github.com/DuronioLab/RNAseq-pipeline)
        - [build_salmon_index.sh](#build_salmon_index)
        - [dmel-all-transcript-r6.55.fasta.gz](#dmel-all-transcript-r6.55.fasta.gz)
        - [dmel-all-chromosome-r6.55.fasta.gz](#dmel-all-chromosome-r6.55.fasta.gz)
        - [decoys.txt](#decoys.txt)
        - [run_salmon.py](#run_salmon.py)
        - [sample_sheet_RNA-seq_wild-type.txt](#sample_sheet_RNA-seq_wild-type.txt)
        - gene percentage overlap bed files generated in [Calcluate H4K20me1 peak gene overlap](#calculate_H4K20me1_peak_gene_overlap)
        - [wt_H4K20me1_gene_expression_correlation.R](#wt_H4K20me1_gene_expression_correlation.R)
    - [Run code](#run_code)
        - [Build salmon index](#build_salmon_index)
        - [Run Salmon](#run_salmon)
    - [Expected outputs](#expected_outputs)
        - Plot in FIGURE 1E
- [Spike-in normalization for H4K20me1 CUT&RUN](#Spike-in_normalization_for_H4K20me1_CUT&RUN)
    - [Required files](#required-files)
        - [H4K20me1 CUT&RUN fastq files (GSE268819)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE268819)
        - [H4K20me1_CnR_sample_sheet](#H4K20me1_CnR_sample_sheet)
        - [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline)
        - [dm6.chrom.sizes](#dm6.chrom.sizes)
        - [h4k20me1_spikeNorm_SRPMC.sh](#h4k20me1_spikeNorm_SRPMC.sh)
    - [Run code](#run_code)
    - [Expected outputs](#expected_outputs)
        - Spike-in normalized CUT&RUN bigwigs
- [Make H4K20me1 spike-in normalized heatmaps](Make_H4K20me1_spike-in_normalized_heatmaps)
    - [Required files](#required-files)
        - [H4K20me1 spike-normalized bigwigs (GSE268819)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE268819)
    - [Run code](#run_code)
    - [Expected outputs](#expected_outputs)
        - Plot in FIGURE 2A
- [Differential expression analysis Set8, H4K20, and l(3)mbt mutants](#differential_expression_analysis)
    - [Required files](#required-files)
        - [Whole larvae RNA-seq fastq files (GSE268821)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE268821)
        - [RNA-seq pipeline](https://github.com/DuronioLab/RNAseq-pipeline)
        - [run_salmon.py](#run_salmon.py)
        - [sample_sheet_RNA-seq_3wl.txt](#sample_sheet_RNA-seq_3wl.txt)
        - [differential_expression_analysis.R](#differential_expression_analysis.R)
    - [Run code](#run_code)
    - [Expected outputs](#expected_outputs)
        - Plots in FIGURES 3 and 5
- [Call GFP-L(3)mbt peaks](#call_L(3)mbt_peaks)
    - [Required files](#required-files)
        - [GFP CUT&RUN fastq files (GSE268820)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE268820)
        - [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline)
        - [GFP-l3mbt_CnR_sample_sheet](#GFP-l3mbt_CnR_sample_sheet)
        - [call_peaks_GFP_CnR.R](#call_peaks_GFP_CnR.R)
    - [Run code](#run_code)
    - [Expected outputs](#expected_outputs)
        - GFP-L3mbt.vs.OregonR.peaks.bed
- [Process wild-type GFP-L(3)mbt CUT&RUN and L(3)mbt ChIP-seq](#Process_wild-type_GFP_CUT&RUN)
    - [Required files](#required-files)
        - [GFP CUT&RUN fastq files (GSE268820)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE268820)
        - [L(3)mbt ChIP-seq fastq files (GSE29206)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29206)
        - [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline)
        - [zNorm.R](https://github.com/snystrom/cutNrun-pipeline/blob/master/scripts/zNorm.r)
    - [Run code](#run_code)
    - [Expected outputs](#expected_outputs)
- [Make H4K20me1 and L(3)mbt heatmaps](#Make_H4K20me1_and_L(3)mbt_heatmap)
    - [Required files](#required-files)
        - k20me1_genes_0.5.bed from [Calcluate H4K20me1 peak gene overlap](#calculate_H4K20me1_peak_gene_overlap)
        - processed bigwigs from [Process H4K20me1 wing disc CUT&RUN and whole larvae ChIP-seq](#Process_H4K20me1_wing_disc_CUT&RUN_and_whole_larvae_ChIP-seq) and [Process wild-type GFP-L(3)mbt CUT&RUN and L(3)mbt ChIP-seq](#Process_wild-type_GFP_CUT&RUN)
    - [Run code](#run_code)
    - [Expected output](#expected_output)
        - Plots in FIGURES 6C and D
- [Spike-in normalization for GFP CUT&RUN](#Spike-in_normalization_for_GFP_CUT&RUN)
    - [Required files](#required-files)
        - [GFP CUT&RUN fastq files (GSE268820)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE268820)
        - [GFP-l3mbt_CnR_sample_sheet](#GFP-l3mbt_CnR_sample_sheet)
        - [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline)
        - [gfp_spikeNorm_SRPMC.sh](#gfp_spikeNorm_SRPMC.sh)
    - [Run code](#run_code)
    - [Expected output](#expected_output)
        - Spike-in normalized bigwigs
- [Make H4K20me1 and L(3)mbt heatmap](#Make_H4K20me1_and_L(3)mbt_heatmap)
    - [Required files](#required-files)
        - [GFP-L3mbt.vs.OregonR.peaks.bed](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE268820)
        - spike-in normalized bigwigs from [Spike-in normalization for GFP CUT&RUN](#Spike-in_normalization_for_GFP_CUT&RUN)
    - [Run code](#run_code)
    - [Expected output](#expected_output)
        - Plots in FIGURE 6F
- [Acknowledgements](#acknowledgements)
    - Markus, Jeanne-Marie, Spencer, Matt

## Introduction:
This repository is provides code to reproduce results in Crain et al. Genes & Development 2024. 
Each Rscript section can be run independently, but some are dependent on outputs from previous sections. These dependencies are noted throughout.

## Quick Start:
Download data from GEO
- GSE268819 - H4K20me1 CUT&RUN
- GSE268820 - GFP-L(3)mbt CUT&RUN
- GSE268821 - Set8, H4K20, and l(3)mbt mutant RNA-seq

## Call H4K20me1 peaks
### Required files
#### H4K20me1 CUT&RUN fastq files (GSE268819)
GSM8299933
GSM8299934
GSM8299935
GSM8299936
GSM8299937
GSM8299938
These files contain sequencing reads from Oregon-R and Oregon-R no primary control can be downloaded from GEO
#### CUT&RUN pipeline
For processing CUT&RUN sequencing files. See documentation here for more information https://github.com/snystrom/cutNrun-pipeline 
#### call_peaks_H4K20me1_CnR.R
This R script takes bam files from the CUT&RUN pipeline as input. Aligned reads are binned into 150bp windows that overlap by 50bp using the csaw package and significant windows (OregonR_H4K20me1 vs. OregonR_no_primary) are determined by edgeR. Significant windows within 1kb are merged to make final H4K20me1 peaks.
### Run code
#### Run CUT&RUN pipeline
#### Move H4K20me1_CnR_sample_sheet.txt, call_peaks_H4K20me1_CnR.R, and OregonR and OregonR_no_primary .bam and .bai files into the same working directory and/or ensure the Path_to_bam matches from H4K20me1_CnR_sample_sheet.txt matches the .bam and .bai file names and paths.
#### Run R script
### Expected output
H4K20me1.vs.no_primary.peaks.bed - a bed file with H4K20me1 peaks. 

## Annotate H4K20me1 peaks
### Make upsetplot
#### Required files
##### H4K20me1.vs.no_primary.peaks.bed (GSE268819)
This is the peak file generated in the previous section and is available to download directly from GEO.
##### make_upsetplot.R
This Rscript uses ChIPseeker to annotate the peak file and produce an upset plot depicting the genomic annotations of all H4K20me1 peaks. 
#### Run code
##### Download or move H4K20me1.vs.no_primary.peaks.bed into working directory.
##### Run Rscript
#### Expected output
##### Plot in FIGURE 1B

### Calculate H4K20me1 peak gene overlap
#### Required files
##### H4K20me1.vs.no_primary.peaks.bed (GSE268819)
This is the peak file generated in the previous section and is available to download directly from GEO.
##### make_gene_bed_from_gtf.R
Rscript that takes a .gtf file, extracts protein-coding genes, and makes the protein_genes_r6.55.bed file
##### H4K20me1_genes.R
Rscript that takes H4K20me1 gene overlap .bed files and generates the nok20me1.bed (<10% H4K20me1 overlap),  low_k20me1_genes.bed (< 50% H4K20me1 overlap), and non-overlapping gene .bed files for each category
#### Run code
Download dmel-all-r6.55.gtf from FlyBase. Run make_gene_bed_from_gtf.R to generate protein_genes_r6.55.bed file.
Then use bedtools to determine overlap of H4K20me1 peaks with genes.
```
module load bedtools
bedtools intersect -a protein_genes_r6.55.bed -b H4K20me1.vs.no_primary.peaks.bed | sort -u -k4,4 > k20me1_genes_anyOverlap.bed
bedtools intersect -a protein_genes_r6.55.bed -b H4K20me1.vs.no_primary.peaks.bed -f 0.1 | sort -u -k4,4 > k20me1_genes_0.1.bed
bedtools intersect -a protein_genes_r6.55.bed -b H4K20me1.vs.no_primary.peaks.bed -f 0.25 | sort -u -k4,4 > k20me1_genes_0.25.bed
bedtools intersect -a protein_genes_r6.55.bed -b H4K20me1.vs.no_primary.peaks.bed -f 0.50 | sort -u -k4,4 > k20me1_genes_0.5.bed
bedtools intersect -a protein_genes_r6.55.bed -b H4K20me1.vs.no_primary.peaks.bed -f 0.75 | sort -u -k4,4 > k20me1_genes_0.75.bed
```
Run H4K20me1_genes.R to generate k20me1_genes_0.75_only, k20me1_genes_0.5_only, k20me1_genes_0.25_only, k20me1_genes_0.1_only, low_k20me1_genes, and no_k20me1_genes .bed files.

#### Expected outputs
##### protein_genes_r6.55.bed
##### k20me1_genes_anyOverlap.bed - genes with > 1bp H4K20me1 overlap
##### k20me1_genes_0.1.bed - genes with > 10% H4K20me1 overlap
##### k20me1_genes_0.25.bed - genes with > 525% H4K20me1 overlap
##### k20me1_genes_0.5.bed - genes with > 50% H4K20me1 overlap
##### k20me1_genes_0.75.bed - genes with > 75% H4K20me1 overlap

##### low_k20me1_genes.bed - genes with < 50% H4K20me1 overlap
##### k20me1_genes_0.75_only.bed - genes with > 75% H4K20me1 overlap
##### k20me1_genes_0.5_only.bed - genes with > 50% and < 75% H4K20me1 overlap
##### k20me1_genes_0.25_only.bed - genes with > 25% and < 50% H4K20me1 overlap
##### k20me1_genes_0.1_only.bed - genes with > 10% and < 25% H4K20me1 overlap
##### no_k20me1_genes.bed - genes with < 10% H4K20me1 overlap

##### Plot in FIGURE 1C

## Process H4K20me1 wing disc CUT&RUN and whole larvae ChIP-Seq
### Required files
#### H4K20me1 ChIP-seq fastqs (GSE47254)
GSM1147213
GSM1147214
GSM1147215
GSM1147216
#### CUT&RUN pipeline
For processing CUT&RUN sequencing files. See documentation here for more information [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline) 
#### zNorm.R
For z-Normalizing bigwig files. See [zNorm.R](https://github.com/snystrom/cutNrun-pipeline/blob/master/scripts/zNorm.r).
### Run code
```
module load deeptools
bigwigAverage -b OregonR_H4K20me1_rep1_allFrags_rpgcNorm.bw OregonR_H4K20me1_rep2_allFrags_rpgcNorm.bw OregonR_H4K20me1_rep3_allFrags_rpgcNorm.bw --scaleFactors '1' --binSize 1 -o OregonR_H4K20me1_allFrags_rpgcNorm_allReps_avg.bw --outFileFormat 'bigwig'
bigwigAverage -b OregonR_no_primary_rep1_allFrags_rpgcNorm.bw OregonR_no_primary_rep2_allFrags_rpgcNorm.bw OregonR_no_primary_rep3_allFrags_rpgcNorm.bw --scaleFactors '1' --binSize 1 -o OregonR_no_primary_allFrags_rpgcNorm_allReps_avg.bw --outFileFormat 'bigwig'
bigwigAverage -b ChIP_H4K20me1_rep1_allFrags_rpgcNorm.bw ChIP_H4K20me1_rep2_allFrags_rpgcNorm.bw --scaleFactors '1' --binSize 1 -o ChIP_H4K20me1_allFrags_rpgcNorm_allReps_avg.bw --outFileFormat 'bigwig'
bigwigAverage -b input_H4K20me1_rep1_allFrags_rpgcNorm.bw input_H4K20me1_rep2_allFrags_rpgcNorm.bw --scaleFactors '1' --binSize 1 -o input_H4K20me1_allFrags_rpgcNorm_allReps_avg.bw --outFileFormat 'bigwig'
bigwigCompare -b1 OregonR_H4K20me1_allFrags_rpgcNorm_allReps_avg.bw -b2 OregonR_no_primary_allFrags_rpgcNorm_allReps_avg.bw --operation ratio --pseudocount 1 1      --scaleFactors '1.0:1.0' --binSize 1 -o OregonR_H4K20me1_allFrags_rpgcNorm_allReps_avg_ratioCtrl.bw --outFileFormat 'bigwig'
bigwigCompare -b1 ChIP_H4K20me1_allFrags_rpgcNorm_allReps_avg.bw -b2 input_H4K20me1_allFrags_rpgcNorm_allReps_avg.bw --operation ratio --pseudocount 1 1 --scaleFactors '1.0:1.0' --binSize 1 -o ChIP_H4K20me1_allFrags_rpgcNorm_allReps_avg_ratioCtrl.bw --outFileFormat 'bigwig'
```
Pass OregonR_H4K20me1_allFrags_rpgcNorm_allReps_avg_ratioCtrl.bw and ChIP_H4K20me1_allFrags_rpgcNorm_allReps_avg_ratioCtrl.bw individually to [zNorm.R](https://github.com/snystrom/cutNrun-pipeline/blob/master/scripts/zNorm.r)
### Expected outputs
#### From [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline) cutNrun-pipeline/BigWig/
    - OregonR_H4K20me1_rep1_allFrags_rpgcNorm.bw
    - OregonR_H4K20me1_rep2_allFrags_rpgcNorm.bw
    - OregonR_H4K20me1_rep3_allFrags_rpgcNorm.bw
    - OregonR_no_primary_rep1_allFrags_rpgcNorm.bw
    - OregonR_no_primary_rep2_allFrags_rpgcNorm.bw
    - OregonR_no_primary_rep3_allFrags_rpgcNorm.bw
    - ChIP_H4K20me1_rep1_allFrags_rpgcNorm.bw
    - ChIP_H4K20me1_rep2_allFrags_rpgcNorm.bw
    - input_H4K20me1_rep1_allFrags_rpgcNorm.bw
    - input_H4K20me1_rep2_allFrags_rpgcNorm.bw
#### From [zNorm.R](https://github.com/snystrom/cutNrun-pipeline/blob/master/scripts/zNorm.r)
    - OregonR_H4K20me1_allFrags_rpgcNorm_allReps_avg_ratioCtrl_zNorm.bw
    - ChIP_H4K20me1_allFrags_rpgcNorm_allReps_avg_ratioCtrl_zNorm.bw

## Make H4K20me1 gene overlap heatmap
### Required files
From [Calcluate H4K20me1 peak gene overlap](#calculate_H4K20me1_peak_gene_overlap)
    - nok20me1_genes.bed - genes with < 50% H4K20me1 overlap
    - k20me1_genes_0.75_only.bed - genes with > 75% H4K20me1 overlap
    - k20me1_genes_0.5_only.bed - genes with > 50% and < 75% H4K20me1 overlap
    - k20me1_genes_0.25_only.bed - genes with > 25% and < 50% H4K20me1 overlap
    - k20me1_genes_0.1_only.bed - genes with > 10% and < 25% H4K20me1 overlap
From [Process H4K20me1 wing disc CUT&RUN and whole larvae ChIP-seq](#Process_H4K20me1_wing_disc_CUT&RUN_and_whole_larvae_ChIP-seq)
    - OregonR_H4K20me1_allFrags_rpgcNorm_allReps_avg_ratioCtrl_zNorm.bw
    - ChIP_H4K20me1_allFrags_rpgcNorm_allReps_avg_ratioCtrl_zNorm.bw
### Run code
```
module load deeptools
computeMatrix scale-regions --regionsFileName k20me1_genes_0.75_only.bed k20me1_genes_0.5_only.bed k20me1_genes_0.25_only.bed k20me1_genes_0.1_only.bed no_k20me1_genes.bed  --scoreFileName OregonR_H4K20me1_allFrags_rpgcNorm_allReps_avg_ratioCtrl_zNorm.bw ChIP_H4K20me1_allFrags_rpgcNorm_allReps_avg_ratioCtrl_zNorm.bw  --samplesLabel 'Oregon-R wing disc CUT&RUN' 'Oregon-R whole larvae ChIP-seq'  --regionBodyLength 1000 --beforeRegionStartLength 1000 --afterRegionStartLength 1000  --unscaled5prime 0 --unscaled3prime 0 --sortRegions 'keep' --sortUsing 'mean' --averageTypeBins 'mean'  --missingDataAsZero --binSize 50 -o h4k20me1_gene_overlap_computeMatrix
```
```
plotHeatmap --matrixFile h4k20me1_gene_overlap_computeMatrix --outFileName h4k20me1_gene_overlap.pdf --plotFileFormat 'pdf' --dpi '200' --sortRegions 'descend' --sortUsing 'mean'  --averageTypeSummaryPlot 'mean'  --plotType 'lines'  --missingDataColor 'black'  --colorMap Greens  --alpha '1.0'   --sortUsingSamples 1  --heatmapWidth 7.5 --heatmapHeight 25.0  --whatToShow 'plot and heatmap'  --startLabel 'TSS' --endLabel 'TES'  --regionsLabel over0.75 over0.5 over0.25 over0.1 nok20me1   --legendLocation 'best'  --labelRotation '0'
```
### Expected outputs
#### Plot in FIGURE 1D

## H4K20me1 gene expression correlation
### Required files
#### Oregon-R whole larvae RNA-seq fastqs (GSE268821)
GSM8299933
GSM8299934
GSM8299935
#### yw wing disc RNA-seq fastqs (GSE141632)
GSM4210275
GSM4210276
GSM4210277
#### [RNA-seq pipeline](https://github.com/DuronioLab/RNAseq-pipeline)
#### buildSalmonIndex.sh
Shell script that builds the required index for executing Salmon 
#### dmel-all-transcript-r6.55.fasta.gz
Drosophila transcripts
#### dmel-all-chromosome-r6.55.fasta.gz
Drososphila genome
#### decoys.txt
Decoy sequences using the entire Drosophila genome (see https://salmon.readthedocs.io/en/latest/salmon.html#using-salmon for more details)
#### run_salmon.py
Python script that runs Salmon quantification
#### sample_sheet_RNA-seq_wild-type.txt
#### From [Calcluate H4K20me1 peak gene overlap](#calculate_H4K20me1_peak_gene_overlap)
k20me1_genes_0.1.bed - genes with > 10% H4K20me1 overlap
k20me1_genes_0.25.bed - genes with > 25% H4K20me1 overlap
k20me1_genes_0.5.bed - genes with > 50% H4K20me1 overlap
k20me1_genes_0.75.bed - genes with > 75% H4K20me1 overlap
#### [wt_H4K20me1_gene_expression_correlation.R](#wt_H4K20me1_gene_expression_correlation.R)
### Run code
#### Download .fastq files from GEO. 
#### Build protein-coding Salmon index 
Download dmel-all-transcript-r6.55.fasta.gz and dmel-all-chromosome-r6.55.fasta.gz from FlyBase and ensure those files, decoys.txt, and buildSalmonIndex.sh are in the same working directory. Run buildSalmonIndex.sh
#### Run RNA-seq pipeline (https://github.com/DuronioLab/RNAseq-pipeline) to perform quality control and trim adapters
#### Run Salmon
Ensure Salmon index, run_salmon.py, and sample_sheet_RNA-seq_wild-type.txt are in the same working directory. Update run_salmon with appropriate index name, sample sheet name, and path to trimmed fastq files. Run run_salmon.py. 
### Expected outputs
#### Plot in FIGURE 1E

## Spike-in normalization for H4K20me1 CUT&RUN
### Required files
#### H4K20me1 CUT&RUN fastqs (GSE268819)
GSM8299933-GSM8299959
#### H4K20me1_CnR_sample_sheet.txt
#### CUT&RUN pipeline
For processing CUT&RUN sequencing files. See documentation here for more information [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline) 
#### h4k20me1_spikeNorm_SRPMC.sh
### Run code
### Expected outputs
Spike-in normalized bigwigs

## Make H4K20me1 spike-in normalized heatmaps
### Required files
#### H4K20me1 spike-normalized bigwigs (GSE268819)
##### OregonR_spikeNorm_K20me1_allReps_avg.bw
##### Set8null_spikeNorm_K20me1_allReps_avg.bw
##### Set8wt_spikeNorm_K20me1_allReps_avg.bw
##### Set8rg_spikeNorm_K20me1_allReps_avg.bw
##### HWT_spikeNorm_K20me1_allReps_avg.bw
##### K20A_spikeNorm_K20me1_allReps_avg.bw
##### K20R_spikeNorm_K20me1_allReps_avg.bw
### Run code
```
module load deeptools
computeMatrix scale-regions --regionsFileName H4K20me1.vs.no_primary.peaks.bed --scoreFileName OregonR_spikeNorm_K20me1_allReps_avg.bw Set8null_spikeNorm_K20me1_allReps_avg.bw Set8wt_spikeNorm_K20me1_allReps_avg.bw Set8rg_spikeNorm_K20me1_allReps_avg.bw HWT_spikeNorm_K20me1_allReps_avg.bw K20A_spikeNorm_K20me1_allReps_avg.bw K20R_spikeNorm_K20me1_allReps_avg.bw --samplesLabel 'Oregon-R' 'Set8null' 'Set8wt' 'Set8rg' 'HWT' 'H4K20A' 'H4K20R' --regionBodyLength 200 --beforeRegionStartLength 1000 --afterRegionStartLength 1000  --unscaled5prime 0 --unscaled3prime 0 --sortRegions 'keep' --sortUsing 'mean' --averageTypeBins 'mean'  --missingDataAsZero --binSize 50 -o h4k20me1_spikeNorm_computeMatrix
```
```
plotHeatmap --matrixFile h4k20me1_spikeNorm_computeMatrix --outFileName h4k20me1_spikeNorm.pdf  --plotFileFormat 'pdf'    --dpi '200'  --sortRegions 'descend'   --sortUsing 'mean'  --averageTypeSummaryPlot 'mean'  --plotType 'lines'  --missingDataColor 'black'  --colorMap Greens  --alpha '1.0'   --sortUsingSamples 1  --xAxisLabel NA --yAxisLabel 'H4K20me1 peaks'  --heatmapWidth 7.5 --heatmapHeight 25.0  --whatToShow 'plot and heatmap'  --startLabel 'Peak_start' --endLabel 'Peak_end'  --samplesLabel or set8null set8wt set8rg hwt k20a k20r --regionsLabel k20me1_peaks   --legendLocation 'best'  --labelRotation '0'
```
### Expected output
#### Plot in FIGURE 2A

## Differential expression analysis Set8, H4K20, and l(3)mbt mutants
### Required files
#### Whole larvae RNA-seq fastq files (GSE268821)
GSM8299968-8300014
#### [RNA-seq pipeline](https://github.com/DuronioLab/RNAseq-pipeline)
#### Salmon protein coding genes index from [H4K20me1 gene expression correlation](#H4K20me1_gene_expression_correlation)
#### run_salmon.py
#### sample_sheet_RNA-seq_3wl.txt
#### [differential_expression_analysis.R](#differential_expression_analysis.R)
### Run code
### Expected outputs
#### Plots in FIGURES 3 and 5

## Call GFP-L(3)mbt peaks
### Required files
#### GFP CUT&RUN fastq files (GSE268820)
GSM8299960-8299963
#### [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline)
#### call_peaks_GFP_CnR.R
### Run code
### Expected outputs
#### GFP-L3mbt.vs.OregonR.peaks.bed

## Process wild-type GFP-L(3)mbt CUT&RUN and L(3)mbt ChIP-seq
### Required files
#### GFP-L(3)mbt CUT&RUN fastq files
GSM8299960 GSM8299961
#### L(3)mbt ChIP-seq fastq files (GSE29206)
GSM722523-722526
#### [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline)
#### [zNorm.R](https://github.com/snystrom/cutNrun-pipeline/blob/master/scripts/zNorm.r)
### Run code
```
deeptools bigwigAverage -b l3mbtGFP_rep1_all_Frags_rpgcNorm.bw l3mbtGFP_rep2_all_Frags_rpgcNorm.bw -bs 1 -o l3mbtGFP_all_Frags_rpgcNorm_allReps_avg.bw
deeptools bigwigAverage -b OregonR_rep1_all_Frags_rpgcNorm.bw OregonR_rep2_all_Frags_rpgcNorm.bw -bs 1 -o OregonR_all_Frags_rpgcNorm_allReps_avg.bw
deeptools bigwigAverage -b l3mbt_ChIP_rep1_all_Frags_rpgcNorm.bw l3mbt_ChIP_rep2_all_Frags_rpgcNorm.bw -bs 1 -o l3mbt_ChIP_all_Frags_rpgcNorm_allReps_avg.bw
deeptools bigwigAverage -b l3mbt_input_rep1_all_Frags_rpgcNorm.bw l3mbt_input_rep2_all_Frags_rpgcNorm.bw -bs 1 -o l3mbt_input_all_Frags_rpgcNorm_allReps_avg.bw
deeptools bigwigCompare -b1 l3mbtGFP_all_Frags_rpgcNorm_allReps_avg.bw -b2 OregonR_all_Frags_rpgcNorm_allReps_avg.bw -bs 1 -o l3mbtGFP_all_Frags_rpgcNorm_allReps_avg_ratioCtrl.bw
deeptools bigwigCompare -b1 l3mbt_ChIP_all_Frags_rpgcNorm_allReps_avg.bw -b2 l3mbt_input_all_Frags_rpgcNorm_allReps_avg.bw -bs 1 -o l3mbt_ChIP_all_Frags_rpgcNorm_allReps_avg_ratioCtrl.bw
```
Pass l3mbtGFP_all_Frags_rpgcNorm_allReps_avg_ratioCtrl.bw and l3mbt_ChIP_all_Frags_rpgcNorm_allReps_avg_ratioCtrl.bw individually to [zNorm.R](https://github.com/snystrom/cutNrun-pipeline/blob/master/scripts/zNorm.r)
### Expected outputs
#### From [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline) cutNrun-pipeline/BigWig/
    - l3mbtGFP_rep1_all_Frags_rpgcNorm.bw
    - l3mbtGFP_rep2_all_Frags_rpgcNorm.bw
    - OregonR_rep1_all_Frags_rpgcNorm.bw
    - OregonR_rep2_all_Frags_rpgcNorm.bw
    - l3mbt_ChIP_rep1_all_Frags_rpgcNorm.bw 
    - l3mbt_ChIP_rep2_all_Frags_rpgcNorm.bw
    - l3mbt_input_rep1_all_Frags_rpgcNorm.bw 
    - l3mbt_input_rep2_all_Frags_rpgcNorm.bw
#### From [zNorm.R](https://github.com/snystrom/cutNrun-pipeline/blob/master/scripts/zNorm.r)
    - l3mbtGFP_all_Frags_rpgcNorm_allReps_avg_ratioCtrl_zNorm.bw
    - l3mbt_ChIP_all_Frags_rpgcNorm_allReps_avg_ratioCtrl_zNorm.bw

## Make H4K20me1 and L(3)mbt heatmaps
### Required files
#### k20me1_genes_0.5.bed from [Calcluate H4K20me1 peak gene overlap](#calculate_H4K20me1_peak_gene_overlap)
#### processed bigwigs from [Process H4K20me1 wing disc CUT&RUN and whole larvae ChIP-seq](#Process_H4K20me1_wing_disc_CUT&RUN_and_whole_larvae_ChIP-seq) and [Process wild-type GFP-L(3)mbt CUT&RUN and L(3)mbt ChIP-seq](#Process_wild-type_GFP_CUT&RUN)
### Run code
```
deeptools computeMatrix scale-regions --regionsFileName 'k20me1_genes_0.5.bed' --scoreFileName 'OregonR_k20me1_rpgc_avg_allReps_ratioCtrl_zNorm.bw' 'l3mbtGFP_all_Frags_rpgcNorm_allReps_avg_ratioCtrl.bw' l3mbt_ChIP_all_Frags_rpgcNorm_allReps_avg_ratioCtrl.bw --samplesLabel 'H4K20me1' 'GFP-L(3)mbt' 'L(3)mbt_Richter_et_al' --regionBodyLength 1000 --beforeRegionStartLength 1000 --afterRegionStartLength 1000  --unscaled5prime 0 --unscaled3prime 0 --sortRegions 'keep' --sortUsing 'mean' --averageTypeBins 'mean'  --missingDataAsZero --binSize 50
deeptools plotProfile
```
```
deeptools computeMatrix reference-point --regionsFileName 'genes_r6.55.bed' --scoreFileName 'OregonR_k20me1_rpgc_avg_allReps_ratioCtrl_zNorm.bw' 'l3mbtGFP_all_Frags_rpgcNorm_allReps_avg_ratioCtrl.bw' l3mbt_ChIP_all_Frags_rpgcNorm_allReps_avg_ratioCtrl.bw --samplesLabel 'H4K20me1' 'GFP-L(3)mbt' 'L(3)mbt_Richter_et_al' --referencePoint TSS --beforeRegionStartLength 3000 --afterRegionStartLength 3000  --unscaled5prime 0 --unscaled3prime 0 --sortRegions 'keep' --sortUsing 'mean' --averageTypeBins 'mean'  --missingDataAsZero --binSize 50
deeptools plotHeatmap
```
### Expected outputs
Plots in FIGURES 6C and D

## Spike-in normalization for GFP CUT&RUN
### Required files
#### GFP CUT&RUN fastq files (GSE268820)
GSM8299960-GSM8299967
#### GFP_SRPMC_scaling_factors.txt
#### [CUT&RUN pipeline](https://github.com/snystrom/cutNrun-pipeline)
#### spikeNorm_SRPMC.sh ##NEED THIS STILL
### Run code
### Expected outputs
Spike-in normalized bigwigs

## Make H4K20 mutant L(3)mbt heatmaps
### Required files
#### GFP-L3mbt.vs.OregonR.peaks.bed
#### spike-in normalized bigwigs from [Spike-in normalization for GFP CUT&RUN](#Spike-in_normalization_for_GFP_CUT&RUN)
### Run code
```
deeptools computeMatrix reference-point --regionsFileName 'GFP-L3mbt.vs.OregonR.peaks.bed' --scoreFileName 'GFP-L3mbt_spikeNorm_GFP_allReps_avg.bw' 'HWT_His4rnull_spikeNorm_GFP_allReps_avg.bw' 'K20R_His4rnull_spikeNorm_GFP_allReps_avg.bw' --samplesLabel 'l3mbtGFP' 'HWT_His4rnull_l3mbtGFP' --referencePoint center --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --unscaled5prime 0 --unscaled3prime 0 --sortRegions 'keep' --sortUsing 'mean' --averageTypeBins 'mean' --missingDataAsZero --binSize 50
deeptools plotHeatmap

deeptools computeMatrix reference-point --regionsFileName 'sorted/filtered regions from l3mbtGFP and HWT_His4rnull_l3mbtGFP' --scoreFileName 'K20R_His4rnull_spikeNorm_GFP_allReps_avg.bw' --samplesLabel 'H4K20R_His4rnull_l3mbtGFP' --referencePoint center --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --unscaled5prime 0 --unscaled3prime 0 --sortRegions 'keep' --sortUsing 'mean' --averageTypeBins 'mean' --missingDataAsZero --binSize 50
deeptools plotHeatmap
```
### Expected outputs
Plots in FIGURE 6F

