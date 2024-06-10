#load required packages
library(tidyverse)
library(DESeq2)
library(tximport)
library(ggpubr)
library(rtracklayer)

#Supplemental FIGURE S1A - Set8null (Crain et al. 2024) vs. Set8null (Bamgbose et al. 2024)
#import sample sheet
samples_set8 <- read_tsv("sample_sheet_RNA-seq_Set8null_comp.txt")

#set file names to match Salmon output
files <- file.path(samples_set8$useName, "quant.sf")

#check if file labels are correct and exist in current working directory
all(file.exists(files))
which(!file.exists(files))

#import gtf file (downloaded from FlyBase)
gtf <- import("dmel-all-r6.55.gtf") %>% 
  as.data.frame()

#transcript name to gene name conversion
tx2gene <- gtf %>% 
  dplyr::select(transcript_id,gene_symbol) %>% 
  na.omit()

txi_set8 <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi_set8,
                                   colData = samples_set8,
                                   design = ~ sample)


ddsTxi.test <- DESeq(ddsTxi)

#Normalize data and make PCA plot
#filter low counts
filt <- rowSums(assay(ddsTxi.test)) > 1

#Normalize count data using variance stabilizing transformation
vsd <- vst(ddsTxi.test[filt,],blind = FALSE)
vsd.df <- as.data.frame(assay(vsd))
colnames(vsd.df) <- samples_set8$useName

#plot PCA
plotPCA(vsd,intgroup="sample")

vsd_set8_means.df <- vsd.df %>% 
  rowwise() %>% 
  mutate(set8null_mean = mean(c_across(Set8null_rep1:Set8null_rep4))) %>% 
  mutate(set8null_polyA_mean = mean(c_across(Set8null_polyA_rep1:Set8null_polyA_rep3)))

#plot mean normalized counts Crain vs. Bamgbose
set8_comp <- ggscatter(vsd_set8_means.df,
                       x = "set8null_mean",
                       y = "set8null_polyA_mean",
                       add = "reg.line",
                       cor.coef = TRUE)
set8_comp
