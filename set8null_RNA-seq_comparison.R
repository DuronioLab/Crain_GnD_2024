#Supplemental FIGURE S1A - Set8null (Crain et al. 2024) vs. Set8null (Bamgbose et al. 2024)
#import sample sheet
samples_set8 <- read.tsv("sample_sheet_RNA-seq_Set8null_comp.txt")

txi_set8 <- tximport(files_set8, type="salmon", tx2gene=tx2gene)

ddsTxi_set8 <- DESeqDataSetFromTximport(txi_set8,
                                        colData = samples_set8,
                                        design = ~ sample)

summary(ddsTxi_set8)
ddsTxi_set8.df <- as.data.frame(assay(ddsTxi_set8))

ddsTxi_set8.test <- DESeq(ddsTxi_set8)
ddsTxi.test_set8.df <- as.data.frame(assay(ddsTxi_set8.test))

#plot mean normalized counts Crain vs. Bamgbose
set8_comp <- ggscatter(vsd_set8_means.df,
                       x = "Set8null_mean",
                       y = "set8null_polyA_mean",
                       add = "reg.line",
                       cor.coef = TRUE)
set8_comp