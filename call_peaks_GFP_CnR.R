##Load required packages
library(csaw)
library(edgeR)
library(tidyverse)

#Call H4K20me1 peaks using csaw and edgeR
##Based on vignette from
#http://bioconductor.org/books/3.18/csawBook/index.html#introduction
#Source: https://github.com/LTLA/csawUsersGuide

set.seed(10)

##GFP-L(3)mbt
#read in annotation file with sample names and paths to bam files - Oregon-R and Oregon-R no_primary ctrl
gfp_data <- read_tsv("GFP-l3mbt_CnR_sample_sheet.txt") %>% 
  filter(!grepl("HWT", Name)) %>% 
  filter(!grepl("K20R", Name))

#load csaw package
library(csaw)

#set parameters to paired-end mode and only consider mapped reads with quality score > 10
param <- readParam(minq=10, pe = "both")
param

#bin counts into 150bp windows with 50bp overlap
win.data <- windowCounts(gfp_data$Path_to_bam, param=param, width=150, spacing = 50)

#bin into larger windows for background thresholding
bins <- windowCounts(gfp_data$Path_to_bam, bin=TRUE, width=10000, param=param)

#compute filter statistics
filter.stat <- filterWindowsGlobal(win.data,bins)

#look at binned data to determine background threshold
hist(filter.stat$filter, main="", breaks=50,
     xlab="Log-fold change from control")
abline(v=log2(2), col="red", lwd=2)

#identify windows that are > log2(2) above background
keep <- filter.stat$filter > log2(2)
summary(keep)

#filter to retain only windows above background
filtered.data <- win.data[keep,]

#calculate compositional bias factors
comp_norm <- normFactors(bins, se.out=win.data)
normfacs <- comp_norm$norm.factors
normfacs

##Convert RangedSummarizedExperiment object into a DGEList for modelling with edgeR
#load edgeR
library(edgeR)
y <- asDGEList(filtered.data,
               norm.factors=normfacs)

##Construct a design matrix for our experimental design
genotype <- gfp_data[1:4,]$Name
genotype[grepl("GFP", genotype)] <- "l3mbtGFP"
genotype[grepl("OregonR", genotype)] <- "OregonR"

genotype <- factor(genotype)
design <- model.matrix(~0+genotype)
colnames(design)[1:2] <- levels(genotype)
design

##Estimate the negative binomial (NB) and quasi-likelihood (QL) dispersions for each window
y <- estimateDisp(y, design)
summary(y$trended.dispersion)

plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$df.prior)

plotQLDisp(fit)

#make PCA plot to check that biological replicates of each genotype cluster together
colorPal <- c("blue","red")
logCPM <- cpm(y, log=TRUE, prior.count=2)

plotMDS(logCPM, top=10000, labels=genotype,
        col=colorPal[as.integer(genotype)])

#Test for DB between conditions in each window using the QL F-test.
contrast <- makeContrasts(l3mbtGFP-OregonR, levels=design)
res <- glmQLFTest(fit, contrast=contrast)
res.df <- as.data.frame(res)

unmerged <- res$table
unmerged.sig.pos <- unmerged %>% 
  filter(logFC > 3 & PValue < 0.05)

#select windows that are significantly enriched (logFC > 3) in l3mbtGFP compared to Oregon-R
unmerged.sig.pos.keep <- unmerged$logFC > 3 & unmerged$PValue < 0.05
summary(unmerged.sig.pos.keep)

#filter for windows with significantly enriched signal in Oregon-R compared to Oregon-R no_primary
filtered.data_unmerged_up <- filtered.data[unmerged.sig.pos.keep,]

#merge significant windows within 250bp
merged <- mergeResultsList(list(filtered.data_unmerged_up), 
                           tab.list=list(unmerged.sig.pos),
                           equiweight=TRUE, tol=250)
merged$regions

out.ranges.df <- data.frame(merged$regions,merged$combined)
out.ranges_sig.df <- out.ranges.df %>% 
  filter(FDR < 0.05 & direction == "up")

gfp_cnr_peaks.bed <- out.ranges_sig.df[1:4]
colnames(gfp_cnr_peaks.bed) <- c("seqnames","start","end","width")
write_tsv(gfp_cnr_peaks.bed,file = "GFP-L3mbt.vs.OregonR.peaks.bed",col_names = FALSE)
