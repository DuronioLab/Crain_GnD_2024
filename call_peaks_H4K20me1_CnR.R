###FIGURE 1: H4K20me1 peak calling and gene expression correlation

#Call H4K20me1 peaks using csaw and edgeR
##Based on vignette from
#http://bioconductor.org/books/3.18/csawBook/index.html#introduction
#Source: https://github.com/LTLA/csawUsersGuide

##H4K20me1
#read in annotation file with sample names and paths to bam files - Oregon-R and Oregon-R no_primary ctrl
k20me1_data <- read_tsv("H4K20me1_CnR_sample_sheet.txt")

#load csaw package
library(csaw)

#set parameters to paired-end mode and only consider mapped reads with quality score > 10
param <- readParam(minq=10, pe = "both")
param

#bin counts into 150bp windows with 50bp overlap
win.data <- windowCounts(k20me1_data$Path, param=param, width=150, spacing = 50)

#bin into larger windows for background thresholding
bins <- windowCounts(k20me1_data$Path, bin=TRUE, width=10000, param=param)

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
genotype <- k20me1_data[1:6,]$Name
genotype[grep("or", genotype)] <- "or"
genotype[grep("ornopr", genotype)] <- "ornopr"

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
contrast <- makeContrasts(or-ornopr, levels=design)
res <- glmQLFTest(fit, contrast=contrast)
res.df <- as.data.frame(res)

#select windows that are significantly enriched (logFC > 2) in Oregon-R compared to Oregon-R no_primary
unmerged.sig.pos.keep <- unmerged$logFC > 2 & unmerged$PValue < 0.05
summary(unmerged.sig.pos.keep)

#filter for windows with significantly enriched signal in Oregon-R compared to Oregon-R no_primary
filtered.data_unmerged_up <- filtered.data[unmerged.sig.pos.keep,]

#merge significant windows within 1kb
merged <- mergeResultsList(list(filtered.data_unmerged_up), 
                           tab.list=list(res$table),
                           equiweight=TRUE, tol=1000)
merged$regions

out.ranges.df <- data.frame(merged$regions,merged$combined)
out.ranges_sig.df <- out.ranges.df %>% 
  filter(FDR < 0.05 & direction == "up")

k20me1_cnr_peaks.bed <- out.ranges_sig.df[1:4]
colnames(k20me1_cnr_peaks.bed) <- c("seqnames","start","end","width")
write_tsv(k20me1_cnr_peaks.bed,file = "H4K20me1.vs.no_primary.peaks.bed",col_names = FALSE)