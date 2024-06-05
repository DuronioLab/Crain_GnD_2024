##FIGURE 1B - H4K20me1 overlap with genomic features
library(ChIPseeker)
h4k20me1_peak <- readPeakFile("H4K20me1.vs.no_primary.peaks.bed")
h4k20me1_peak.df <- as.data.frame(h4k20me1_peak)

h4k20me1_peakAnno <- annotatePeak(h4k20me1_peak, tssRegion=c(-300, 300),
                                  TxDb=txdb, annoDb="org.Dm.eg.db"
)

h4k20me1_peakAnno.df <- as.data.frame(h4k20me1_peakAnno)

upsetplot(peakAnno,vennpie=TRUE)