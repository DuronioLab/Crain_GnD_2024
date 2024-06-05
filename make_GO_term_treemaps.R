#Supplemental Figures S2-S5 - Semantic similarity clustering of GO terms for clusters in FIGURE 3E
library(rrvgo)

go_terms_treemap <- function(go_terms.df,ont,method){
  go_filt <- go_terms.df %>% 
    filter(grepl(paste(ont),source))
  
  simMatrix <- calculateSimMatrix(go_filt$term_id,
                                  orgdb="org.Dm.eg.db",
                                  ont="BP",
                                  method="Rel")
  
  scores <- setNames(-log10(go_filt$p_value),
                     go_filt$term_id)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Dm.eg.db")
  
  treemapPlot(reducedTerms,size = "score")
  return(reducedTerms)
}

go_terms_treemap(k20me1_normCounts_cluster1_go,"BP","Rel")
go_terms_treemap(k20me1_normCounts_cluster2_go,"BP","Rel")
go_terms_treemap(k20me1_normCounts_cluster3_go,"BP","Rel")
go_terms_treemap(k20me1_normCounts_cluster4_go,"BP","Rel")
go_terms_treemap(k20me1_normCounts_cluster5_go,"BP","Rel")
go_terms_treemap(k20me1_normCounts_cluster6_go,"BP","Rel")

go_terms_treemap(nok20me1_normCounts_cluster1_go,"BP","Rel")
go_terms_treemap(nok20me1_normCounts_cluster2_go,"BP","Rel")
go_terms_treemap(nok20me1_normCounts_cluster3_go,"BP","Rel")
go_terms_treemap(nok20me1_normCounts_cluster4_go,"BP","Rel")
go_terms_treemap(nok20me1_normCounts_cluster5_go,"BP","Rel")
go_terms_treemap(nok20me1_normCounts_cluster6_go,"BP","Rel")

go_terms_treemap(all_normCounts_cluster1_go,"BP","Rel")
go_terms_treemap(all_normCounts_cluster2_go,"BP","Rel")
go_terms_treemap(all_normCounts_cluster3_go,"BP","Rel")
go_terms_treemap(all_normCounts_cluster4_go,"BP","Rel")
go_terms_treemap(all_normCounts_cluster5_go,"BP","Rel")
go_terms_treemap(all_normCounts_cluster6_go,"BP","Rel")