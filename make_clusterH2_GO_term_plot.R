#FIGURE 3F - GO enrichment analysis of Cluster H2 from FIGURE 3E
library(simplifyEnrichment)
go_terms_htmap <- function(go_terms.df,ont){
  go_filt <- go_terms.df %>% 
    filter(grepl(paste(ont),source)) %>% 
    filter(!grepl("process",term_name)) %>% 
    filter(!grepl("involved",term_name)) %>% 
    filter(!grepl("component",term_name)) %>%
    filter(!grepl("regulation",term_name)) %>%
    filter(p_value < 1e-2)
  
  mat <- GO_similarity(go_filt$term_id,
                       ont = paste(ont)
  )
  df <- simplifyGO(mat)
  colnames(df) <- c("term_id","cluster")
  
  go_cluster <- go_filt %>% 
    left_join(df) %>%
    group_by(cluster) %>%
    filter(n() >= 2) %>%
    ungroup() %>% 
    arrange(cluster,desc(term_size))
  
  htmap <- go_cluster %>% 
    dplyr::select(p_value,term_name) %>% 
    column_to_rownames("term_name")
  htmap[lapply(htmap$parents, length) == 1, , drop = FALSE]
  htmap <- as.matrix(htmap)
  # Get unique values in the "cluster" column
  unique_clusters <- unique(go_cluster$cluster)
  
  # Generate a palette of distinct colors
  distinct_colors <- brewer.pal(n = length(unique_clusters), name = "Paired")
  distinct_colors <- sample(distinct_colors, length(unique_clusters))
  clust_color <- setNames(distinct_colors, unique_clusters)
  
  color_list <- circlize::colorRamp2(c(-1e-40, 1e-30,-1e-20, 1e-10,1e-5), c("#032C3A","#065774", "#0982AE","#8BDCF9","#C5EEFC"))
  
  cluster_anno = rowAnnotation(cluster = as.factor(go_cluster$cluster),
                               col = list(cluster = clust_color))
  
  go_htmap <- Heatmap(htmap,
                      col = color_list,
                      cluster_rows = FALSE,
                      show_row_names = TRUE,
                      row_names_gp = gpar(fontsize = 4),
                      heatmap_legend_param = list(at = c(1e-40, 1e-30,1e-20, 1e-10,1e-5),
                                                  break_dist = 1),
                      right_annotation = cluster_anno,
                      row_split = as.factor(go_cluster$cluster),
                      row_gap = unit(3, "mm"),
                      heatmap_width = unit(3, "in"),
                      heatmap_height = unit(5, "in")
  )
  
  draw(go_htmap)
}

go_terms_htmap(k20me1_normCounts_cluster2_go,"BP")