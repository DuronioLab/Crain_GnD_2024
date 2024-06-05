#FIGURE 4E - Quantification of eye clones 
##Figure 4
#eye clones boxplot
eye_clone_data <- read_tsv("HWT_K20A_K20R_mosaic_eye_quant.txt")
eye_clone_data_long <- eye_clone_data %>% 
  pivot_longer(cols = !Sex, names_to = "genotype") %>% 
  mutate(genotype_sex = paste0(genotype, "_", Sex)) %>% 
  filter(!grepl("NA",genotype_sex))

eye_clone_stat.test <- compare_means(value ~ genotype,
                                     data = eye_clone_data_long,
                                     method = "t.test")

eye_clone_graph <- ggboxplot(eye_clone_data_long,
                             x = "genotype",
                             y = "value",
                             add = "jitter",
                             fill = "genotype",
                             palette = "jco")
eye_clone_graph + 
  rremove("legend") +
  ylim(0,80)