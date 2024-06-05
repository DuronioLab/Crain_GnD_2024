##FIGURE 1E
#Correlation of genes with H4K20me1 and gene expression
##For Salmon transcript quantification - see Github
#transfer Salmon output (entire salmon_quant directory) to current working directory

#import sample sheet
samples_wt <- read.tsv("Salmon_wt_sample_sheet.txt")

#set file names to match Salmon output
files_wt <- file.path(samples_wt$useName, "quant.sf")

#check if file labels are correct and exist in current working directory
all(file.exists(files_wt))
which(!file.exists(files_wt))

#import gtf file (downloaded from FlyBase)
gtf <- import("dmel-all-r6.55.gtf") %>% 
  as.data.frame()

#transcript name to gene name conversion
tx2gene <- gtf %>% 
  dplyr::select(transcript_id,gene_symbol) %>% 
  na.omit()

#import Salmon quant data for differential expression analysis with DESeq2
txi_wt <- tximport(files_wt, type="salmon", tx2gene=tx2gene)

#make design matrix for DESeq2
ddsTxi_wt <- DESeqDataSetFromTximport(txi_wt,
                                      colData = samples_wt,
                                      design = ~ sample)

#Run DESeq2
ddsTxi_wt.test <- DESeq(ddsTxi_wt)

#Normalize count data with variance stabilizing transformation
vsd_wt <- vst(ddsTxi_wt.test,blind = FALSE)

vsd_wt.df <- as.data.frame(assay(vsd_wt))
colnames(vsd_wt.df) <- samples_wt$useName
vsd_wt.df <- vsd_wt.df %>% 
  rownames_to_column("gene_symbol")

#calculate mean and sd abundance for each gene across reps  
vsd_wt_avg <- vsd_wt.df %>% 
  rowwise() %>% 
  mutate(or_3wl_avg = mean(c_across(contains("OregonR")), na.rm = TRUE),
         or_3wl_sd = sd(c_across(contains("OregonR")), na.rm = TRUE)) %>% 
  mutate(yw_wing_avg = mean(c_across(contains("yw")), na.rm = TRUE),
         yw_wing_sd = sd(c_across(contains("yw")), na.rm = TRUE))

#label genes in RNA abundance data with amount of h4k20me1 overlap
vsd_or_3wl_avg_k20me1_0.75 <- vsd_wt_avg %>% 
  filter(gene_symbol %in% k20me1_0.75_genes$gene_symbol) %>% 
  mutate(k20me1_status = "whole_k20me1_0.75") %>% 
  select(gene_symbol,or_3wl_avg,or_3wl_sd,k20me1_status)
colnames(vsd_or_3wl_avg_k20me1_0.75) <- c("gene_symbol","RNA_avg","RNA_sd","k20me1_status")

vsd_or_3wl_avg_k20me1_0.5 <- vsd_wt_avg %>% 
  filter(gene_symbol %in% k20me1_0.5_genes$gene_symbol) %>% 
  mutate(k20me1_status = "whole_k20me1_0.5") %>% 
  select(gene_symbol,or_3wl_avg,or_3wl_sd,k20me1_status)
colnames(vsd_or_3wl_avg_k20me1_0.5) <- c("gene_symbol","RNA_avg","RNA_sd","k20me1_status")

vsd_or_3wl_avg_k20me1_0.25 <- vsd_wt_avg %>% 
  filter(gene_symbol %in% k20me1_0.25_genes$gene_symbol) %>% 
  mutate(k20me1_status = "whole_k20me1_0.25") %>% 
  select(gene_symbol,or_3wl_avg,or_3wl_sd,k20me1_status)
colnames(vsd_or_3wl_avg_k20me1_0.25) <- c("gene_symbol","RNA_avg","RNA_sd","k20me1_status")

vsd_or_3wl_avg_k20me1_0.1 <- vsd_wt_avg %>% 
  filter(gene_symbol %in% k20me1_0.1_genes$gene_symbol) %>% 
  mutate(k20me1_status = "whole_k20me1_0.1") %>% 
  select(gene_symbol,or_3wl_avg,or_3wl_sd,k20me1_status)
colnames(vsd_or_3wl_avg_k20me1_0.1) <- c("gene_symbol","RNA_avg","RNA_sd","k20me1_status")

vsd_or_3wl_avg_nok20me1 <- vsd_wt_avg %>% 
  filter(!gene_symbol %in% k20me1_0.1_genes$gene_symbol) %>% 
  mutate(k20me1_status = "whole_nok20me1") %>% 
  select(gene_symbol,or_3wl_avg,or_3wl_sd,k20me1_status)
colnames(vsd_or_3wl_avg_nok20me1) <- c("gene_symbol","RNA_avg","RNA_sd","k20me1_status")

###wing
vsd_yw_wing_avg_k20me1_0.75 <- vsd_wt_avg %>% 
  filter(gene_symbol %in% k20me1_0.75_genes$gene_symbol) %>% 
  mutate(k20me1_status = "wing_k20me1_0.75") %>%
  select(gene_symbol,yw_wing_avg,yw_wing_sd,k20me1_status)
colnames(vsd_yw_wing_avg_k20me1_0.75) <- c("gene_symbol","RNA_avg","RNA_sd","k20me1_status")

vsd_yw_wing_avg_k20me1_0.5 <- vsd_wt_avg %>% 
  filter(gene_symbol %in% k20me1_0.5_genes$gene_symbol) %>% 
  mutate(k20me1_status = "wing_k20me1_0.5") %>%
  select(gene_symbol,yw_wing_avg,yw_wing_sd,k20me1_status)
colnames(vsd_yw_wing_avg_k20me1_0.5) <- c("gene_symbol","RNA_avg","RNA_sd","k20me1_status")

vsd_yw_wing_avg_k20me1_0.25 <- vsd_wt_avg %>% 
  filter(gene_symbol %in% k20me1_0.25_genes$gene_symbol) %>% 
  mutate(k20me1_status = "wing_k20me1_0.25") %>% 
  select(gene_symbol,yw_wing_avg,yw_wing_sd,k20me1_status)
colnames(vsd_yw_wing_avg_k20me1_0.25) <- c("gene_symbol","RNA_avg","RNA_sd","k20me1_status")

vsd_yw_wing_avg_k20me1_0.1 <- vsd_wt_avg %>% 
  filter(gene_symbol %in% k20me1_0.1_genes$gene_symbol) %>% 
  mutate(k20me1_status = "wing_k20me1_0.1") %>% 
  select(gene_symbol,yw_wing_avg,yw_wing_sd,k20me1_status)
colnames(vsd_yw_wing_avg_k20me1_0.1) <- c("gene_symbol","RNA_avg","RNA_sd","k20me1_status")

vsd_yw_wing_avg_nok20me1 <- vsd_wt_avg %>% 
  filter(!gene_symbol %in% k20me1_0.1_genes$gene_symbol) %>% 
  mutate(k20me1_status = "wing_nok20me1") %>% 
  select(gene_symbol,yw_wing_avg,yw_wing_sd,k20me1_status)
colnames(vsd_yw_wing_avg_nok20me1) <- c("gene_symbol","RNA_avg","RNA_sd","k20me1_status")

#generate random set of genes that is the same size as genes with >50% H4K20me1 overlap and <50% overlap
vsd_wt_avg_select <- vsd_wt_avg %>% 
  select(gene_symbol,or_3wl_avg,or_3wl_sd)
colnames(vsd_wt_avg_select) <- c("gene_symbol","RNA_avg","RNA_sd")

set.seed(123)  # for reproducibility

random_genes_vsd_yes <- sample(vsd_wt_avg_select$gene_symbol, nrow(k20me1_genes), replace = FALSE)
remaining_genes <- setdiff(vsd_wt_avg_select$gene_symbol, random_genes_vsd_yes)
random_genes_vsd_no <- sample(remaining_genes, nrow(vsd_or_3wl_avg_nok20me1), replace = FALSE)
vsd_wt_avg_random <- vsd_wt_avg_select %>% 
  mutate(k20me1_status = case_when(gene_symbol %in% random_genes_vsd_yes ~ "yes",
                                   gene_symbol %in% random_genes_vsd_no ~ "no",
                                   TRUE ~ NA)) 

table(vsd_wt_avg_random$k20me1_status)

k20me1_RNA_boxplot_vsd.df <- rbind(vsd_or_3wl_avg_k20me1_0.75,
                                   vsd_or_3wl_avg_k20me1_0.5,
                                   vsd_or_3wl_avg_k20me1_0.25,
                                   vsd_or_3wl_avg_k20me1_0.1,
                                   vsd_or_3wl_avg_nok20me1,
                                   vsd_yw_wing_avg_k20me1_0.75,
                                   vsd_yw_wing_avg_k20me1_0.5,
                                   vsd_yw_wing_avg_k20me1_0.25,
                                   vsd_yw_wing_avg_k20me1_0.1,
                                   vsd_yw_wing_avg_nok20me1,
                                   vsd_wt_avg_random
) %>% 
  na.omit()
table(k20me1_RNA_boxplot_vsd.df$k20me1_status)

test <- k20me1_genes %>% 
  filter(!gene_symbol %in% vsd_or_3wl_avg_k20me1_0.5$gene_symbol)

test2 <- k20me1_0.5_genes %>% 
  filter(gene_symbol %in% resLFC_Set8null_ashr.df2$gene_symbol)

#FIGURE 1E - boxplot of abundance
#wilcox.test compares medians even though it is called within compare_means()
k20me1_wt_stat.test <- compare_means(RNA_avg ~ k20me1_status,
                                     data = k20me1_RNA_boxplot_vsd.df,
                                     method = "wilcox.test",
                                     p.adjust.method = "BH")
k20me1_RNA_stat.test_vsd.means <- k20me1_RNA_boxplot_vsd.df %>% 
  group_by(k20me1_status) %>% 
  summarize(mean = mean(RNA_avg), median = median(RNA_avg))

k20me1_RNA_stat.test_vsd.means <- k20me1_RNA_boxplot_vsd.df %>% 
  group_by(k20me1_status) %>% 
  summarize(mean = mean(RNA_avg), median = median(RNA_avg), wilcox = wilcox.test())

boxplot_colors <- c("#CCFAE5","#00E67C","#009953","#006838","#002C17",
                    "#CCFAE5","#00E67C","#009953","#006838","#002C17",
                    "#808285","#BCBEC0")
#Make Figure 1E
k20me1_RNA_boxplot_vsd <- ggboxplot(k20me1_RNA_boxplot_vsd.df,
                                    x = "k20me1_status",
                                    y = "RNA_avg",
                                    order = c("whole_nok20me1",
                                              "whole_k20me1_0.1",
                                              "whole_k20me1_0.25",
                                              "whole_k20me1_0.5",
                                              "whole_k20me1_0.75",
                                              "wing_nok20me1",
                                              "wing_k20me1_0.1",
                                              "wing_k20me1_0.25",
                                              "wing_k20me1_0.5",
                                              "wing_k20me1_0.75",
                                              "yes",
                                              "no"),
                                    fill = "k20me1_status",
                                    palette = boxplot_colors,
                                    alpha = 0.5,
                                    ylim = c(0,20))

k20me1_RNA_boxplot_vsd +
  geom_hline(yintercept=median(vsd_or_3wl_avg_k20me1$RNA_avg), 
             linetype="dashed", color = "red",linewidth=1) +
  rremove("legend")