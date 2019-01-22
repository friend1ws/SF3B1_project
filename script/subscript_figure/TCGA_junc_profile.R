library(tidyverse)
source("../conf/plot_config.R")

param_file <- "../output/recount2/TCGA/param_matrix.recount2.zibb.txt"
junc_file <- "../output/recount2/TCGA/TCGA_recount_SF3B1_junction_summary.txt"

target_junc_list <- read_tsv(param_file, col_names = FALSE) %>% filter(X4 <= 0.011 & X7 >= 0.20) %>% .$X1



D <- read_tsv(junc_file) %>% filter(Splicing_Key %in% target_junc_list)

total_count <- D$Read_Count1 + D$Read_Count2
D[total_count >= 200, "Read_Count1"] <- 200 * D$Read_Count1[total_count >= 200] / total_count[total_count >= 200]
D[total_count >= 200, "Read_Count2"] <- 200 * D$Read_Count2[total_count >= 200] / total_count[total_count >= 200]

# sample_list <- unique(D$Sample_Name)

filtD <- left_join(expand.grid(Splicing_Key = target_junc_list, Sample_Name = sample_list, stringsAsFactors = FALSE), 
                   D, 
                   by = c("Splicing_Key", "Sample_Name"))

filtD$Read_Count2[is.na(filtD$Read_Count2)] <- 0
filtD$Read_Count1[is.na(filtD$Read_Count1)] <- 0

filtD <- filtD %>% 
  mutate(Splicing_Key2 = unlist(lapply(strsplit(Splicing_Key, split = ";"), function(x) {paste(x[1], x[2], x[3], x[4], sep = ",")}))) %>%
  select(Splicing_Key2, Sample_Name, Read_Count1, Read_Count2) %>% 
  gather(key = Splicing_Type, value = Count, Read_Count1, Read_Count2)

filtD$Splicing_Type = factor(filtD$Splicing_Type,
                             levels = c("Read_Count1", "Read_Count2"),
                             labels = c("Abnormal", "Normal"))

filtD$Splicing_Key2 <- factor(filtD$Splicing_Key2)


sample_list <- c("TCGA-VD-A8KH", "TCGA-D8-A27W", "TCGA-55-7576")

filtD2 <- filtD %>% filter(Sample_Name %in% sample_list)
filtD2$Sample_Name <- factor(filtD2$Sample_Name,
                             levels = c("TCGA-VD-A8KH", "TCGA-D8-A27W", "TCGA-55-7576"),
                             labels= c("TCGA-VD-A8KH (UVM, SF3B1:R625)", "TCGA-D8-A27W (BRCA, SF3B1:WT)", "TCGA-55-7576 (LUAD, SF3B1:WT)"))

                           
ggplot(filtD2, aes(x = Splicing_Key2, y = Count, fill = Splicing_Type)) + 
  geom_bar(stat = "identity") + 
  facet_grid(.~Sample_Name) +
  my_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "SF3B1 specific alternative 3'SS loci", y = "Count", fill = "") +
  scale_fill_manual(values = c("#f0027f", "#386cb0"))
  
ggsave("../output/recount2/TCGA/TCGA_junc_profile.tiff", width = 18, height = 6, units = "cm", dpi = 600)
  


