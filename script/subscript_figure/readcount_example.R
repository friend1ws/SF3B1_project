library(tidyverse)

source("../conf/plot_config.R")

D <- read_tsv("../output/recount2/TCGA/TCGA_recount_SF3B1_junction_summary.txt")

splicing_key_list <- D %>% group_by(Splicing_Key) %>% 
  summarize(Max_Read_Count1 = max(Read_Count1), 
            Max_Read_Count2 = max(Read_Count2), 
            Max_Read_Ratio1 = max(Read_Count1 / (Read_Count1 + Read_Count2 + 1)),
            Mean_Read_Ratio1 = mean(Read_Count1 / (Read_Count1 + Read_Count2 + 1))) %>%
  filter(Max_Read_Count1 >= 3, Max_Read_Count2 >= 3, Max_Read_Ratio1 >= 0.05) %>% .$Splicing_Key

target_sample <- D %>% 
  filter(grepl("SF3B1:K700|SF3B1:K666|SF3B1:H662|SF3B1:N626|SF3B1:R625|SF3B1:E622|SF3B1:G740|SF3B1:K741|SF3B1:G742|None", Mutation_Info)) %>%
  .$Sample_Name %>% unique()

print(splicing_key_list[1])
D_target <- D %>% filter(Sample_Name %in% target_sample) %>%
  filter(Splicing_Key %in% splicing_key_list[1]) %>% 
  mutate(status = ifelse(Mutation_Info != "None", "SF3B1 mt", "SF3B1 wt"))
D_target$status <- factor(D_target$status, levels = c("SF3B1 wt", "SF3B1 mt"))

ggplot(D_target, aes(x = Read_Count1 + Read_Count2, y = Read_Count1, colour = as.character(status))) +
  geom_point(size = 0.7) +
  my_theme() +
  facet_grid(.~status) +
  xlim(c(0, 400)) +
  labs(x = "Total read count", y = "Abberant read count") +
  guides(colour = FALSE) 


ggsave("../figure/readcount_example.tiff", width = 12, height = 6, dpi = 600, units = "cm")


