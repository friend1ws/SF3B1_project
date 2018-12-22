library(ggplot2)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
junc_file <- args[1]
param_file <- args[2]
out_dir <- args[3]

P <- read.table(param_file, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(V4 <= 0.011 & V7 >= 0.15)

target_junc_list <- P$V1


D <- read.table(junc_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE) %>%
  filter(Splicing_Key %in% target_junc_list)

total_count <- D$Read_Count1 + D$Read_Count2
D[total_count >= 200, "Read_Count1"] <- 200 * D$Read_Count1[total_count >= 200] / total_count[total_count >= 200]
D[total_count >= 200, "Read_Count2"] <- 200 * D$Read_Count2[total_count >= 200] / total_count[total_count >= 200]

sample_list <- unique(D$Sample_Name)

filtD <- left_join(expand.grid(Splicing_Key = target_junc_list, Sample_Name = sample_list, stringsAsFactors = FALSE), D, by = c("Splicing_Key", "Sample_Name"))

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

for (sample in sample_list) {
  print(sample)
  tD <- filtD %>% filter(Sample_Name == sample)
  ggplot(tD, aes(x = Splicing_Key2, y = Count, fill = Splicing_Type)) + 
    geom_bar(stat = "identity") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = "", y = "Count", fill = "")
  

  ggsave(paste(out_dir, "/", sample, ".tiff", sep = ""), width = 16, height = 8, dpi = 300)
  
}
