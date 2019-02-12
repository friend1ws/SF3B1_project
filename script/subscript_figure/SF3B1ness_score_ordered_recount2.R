library(tidyverse)
source("../conf/plot_config.R")

D <- read_tsv("../output/sra/hotspot/SF3B1.hotspot.result.info2.txt") %>% filter(Bases >= 500000000) 

scores <- c()
run_id <- c()
mut_type <- c()
for(i in 1:nrow(D)) {
  
  if (D$Mutation_Info[i] == "---") {
    tmut_type <- "WT"
  } else {

    mut_infos <- unlist(strsplit(D$Mutation_Info[i], ';'))
    cosm_count <- 0
    for (j in 1:length(mut_infos)) {
      t_cosm_count <- as.numeric(unlist(strsplit(mut_infos[j], ','))[2])
      if (t_cosm_count > cosm_count) {
        cosm_count <- t_cosm_count
      }
    }
    
    if (cosm_count >= 5) {
      tmut_type <- "SF3B1: #COSMIC>=5"
    } else {
      tmut_type <- "SF3B1: #COSMIC<5"
    }
    
  }
  
  run_id <- c(run_id, D$Run_Accession[i])
  scores <- c(scores, as.numeric(D$SF3B1ness_Score[i]))
  mut_type <- c(mut_type, tmut_type)
  
}

mut_type <- factor(mut_type, levels = c("SF3B1: #COSMIC>=5", "SF3B1: #COSMIC<5", "WT"))

D2 <- data.frame(Run_ID = run_id, Score = scores, Mut_Type = mut_type)


ggplot(D2, aes(x = reorder(Run_ID, desc(Score)), y = Score, fill = Mut_Type)) + 
  geom_bar(stat = "identity", width = 0.6) +
  labs(x = "Sample", y = "SF3B1ness score", fill = "") +
  my_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = c("#7570b3", "#66a61e", "grey70"))

ggsave("../figure/SF3B1ness_score_ordered_recount2.tiff", width = 16, height = 6, dpi = 600, units = "cm")


