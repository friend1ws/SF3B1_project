library(ggplot2)
library(RColorBrewer)

source("../conf/plot_config.R")

A <- read.table("../output/recount2/TCGA/TCGA.pred.result.recount2.zibb.txt", sep = "\t", stringsAsFactors = FALSE) %>%
  select(Sample_Name = V1, Cancer_Type = V2, Mutation_Info = V3, Score = V4) %>%
  filter(Mutation_Info != "None")

A$Mutation_Pos <- as.numeric(str_sub(A$Mutation_Info, start = 8, end = -2))


A[!A$Cancer_Type %in% c("UVM", "BRCA", "SKCM", "LIHC", "LUAD", "BLCA"), "Cancer_Type"] <- "Other"
A$Cancer_Type <- factor(A$Cancer_Type,
                        levels = c("BLCA","BRCA", "LIHC", "LUAD", "SKCM", "UVM", "Other"))


ggplot(A %>% filter(Mutation_Pos >= 580 & Mutation_Pos <= 820), 
       aes(x = Mutation_Pos, y = Score, colour = Cancer_Type)) + 
  geom_point(size = rel(0.8), alpha = 0.9) +
  my_theme() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red")+
  scale_x_continuous(breaks=c(600,625,666,700, 740, 781), 
                     labels=c("600", "625", "666", "700", "740", "781")) +
  labs(x = "AA position", y = "SF3B1ness score", colour = "") +
  scale_colour_manual(values = c(brewer.pal(6, "Set2"), "grey60")) +
  theme(legend.position = "bottom")

                      
ggsave("../figure/SF3B1ness_score_aa_position.tiff", width = 20, height = 6, dpi = 600, units = "cm")


                      
  






