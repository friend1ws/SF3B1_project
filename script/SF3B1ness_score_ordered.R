library(ggplot2)


source("../conf/plot_config.R")

A <- read.table("../output/recount2/TCGA/TCGA.pred.result.recount2.zibb.txt", sep = "\t", stringsAsFactors = FALSE) %>%
  select(Sample_Name = V1, Cancer_Type = V2, Mutation_Info = V3, Score = V4)


# A[grepl("SF3B1", A$Mutation_Info) & !grepl("K700|K666|H662|N626|R625|E622|G740|K741|G742", A$Mutation_Info), "Mutation_Info"] <- "SF3B1:Other"

# A[grepl("SF3B1:K700", A$Mutation_Info), "Mutation_Info"] <- "SF3B1:K700"
# A[grepl("SF3B1:K666", A$Mutation_Info), "Mutation_Info"] <- "SF3B1:K666"
# A[grepl("SF3B1:H662", A$Mutation_Info), "Mutation_Info"] <- "SF3B1:H662"
# A[grepl("SF3B1:E622", A$Mutation_Info), "Mutation_Info"] <- "SF3B1:E622"
# A[grepl("SF3B1:R625", A$Mutation_Info), "Mutation_Info"] <- "SF3B1:R625"
# A[grepl("SF3B1:N626", A$Mutation_Info), "Mutation_Info"] <- "SF3B1:N626"
# A[grepl("SF3B1:K741", A$Mutation_Info), "Mutation_Info"] <- "SF3B1:K741"
# A[grepl("SF3B1:G740", A$Mutation_Info), "Mutation_Info"] <- "SF3B1:G740"
# A[grepl("SF3B1:G742", A$Mutation_Info), "Mutation_Info"] <- "SF3B1:G742"

A[grepl("K700|K666|H662|N626|R625|E622|G740|K741|G742", A$Mutation_Info), "Mutation_Info"] <- "SF3B1:Hotspot"

A[grepl("SF3B1", A$Mutation_Info) & !grepl("SF3B1:Hotspot", A$Mutation_Info), "Mutation_Info"] <- "SF3B1:Other"

A$Mutation_Info <- factor(A$Mutation_Info, levels = c("SF3B1:Hotspot", "SF3B1:Other", "None"))

B <- A %>% filter(Score > -200)



ggplot(B, aes(x = reorder(Sample_Name, desc(Score)), y = Score, fill = Mutation_Info)) + 
  geom_bar(stat = "identity", width = 0.6) +
  labs(x = "Sample", y = "SF3B1ness score", fill = "") +
  my_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = c("#F8766D", "#00BA38", "grey70"))

ggsave("../figure/SF3B1ness_score_ordered.tiff", width = 20, height = 6, dpi = 600, units = "cm")

