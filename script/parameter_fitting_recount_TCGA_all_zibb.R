library(tidyr)
library(dplyr)
library(ggplot2)

source("subscript_zibb/zibb_functions.R")

D <- read.table("../output/recount2/TCGA/TCGA_recount_SF3B1_junction_summary.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

tD <- D %>% group_by(Splicing_Key) %>% 
  summarize(Max_Read_Count1 = max(Read_Count1), 
            Max_Read_Count2 = max(Read_Count2), 
            Max_Read_Ratio1 = max(Read_Count1 / (Read_Count1 + Read_Count2 + 1)),
            Mean_Read_Ratio1 = mean(Read_Count1 / (Read_Count1 + Read_Count2 + 1)))

tD2 <- tD %>% filter(Max_Read_Count1 >= 3, Max_Read_Count2 >= 3, Max_Read_Ratio1 >= 0.05)
splicing_key_list <- tD2$Splicing_Key
all_sample_name <- unique(D$Sample_Name)

Pred <- data.frame(sample = c(), cancer_type = c(), mutation_info = c(), score = c())
Params <- c()
  
target_D <- D %>% filter(grepl("SF3B1:K700|SF3B1:K666|SF3B1:H662|SF3B1:N626|SF3B1:R625|SF3B1:E622|SF3B1:G740|SF3B1:K741|SF3B1:G742|None", Mutation_Info))

# test_sample_name <- all_sample_name[seq(k, length(all_sample_name), 2)]
# training_sample_name <- setdiff(all_sample_name, test_sample_name) 
training_sample_name <- all_sample_name

# removing those having SF3B1 VUS mutations
training_sample_name <- intersect(training_sample_name, unique(target_D$Sample_Name))
  
  
params <- c()
junc <- c()
for (i in 1:length(splicing_key_list)) {
  
  tD <- D %>% filter(Splicing_Key == splicing_key_list[i], Sample_Name %in% training_sample_name)
  tD$status <- rep(0, nrow(tD))
  tD$status[tD$Mutation_Info != "None"] <- 1
  
  ttD <- tD %>% group_by(status) %>% summarize(Ratio = mean(Read_Count1 / (Read_Count1 + Read_Count2 + 1)))
  if (sum(tD$status) < 3) next
  if (ttD$Ratio[ttD$status == 1] < ttD$Ratio[ttD$status == 0] * 3) next
  if (quantile(tD$Read_Count2)[4] > 1000) next

  print(i)
  
  class_ind <- tD$status
  vector_n <- tD$Read_Count1[class_ind == 0] + tD$Read_Count2[class_ind == 0]
  vector_k <- tD$Read_Count1[class_ind == 0]
  cret0 <- zibb_optim(vector_n, vector_k)
  
  vector_n <- tD$Read_Count1[class_ind == 1] + tD$Read_Count2[class_ind == 1]
  vector_k <- tD$Read_Count1[class_ind == 1]
  cret1 <- zibb_optim(vector_n, vector_k)
  
  params <- rbind(params, c(cret0$par, cret1$par))
  junc <- c(junc, splicing_key_list[i])

}

rownames(params) <- junc


write.table(params, "../output/recount2/TCGA/param_matrix.recount2.zibb.txt", sep = "\t", col.names = FALSE, quote = FALSE)




