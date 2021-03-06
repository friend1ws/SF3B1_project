library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
model <- args[1]

source("subscript_model/model_functions.R")

# D <- read.table("../output/recount2/TCGA/TCGA_recount_SF3B1_junction_summary.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
D <- read_tsv("../output/recount2/TCGA/TCGA_recount_SF3B1_junction_summary.txt")

splicing_key_list <- D %>% group_by(Splicing_Key) %>% 
  summarize(Max_Read_Count1 = max(Read_Count1), 
            Max_Read_Count2 = max(Read_Count2), 
            Max_Read_Ratio1 = max(Read_Count1 / (Read_Count1 + Read_Count2 + 1)),
            Mean_Read_Ratio1 = mean(Read_Count1 / (Read_Count1 + Read_Count2 + 1))) %>%
  filter(Max_Read_Count1 >= 3, Max_Read_Count2 >= 3, Max_Read_Ratio1 >= 0.05) %>% .$Splicing_Key


# tD2 <- tD %>% filter(Max_Read_Count1 >= 3, Max_Read_Count2 >= 3, Max_Read_Ratio1 >= 0.05)
# splicing_key_list <- tD2$Splicing_Key
all_sample_name <- unique(D$Sample_Name)

target_D <- D %>% filter(grepl("SF3B1:K700|SF3B1:K666|SF3B1:H662|SF3B1:N626|SF3B1:R625|SF3B1:E622|SF3B1:G740|SF3B1:K741|SF3B1:G742|None", Mutation_Info))
training_sample_name <- all_sample_name

# removing those having SF3B1 VUS mutations
training_sample_name <- intersect(training_sample_name, unique(target_D$Sample_Name))
  
D_train <- D %>% filter(Sample_Name %in% training_sample_name) %>% 
	mutate(status = ifelse(Mutation_Info != "None", 1, 0))

params <- c()
junc <- c()

for (i in 1:length(splicing_key_list)) {
  
  D_train_target <- D_train %>% filter(Splicing_Key == splicing_key_list[i])
  if (sum(D_train_target$status) < 3) next
  if (quantile(D_train_target$Read_Count2)[4] > 1000) next

  D_train_target_sum <- D_train_target %>% group_by(status) %>% summarize(Ratio = mean(Read_Count1 / (Read_Count1 + Read_Count2 + 1)))
  if (D_train_target_sum$Ratio[D_train_target_sum$status == 1] < D_train_target_sum$Ratio[D_train_target_sum$status == 0] * 3) next

  print(i)

  class_ind <- D_train_target$status
  vector_n_0 <- D_train_target$Read_Count1[class_ind == 0] + D_train_target$Read_Count2[class_ind == 0]
  vector_k_0 <- D_train_target$Read_Count1[class_ind == 0]
  
  vector_n_1 <- D_train_target$Read_Count1[class_ind == 1] + D_train_target$Read_Count2[class_ind == 1]
  vector_k_1 <- D_train_target$Read_Count1[class_ind == 1]

  tpar <- switch(model,
                 "p" = c(mean(vector_k_0), mean(vector_k_1)),
                 "b" = c(sum(vector_k_0) / sum(vector_n_0), sum(vector_k_1) / sum(vector_n_1)),
                 "bb" = c(bb_optim(vector_n_0, vector_k_0)$par, bb_optim(vector_n_1, vector_k_1)$par),
                 "zib" = c(zib_optim(vector_n_0, vector_k_0)$par, zib_optim(vector_n_1, vector_k_1)$par),
                 "zibb" = c(zibb_optim(vector_n_0, vector_k_0)$par, zibb_optim(vector_n_1, vector_k_1)$par)
  )
        
  params <- rbind(params, tpar) 
  junc <- c(junc, splicing_key_list[i])

}

warnings()

rownames(params) <- junc


write.table(params, paste("../output/recount2/TCGA/param_matrix.recount2.", model, ".txt", sep = ""), sep = "\t", col.names = FALSE, quote = FALSE)




