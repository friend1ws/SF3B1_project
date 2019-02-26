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
for (k in 1:2) {
  
  target_D <- D %>% filter(grepl("SF3B1:K700|SF3B1:K666|SF3B1:H662|SF3B1:N626|SF3B1:R625|SF3B1:E622|SF3B1:G740|SF3B1:K741|SF3B1:G742|None", Mutation_Info))
  
  test_sample_name <- all_sample_name[seq(k, length(all_sample_name), 2)]
  training_sample_name <- setdiff(all_sample_name, test_sample_name) 

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
    # cret0 <- zibb_optim(vector_n[class_ind == 0], vector_k[class_ind == 0])
    # cret1 <- zibb_optim(vector_n[class_ind == 1], vector_k[class_ind == 1])

  }

  rownames(params) <- junc


  scores <- c()
  for(n in 1:length(test_sample_name)) {
    tD <- D %>% filter(Sample_Name %in% test_sample_name[n])
  
    probs <- c()
    for (j in 1:length(junc)) {
      ttD <- tD %>% filter(Splicing_Key == rownames(params)[j])
      if (nrow(ttD) == 1) {
        prob0 <- get_prob_for_zibb(params[j, 1:3], ttD$Read_Count1 + ttD$Read_Count2, ttD$Read_Count1)
        prob1 <- get_prob_for_zibb(params[j, 4:6], ttD$Read_Count1 + ttD$Read_Count2, ttD$Read_Count1)
      
        probs <- rbind(probs, c(prob0, prob1))
      }
    }
  
    tlratios <- log(probs[,2]) - log(probs[,1])
    tlratios[tlratios < -20] <- -20
    tlratios[tlratios > 20] <- 20
    scores <- c(scores, sum(tlratios))
  
    # scores <- c(scores, colSums(log(probs))[2] - colSums(log(probs))[1])
    if (scores[n] > 0 | tD[1,]$Mutation_Info != "None") print(as.character(c(tD[1,1:4], scores[n])))
  }

  
  ssD <- data.frame(Score = scores)
  ssD$Sample_Name <- as.character(test_sample_name)

  sD <- D %>% select(Cancer_Type, Sample_Name, Mutation_Info) %>% unique() %>% filter(Sample_Name %in% test_sample_name) %>% left_join(ssD, by = "Sample_Name")

  tA <- data.frame(sample = sD$Sample_Name,
                  cancer_type = sD$Cancer_Type,
                  mutation_info = sD$Mutation_Info,
                  score = sD$Score)

  Pred <- rbind(Pred, tA)
}


write.table(params, "../output/recount2/TCGA/param_matrix.recount2.zibb.txt", sep = "\t", col.names = FALSE, quote = FALSE)

write.table(Pred, "../output/recount2/TCGA/TCGA.pred.result.recount2.zibb.txt", sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)



