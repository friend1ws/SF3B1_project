library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
model <- args[1]

source("subscript_zibb/zibb_functions.R")


D <- read_tsv("../output/recount2/TCGA/TCGA_recount_SF3B1_junction_summary.txt")

splicing_key_list <- D %>% group_by(Splicing_Key) %>% 
    summarize(Max_Read_Count1 = max(Read_Count1), 
              Max_Read_Count2 = max(Read_Count2), 
              Max_Read_Ratio1 = max(Read_Count1 / (Read_Count1 + Read_Count2 + 1)),
              Mean_Read_Ratio1 = mean(Read_Count1 / (Read_Count1 + Read_Count2 + 1))) %>%
    filter(Max_Read_Count1 >= 3, Max_Read_Count2 >= 3, Max_Read_Ratio1 >= 0.05) %>% .$Splicing_Key

all_sample_name <- unique(D$Sample_Name)

Pred <- data.frame(sample = c(), cancer_type = c(), mutation_info = c(), score = c())
Params <- c()
for (k in 1:2) {
  
    target_sample <- D %>% 
        filter(grepl("SF3B1:K700|SF3B1:K666|SF3B1:H662|SF3B1:N626|SF3B1:R625|SF3B1:E622|SF3B1:G740|SF3B1:K741|SF3B1:G742|None", Mutation_Info)) %>%
        .$Sample_Name %>% unique()

    test_sample_name <- all_sample_name[seq(k, length(all_sample_name), 2)]
    training_sample_name <- setdiff(all_sample_name, test_sample_name) 

    # removing those having SF3B1 VUS mutations
    training_sample_name <- intersect(training_sample_name, target_sample)
  
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

  
    P <- as_tibble(params) %>% set_names(str_c("Param", 1:ncol(params), sep = "_")) %>% add_column(Splicing_Key = junc)
  
    D_param <- D %>% filter(Splicing_Key %in% junc) %>%
        select(Sample_Name, Cancer_Type, Splicing_Key, Read_Count1, Read_Count2, Mutation_Info) %>% 
        left_join(P, by = "Splicing_Key")
  
    scores <- c()
    for(n in 1:length(test_sample_name)) {
    
        D_test <- D_param %>% filter(Sample_Name %in% test_sample_name[n])
        sample_name <- unique(D_test %>% .$Sample_Name)
        mutation_info <- unique(D_test %>% .$Mutation_Info)
        cancer_type <- unique(D_test %>% .$Cancer_Type)
    
        probs <- switch(model,
                        "p" = cbind(unlist(map(D_test %>% transpose(), function(x) {dpois(x$Read_Count1, x$Param_1)})),
                                    unlist(map(D_test %>% transpose(), function(x) {dpois(x$Read_Count1, x$Param_2)}))
                                    ),
                        "b" = cbind(unlist(map(D_test %>% transpose(), function(x) {dbinom(x$Read_Count1, x$Read_Count1 + x$Read_Count2, x$Param_1)})),
                                    unlist(map(D_test %>% transpose(), function(x) {dbinom(x$Read_Count1, x$Read_Count1 + x$Read_Count2, x$Param_2)}))
                                    ),
                        "bb" = cbind(unlist(map(D_test %>% transpose(), function(x) {get_prob_for_bb(c(x$Param_1, x$Param_2), x$Read_Count1 + x$Read_Count2, x$Read_Count1)})),
                                     unlist(map(D_test %>% transpose(), function(x) {get_prob_for_bb(c(x$Param_3, x$Param_4), x$Read_Count1 + x$Read_Count2, x$Read_Count1)}))
                                     ),
                        "zib" = cbind(unlist(map(D_test %>% transpose(), function(x) {get_prob_for_zib(c(x$Param_1, x$Param_2), x$Read_Count1 + x$Read_Count2, x$Read_Count1)})),
                                      unlist(map(D_test %>% transpose(), function(x) {get_prob_for_zib(c(x$Param_3, x$Param_4), x$Read_Count1 + x$Read_Count2, x$Read_Count1)}))
                                    ),
                        "zibb" = cbind(unlist(map(D_test %>% transpose(), function(x) {get_prob_for_zibb(c(x$Param_1, x$Param_2, x$Param_3), x$Read_Count1 + x$Read_Count2, x$Read_Count1)})),
                                       unlist(map(D_test %>% transpose(), function(x) {get_prob_for_zibb(c(x$Param_4, x$Param_5, x$Param_6), x$Read_Count1 + x$Read_Count2, x$Read_Count1)}))
                                    )
        )
                            
        probs <- probs / rowSums(probs)

        tlratios <- log(probs[,2]) - log(probs[,1])
        tlratios[tlratios < -20] <- -20
        tlratios[tlratios > 20] <- 20
        tlratios[is.na(tlratios)] <- 0
        scores <- c(scores, sum(tlratios))
  
        if (scores[n] > 0 | mutation_info != "None") print(c(n, cancer_type, sample_name, mutation_info, scores[n]))
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


write.table(params, paste("../output/recount2/TCGA/param_matrix.recount2.", model, ".txt", sep = ""), sep = "\t", col.names = FALSE, quote = FALSE)

write.table(Pred, paste("../output/recount2/TCGA/TCGA.pred.result.recount2.", model, ".txt", sep = ""), sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)



