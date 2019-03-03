library(VGAM)
library(tidyverse)

source("subscript_model/model_functions.R")

args <- commandArgs(trailingOnly = TRUE)
junc_file <- args[1]
pheno_file <- args[2]
param_file <- args[3]
out_file <- args[4]
model <- "zibb"

P <- read_tsv(param_file, col_names = FALSE)
P <- P %>% set_names("Splicing_Key", str_c("Param", 1:(ncol(P) - 1), sep = "_"))

D <- read_tsv(junc_file) 
D_param <- D %>% filter(Splicing_Key %in% P$Splicing_Key) %>%
    left_join(P, by = "Splicing_Key")
 
phenotype <- read_tsv(pheno_file)

scores <- c()

for(n in 1:length(phenotype$run)) {

  run_name <- phenotype$run[n]
  D_test<- D_param %>% filter(Sample_Name == run_name) 
  print(run_name)

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
                  "zibb" = cbind(unlist(map(D_test %>% transpose(), function(x) {get_prob_for_zibb(c(x$Param_1, x$Param_2, x$Param_3), x$Read_Count1 + x$Read_Count2, x$Read_Count1)})),
                                 unlist(map(D_test %>% transpose(), function(x) {get_prob_for_zibb(c(x$Param_4, x$Param_5, x$Param_6), x$Read_Count1 + x$Read_Count2, x$Read_Count1)}))
                               )
  )

  probs[probs < 1e-100] <- 1e-100
  probs <- probs / rowSums(probs)

  tlratios <- log(probs[,2]) - log(probs[,1])
  tlratios[tlratios < -20] <- -20
  tlratios[tlratios > 20] <- 20
  scores <- c(scores, sum(tlratios))
 
}

write.table(data.frame(run_name = phenotype$run, score = scores), out_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)




