library(VGAM)
library(tidyverse)

source("subscript_model/model_functions.R")

args <- commandArgs(trailingOnly = TRUE)
junc_file <- args[1]
pheno_file <- args[2]
param_file <- args[3]
out_file <- args[4]
model <- "zibb"

# tparams <- read.table(param_file, sep = "\t")
# tparams <- read_tsv(param_file)

# print(dim(tparams))
# params <- tparams[,2:7]
# row.names(params) <- tparams[,1]
P <- read_tsv(param_file, col_names = FALSE)
P <- P %>% set_names("Splicing_Key", str_c("Param", 1:(ncol(P) - 1), sep = "_"))


# phenotype <- read.table("SRP033115.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# D <- read.table("SRP033115.junction_coverage.filt.txt", sep = "\t", header = TRUE)

# D <- read.table(junc_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# phenotype <- read.table(pheno_file, sep = "\t", quote = "", header = TRUE, comment.char = "@", stringsAsFactors = FALSE)
D <- read_tsv(junc_file) 
D_param <- D %>% filter(Splicing_Key %in% P$Splicing_Key) %>%
    left_join(P, by = "Splicing_Key")
 
phenotype <- read_tsv(pheno_file)

scores <- c()

for(n in 1:length(phenotype$run)) {

  run_name <- phenotype$run[n]
  D_test<- D_param %>% filter(Sample_Name == run_name) 
  print(run_name)

  # counts1 <- rep(0, nrow(params))
  # counts2 <- rep(0, nrow(params))
  # names(counts1) <- rownames(params)
  # names(counts2) <- rownames(params)
  # for(j in 1:nrow(tD)) {
  #   counts1[tD$Splicing_Key] <- tD$Read_Count1
  #   counts2[tD$Splicing_Key] <- tD$Read_Count2
  # }

  # probs <- c()
  # for (j in 1:length(counts1)) {
  #   tparams <- params[names(counts1)[j],]
  #   tparams1 <- c(as.numeric(tparams[1]), as.numeric(tparams[2]), as.numeric(tparams[3]))
  #   tparams2 <- c(as.numeric(tparams[4]), as.numeric(tparams[5]), as.numeric(tparams[6]))
  #   prob0 <- get_prob_for_zibb(tparams1, counts1[j] + counts2[j], counts1[j])
  #   prob1 <- get_prob_for_zibb(tparams2, counts1[j] + counts2[j], counts1[j])
  #   probs <- rbind(probs, c(prob0, prob1))
  # }
  # if (run_name == "SRR2313107") {
  #   print(dim(probs))
  #   print(cbind(as.numeric(probs[,1]), as.numeric(probs[,2])))
  # }

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




