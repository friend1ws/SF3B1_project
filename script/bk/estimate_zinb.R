library(pscl)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
junc_file <- args[1]
pheno_file <- args[2]
out_file <- args[3]


get_prob_for_znb <- function(x, param, weight) {
  if (x == 0) {
    prob0 <- param$zero0 + (1 - param$zero0) * dnbinom(0, mu = weight * param$mu0, size = param$theta)
    prob1 <- param$zero1 + (1 - param$zero1) * dnbinom(0, mu = weight * param$mu1, size = param$theta)
  } else {
    prob0 <- (1 - param$zero0) * dnbinom(x, mu = weight * param$mu0, size = param$theta)
    prob1 <- (1 - param$zero1) * dnbinom(x, mu = weight * param$mu1, size = param$theta)
  }
  return(c(prob0, prob1))
}


tparams <- read.table("param_matrix.recount2.ninb.txt", sep = "\t")
params <- tparams[,2:6]
rownames(params) <- tparams[,1]

# phenotype <- read.table("SRP033115.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# D <- read.table("SRP033115.junction_coverage.filt.txt", sep = "\t", header = TRUE)
D <- read.table(junc_file, sep = "\t", header = TRUE)
phenotype <- read.table(pheno_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)



scores <- c()

for(n in 1:length(phenotype$run)) {
  print(n)
  run_name <- phenotype$run[n]
  tD <- D %>% filter(Sample_Name == run_name) 
  tweight <- phenotype$mapped_read_count[n] / 10000000
  
  counts <- rep(0, nrow(params))
  names(counts) <- rownames(params)
  for(j in 1:nrow(tD)) {
    counts[tD$Splicing_Key] <- tD$Read_Count
  }
  
  probs <- c()
  for (j in 1:length(counts)) {
    tparams <- params[names(counts)[j],]
    ttparam <- data.frame(mu0 = as.numeric(tparams[1]), mu1 = as.numeric(tparams[2]), zero0 = as.numeric(tparams[3]), zero1 = as.numeric(tparams[4]), theta = as.numeric(tparams[5]))
    probs <- rbind(probs, get_prob_for_znb(counts[j], ttparam, tweight))
  }
  scores <- c(scores, colSums(log(probs))[2] - colSums(log(probs))[1])
  
}

write.table(data.frame(run_name = phenotype$run, score = scores), out_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)




