library(VGAM)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
junc_file <- args[1]
pheno_file <- args[2]
param_file <- args[3]
out_file <- args[4]


get_prob_for_zibb <- function(params, s_n, s_k) {
  
  p_alpha <- params[1]
  p_beta <- params[2]
  p_pi <- params[3]

  if (s_k == 0) {
    return(p_pi + (1 - p_pi) * dbetabinom.ab(0, s_n, p_alpha, p_beta))
  } else {
    return((1 - p_pi) * dbetabinom.ab(s_k, s_n, p_alpha, p_beta))
  }
  
}


tparams <- read.table(param_file, sep = "\t")

params <- tparams[,2:7]
rownames(params) <- tparams[,1]


# phenotype <- read.table("SRP033115.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# D <- read.table("SRP033115.junction_coverage.filt.txt", sep = "\t", header = TRUE)
D <- read.table(junc_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
phenotype <- read.table(pheno_file, sep = "\t", quote = "", header = TRUE, comment.char = "@", stringsAsFactors = FALSE)


scores <- c()

for(n in 1:length(phenotype$run)) {
  
  run_name <- phenotype$run[n]
  tD <- D %>% filter(Sample_Name == run_name) 
  print(run_name)

  counts1 <- rep(0, nrow(params))
  counts2 <- rep(0, nrow(params))
  names(counts1) <- rownames(params)
  names(counts2) <- rownames(params)
  for(j in 1:nrow(tD)) {
    counts1[tD$Splicing_Key] <- tD$Read_Count1
    counts2[tD$Splicing_Key] <- tD$Read_Count2
  }

  probs <- c()
  for (j in 1:length(counts1)) {
    tparams <- params[names(counts1)[j],]
    tparams1 <- c(as.numeric(tparams[1]), as.numeric(tparams[2]), as.numeric(tparams[3]))
    tparams2 <- c(as.numeric(tparams[4]), as.numeric(tparams[5]), as.numeric(tparams[6]))
    prob0 <- get_prob_for_zibb(tparams1, counts1[j] + counts2[j], counts1[j])
    prob1 <- get_prob_for_zibb(tparams2, counts1[j] + counts2[j], counts1[j])
    probs <- rbind(probs, c(prob0, prob1))
  }
  # if (run_name == "SRR2313107") {
  #   print(dim(probs))
  #   print(cbind(as.numeric(probs[,1]), as.numeric(probs[,2])))
  # }
  probs[probs < 1e-100] <- 1e-100

  tlratios <- log(probs[,2]) - log(probs[,1])
  tlratios[tlratios < -20] <- -20
  tlratios[tlratios > 20] <- 20
  scores <- c(scores, sum(tlratios))
  # scores <- c(scores, colSums(log(probs))[2] - colSums(log(probs))[1])
 
}

write.table(data.frame(run_name = phenotype$run, score = scores), out_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)




