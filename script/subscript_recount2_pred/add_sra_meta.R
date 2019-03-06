library(SRAdb)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

sra_dbname <- "~/db/SRAmetadb.sqlite"
sra_con <- dbConnect(dbDriver("SQLite"), sra_dbname)
D <- read_tsv(input_file, col_names = TRUE)

study_title <- c()
sample_attribute <- c() 
bases <- c()
sample_num <- c()
for(i in 1:nrow(D)) {
  rs <- getSRA(search_terms =D$Run[i], out_types=c("study", "sample", "run"), sra_con=sra_con)[1,]
  study_title <- c(study_title, rs$study_title)
  sample_attribute <- c(sample_attribute, rs$sample_attribute)
  bases <- c(bases, rs$bases)

  if (i %% 1000 == 0) print(i)
  # rs2 <- getSRA (search_terms =D$Study[i], out_types=c("study", "run"), sra_con=sra_con)
  # sample_num <- c(sample_num, nrow(rs2))
}

D2 <- cbind(D, bases, study_title, sample_attribute)
D2[,3] <- round(D2[,3], 3)

colnames(D2) <- c("Study_Accession", "Run_Accession", "SF3B1ness_Score", "Bases", "Study_Title", "Sample_Attribute")

write.table(D2 %>% filter(Bases >= 100000000), output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

