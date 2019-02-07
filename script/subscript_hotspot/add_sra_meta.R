library(SRAdb)
library(tidyverse)

sra_dbname <- "~/db/SRAmetadb.sqlite"
sra_con <- dbConnect(dbDriver("SQLite"), sra_dbname)
D <- read_tsv("../output/sra/hotspot/SF3B1.hotspot.result.info.txt", col_names=FALSE)

study_title <- c()
sample_attribute <- c() 
bases <- c()
sample_num <- c()
for(i in 1:nrow(D)) {
  rs <- getSRA(search_terms =D$X2[i], out_types=c("study", "sample", "run"), sra_con=sra_con)[1,]
  study_title <- c(study_title, rs$study_title)
  sample_attribute <- c(sample_attribute, rs$sample_attribute)
  bases <- c(bases, rs$bases)

  rs2 <- getSRA (search_terms =D$X1[i], out_types=c("study", "run"), sra_con=sra_con)
  sample_num <- c(sample_num, nrow(rs2))
}

D2 <- cbind(D, bases, study_title, sample_num, sample_attribute)
D2[,3] <- round(D2[,3], 3)

colnames(D2) <- c("Study_Accession", "Run_Accession", "SF3B1ness_Score", "Mutation_Info", "Bases", "Study_Title", "Study_Sample_Number", "Sample_Attribute")

write.table(D2, "../output/sra/hotspot/SF3B1.hotspot.result.info2.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

