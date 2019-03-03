library(SF3B1ness)

args <- commandArgs(trailingOnly = TRUE)
junc_file <- args[1]
out_file <- args[2]

SF3B1ness_recout2 <- SF3B1ness_recount2(junc_file)
write.table(SF3B1ness_recout2, out_file, quote = FALSE, sep = "\t", row.names = FALSE)


