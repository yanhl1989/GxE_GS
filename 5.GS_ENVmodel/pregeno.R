args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript pregeno.R <file_path>")
}

sample_list <- args[1]
geno_path <- args[2]
env_path <- args[3]

library("dplyr")

sample <- read.table(sample_list, header = T, sep = "\t")
geno <- read.table(geno_path, header=T, sep = ",", check.names = FALSE)
# 判断列名是否为 "NA"，如果是，则排除该列
geno <- geno[, !grepl("^NA$", colnames(geno))]
colnames(geno)[1] <- c("Sample")
env <- read.table(env_path, header = T, sep = "\t")
colnames(env)[1] <- c("ENVs")

result <- inner_join(sample, geno, by = "Sample")
result <- inner_join(result, env, by = "ENVs")
result$ENVs_Sample <- paste(result$ENVs, result$Sample, sep = "_")
result <- result[, c(ncol(result),3:(ncol(result)-1))]

write.table(result, file = "genotype.012", sep = ",", col.names = T, row.names = F, quote = F)
