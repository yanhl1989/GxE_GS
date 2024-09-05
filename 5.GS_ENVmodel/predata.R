args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript predata.R <file_path>")
}

file_path <- args[1]


data <- read.table(file_path, header = T , sep = "\t")

sample <- data[,c(1,2)]

write.table(sample, file = "ENVs_Sample.list", sep = "\t", col.names = T, row.names = F, quote = F)

data$ENVs_Sample <- paste(data$ENVs, data$Sample, sep = "_")

data <- data[, -c(1,2)]

for( i in 1:(ncol(data)-1) ){

    filename <- colnames(data)[i]

    trait <- data[, c(ncol(data), i) ]

    write.table(trait, file = paste(filename, "trait", sep = "."), sep="\t", col.names = T, row.names = F, quote = F)

} 