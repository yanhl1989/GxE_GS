# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 检查是否有足够的参数传递给脚本
if (length(args) < 1) {
  stop("Usage: Rscript extract_ENV <file_env>")
}

library("dplyr")
# 获取第一个变量的值
file_env <- args[1]


File_env <- read.table(file_env, header = FALSE, sep = "\t")

result_df <- data.frame(env_code = c("GY10","QZ09","ZZ12","ZZ10","AY09","QZ08","ALE14","AY11","AY10","AY13","AY07","AY08","SHZ14","AY15","ALE15","LQ08","KEL14","SHZ15","AKS09","AY14","KEL15"))

for (i in 1:nrow(File_env)) {
   file_path <- File_env[i, 1]
   file <- read.table(file_path, header = TRUE, sep = "\t")
   env <- File_env[i, 2]
   if (!grepl("Dim", env)) {
     result_list <- data.frame(env_code = file$env_code, env = file[, env])
     colnames(result_list)[2] <- env
     result_df <- inner_join(result_df, result_list, by = "env_code")
   }
}


write.table(result_df, "output.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
