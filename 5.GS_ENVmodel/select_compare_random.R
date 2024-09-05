library(getopt,quietly = T)
command=matrix(c( 
  'help', 'h', 0,'loical', 'help information',
  'genotype', 'G', 1, 'character', 'Genotype file encoded with 012, Maker in line, sample in column.',
  'phenptype', 'P', 1, 'character', 'Phenotype file, sample in line, sample name should be same as in genotype file.',
  'SNPcontribution', 'C', 1, 'character', 'SNP contribution list.' 
),
byrow=T,ncol=5)
args=getopt(command)

if (!is.null(args$help) || is.null(args$genotype) || is.null(args$phenptype) || is.null(args$SNPcontribution)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}


library(G2P,quietly = T)
library(lightgbm,quietly = T)
library(rrBLUP,quietly = T)
library(spls,quietly = T)
library(randomForest,quietly = T)
library(e1071,quietly = T)
library(brnn,quietly = T)
library(BGLR,quietly = T)
library(glmnet,quietly = T)
library(snowfall,quietly = T)
library(pROC,quietly = T)
library(PRROC,quietly = T)
library(ggplot2,quietly = T)
library(plotly,quietly = T)
library(data.table,quietly = T)
library(reshape,quietly = T)
library(dplyr)
source("/public/home/zzuaga04/hanxv/newNP600/5.GSmodel/pam/function.r")

# setDTthreads(threads = 4)

# genotype_file <- "z70.heapmap.012"
# phenotype_file <- "pheno_blup_FS.txt"

genotype_file <- args$genotype
phenotype_file <- args$phenptype
contribution_file <- args$SNPcontribution

# Markers<-read.delim(file = genotype_file,header = T,col.names = 1)
# Markers<-as.matrix(Markers)

genotype <- fread(file = genotype_file,header = T)
colnames(genotype)[1] <- "Sample"
phenotype <- fread(file = phenotype_file,header = T)
colnames(phenotype)[1] <- "Sample"

# 找到 "Sample" 列的交集
common_samples <- intersect(genotype$Sample, phenotype$Sample)
  
# 根据交集子集重新赋值
genotype <- genotype[genotype$Sample %in% common_samples, ]
genotype <- as.data.frame(genotype)
phenotype <- phenotype[phenotype$Sample %in% common_samples, ]

trait_name<-colnames(phenotype)[2]
trait_dir<-paste0("./",trait_name,"/")
if (!dir.exists(trait_dir))  { dir.create(trait_dir, recursive= T)}

# 输入SNP贡献列表
contribution <- fread(file = contribution_file,header=T)
contribution <- as.data.frame(contribution[,c(1,2)])
colnames(contribution) <- c("featureid","featureGain_sum")
contribution$featureGain_addup <- NA
contribution[1,3] <- contribution[1,2]
for ( i in 2:nrow(contribution)){
  contribution[i,3] <- contribution[i-1,3]+ contribution[i,2]
}
contribution$num <- rownames(contribution)
K <- contribution[which(contribution$featureGain_sum>(0.01*max(contribution$featureGain_sum))),]
A <- max(K$featureGain_addup)
S <- min(K$featureGain_sum)

contribution$num <- as.numeric(contribution$num)

# 绘制贡献累计图
p <- ggplot(data = contribution, aes(x = num, y = featureGain_addup)) +  
  geom_line() +  # 添加折线  
  geom_hline(yintercept = A, linetype = "dashed", color = "#3e90e2") + # 添加虚线
  labs(x = "Num", y = "Feature Gain Addup", title = "Feature Gain Addup over Num") +  
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  scale_x_continuous(breaks = c(min(contribution$num), max(contribution$num))) +  # 设置 x 轴显示的刻度
  geom_point(color = "#3e90e2", pch = 1, alpha = 0)   # 移除数据点，设置 alpha 为 0 可以隐藏数据点

# 保存图形
ggsave(filename = paste0(trait_dir, "feature_gain_addup_plot.png"), plot = p, width = 8, height = 6, dpi = 600)

# 绘制贡献分布图
g <- ggplot(data = contribution, aes(x = num, y = featureGain_sum)) +  
  geom_line() +  # 添加折线  
  geom_hline(yintercept = S, linetype = "dashed", color = "#3e90e2") + # 添加虚线
  labs(x = "Num", y = "Feature Gain Sum", title = "Feature Gain Sum over Num") +  
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  scale_x_continuous(breaks = c(min(contribution$num), max(contribution$num))) +  # 设置 x 轴显示的刻度
  geom_point(color = "#3e90e2", pch = 1, alpha = 0)   # 移除数据点，设置 alpha 为 0 可以隐藏数据点

# 保存图形
ggsave(filename = paste0(trait_dir, "feature_gain_sum_plot.png"), plot = g, width = 8, height = 6, dpi = 600)

# 设置cropgbm筛选snp的递增数量以及，随机筛选snp的递增数量
c_n_m <- nrow(contribution)
select_n_m<-c(1,2,4,8,16,32,64,128,200,300,500,1000)
select_n_m <- select_n_m[select_n_m < c_n_m]

g_n_m <- ncol(genotype)-1
random_n_m<-c(50,100,200,500,1000,3000,5000,10000)
random_n_m<-random_n_m[random_n_m < g_n_m]  

#用十折测试snp数量对精度的影响
nfold=10

cvSampleList <- cvSampleIndex(sampleNum = nrow(phenotype), cross=nfold, seed=1)
Sample <- data.frame(Sample=phenotype$Sample)
  
#lightgbm的相关参数
params <- list(  
  learning_rate = 0.1,  
  num_leaves = 10,  
  min_data_in_leaf = 1,  
  objective = "regression",  
  max_depth = -1,  
  feature_fraction = 1,  
  nrounds = 100, 
  early_stopping_rounds = 20 
 )

# SNP筛选后，数量对于模型精度的影响
select_result_sum <- data.frame(SNP_num = as.numeric(), pearson = as.numeric())

for ( i in select_n_m){

  select_list <- contribution[1:i,1]
  geno <- genotype[,c("Sample",select_list)]

  CVlist <- list()

   for (x in 1: nfold){
    trainIdx <- cvSampleList[[x]]$trainIdx
    testIdx <- cvSampleList[[x]]$testIdx

    Sample_train <- Sample[trainIdx, ]
    Sample_train <- as.data.frame(Sample_train)
    colnames(Sample_train) <- c("Sample")

    Sample_test <- Sample[testIdx, ]
    Sample_test <- as.data.frame(Sample_test)
    colnames(Sample_test) <- c("Sample")

    #构建lightgbm所用数据

    gtrain <- as.data.frame(inner_join(geno, Sample_train, by = "Sample"))
    gdtrain <- as.matrix(gtrain[,-1])
    row.names(gdtrain) <- gtrain[[1]]

    ptrain <- as.data.frame(inner_join(phenotype, Sample_train, by = "Sample"))
    pdtrain <- as.matrix(ptrain[,-1])
    row.names(pdtrain) <- ptrain[[1]]

    dtrain <- lgb.Dataset(data = as.matrix(gdtrain), label = pdtrain[,1])

    gtest <- as.data.frame(inner_join(geno, Sample_test, by = "Sample"))
    gdtest <- as.matrix(gtest[,-1])
    row.names(gdtest) <- gtest[[1]]

    ptest <- as.data.frame(inner_join(phenotype, Sample_test, by = "Sample"))
    pdtest <- as.matrix(ptest[,-1])
    row.names(pdtest) <- ptest[[1]]

    dtest <- lgb.Dataset(data = as.matrix(gdtest), label = pdtest[,1])

    model <- lgb.train(
             params,
             data = dtrain,
             valids = list(test = dtest)
            )

    pred <- as.data.frame(predict(model, as.matrix(gdtest)))
    pred$Sample <- rownames(pred)
    colnames(pred)[1] <- "lightgbm"

    CV <- inner_join(phenotype, pred, by="Sample")

    CVlist[[x]] <- CV    

  }   
  
  CVlist <- do.call(rbind, CVlist)
  CVlist <- as.data.frame(CVlist)
  CVeval <- evaluateGS(realScores = CVlist[,2], 
                      predScores = CVlist[,3], 
                      evalMethod = c( "pearson", "kendall","spearman", "RE",
                                           "Kappa",
                                           #"AUC", "AUCpr", "NDCG", "meanNDCG",
                                           "MSE", "R2", "F1", "accuracy"), 
                      topAlpha = 1:90)
  
  select_result <- data.frame(SNP_num=i, pearson= CVeval$corMethosds["pearson",1])
  rownames(select_result) <- select_result$SNP_num

  select_result_sum <- rbind(select_result_sum, select_result)

}


# 找到最大r值的对应的数量
max <- select_result_sum[which.max(select_result_sum$pearson), 1]

# 细化峰值附近的情况
select_num <- c(select_n_m, seq(max-5,max+5, by=1 ))
select_num <- sort(unique(select_num))
select_num <- select_num[select_num > 0 & select_num < c_n_m]

# SNP筛选细化后，数量对于模型精度的影响
select_result_thin_sum <- data.frame(SNP_num = as.numeric(), pearson = as.numeric())

for ( i in select_num){

  select_list <- contribution[1:i,1]
  geno <- genotype[,c("Sample",select_list)]

  CVlist <- list()

   for (x in 1: nfold){
    trainIdx <- cvSampleList[[x]]$trainIdx
    testIdx <- cvSampleList[[x]]$testIdx

    Sample_train <- Sample[trainIdx, ]
    Sample_train <- as.data.frame(Sample_train)
    colnames(Sample_train) <- c("Sample")

    Sample_test <- Sample[testIdx, ]
    Sample_test <- as.data.frame(Sample_test)
    colnames(Sample_test) <- c("Sample")

    #构建lightgbm所用数据

    gtrain <- as.data.frame(inner_join(geno, Sample_train, by = "Sample"))
    gdtrain <- as.matrix(gtrain[,-1])
    row.names(gdtrain) <- gtrain[[1]]

    ptrain <- as.data.frame(inner_join(phenotype, Sample_train, by = "Sample"))
    pdtrain <- as.matrix(ptrain[,-1])
    row.names(pdtrain) <- ptrain[[1]]

    dtrain <- lgb.Dataset(data = as.matrix(gdtrain), label = pdtrain[,1])

    gtest <- as.data.frame(inner_join(geno, Sample_test, by = "Sample"))
    gdtest <- as.matrix(gtest[,-1])
    row.names(gdtest) <- gtest[[1]]

    ptest <- as.data.frame(inner_join(phenotype, Sample_test, by = "Sample"))
    pdtest <- as.matrix(ptest[,-1])
    row.names(pdtest) <- ptest[[1]]

    dtest <- lgb.Dataset(data = as.matrix(gdtest), label = pdtest[,1])

    model <- lgb.train(
             params,
             data = dtrain,
             valids = list(test = dtest)
            )

    pred <- as.data.frame(predict(model, as.matrix(gdtest)))
    pred$Sample <- rownames(pred)
    colnames(pred)[1] <- "lightgbm"

    CV <- inner_join(phenotype, pred, by="Sample")

    CVlist[[x]] <- CV    

  }   
  
  CVlist <- do.call(rbind, CVlist)
  CVlist <- as.data.frame(CVlist)
  CVeval <- evaluateGS(realScores = CVlist[,2], 
                      predScores = CVlist[,3], 
                      evalMethod = c( "pearson", "kendall","spearman", "RE",
                                           "Kappa",
                                           #"AUC", "AUCpr", "NDCG", "meanNDCG",
                                           "MSE", "R2", "F1", "accuracy"), 
                      topAlpha = 1:90)
  
  select_result_thin <- data.frame(SNP_num=i, pearson= CVeval$corMethosds["pearson",1])
  rownames(select_result_thin) <- select_result_thin$SNP_num

  select_result_thin_sum <- rbind(select_result_thin_sum, select_result_thin)

}

write.table(select_result_thin_sum, file=paste0(trait_dir, "Improved_SNP_Cropgbm_select_accuracy.txt"),
            sep="\t", row.names = F, col.names=T, quote=F)

# 找到最优snp组合
select_best_list <- contribution[1:select_result_thin_sum[which.max(select_result_thin_sum$pearson), 1],1]
# geno_best <- genotype[,c("Sample",select_best_list)]
write.table(select_best_list, file = paste0(trait_dir, "genotype.list"), col.names=F, row.names=F, sep=",", quote=F)

# 找到最大值的行
select_max_row <- select_result_thin_sum[which.max(select_result_thin_sum$pearson), ]

# 创建绘图对象
s <- ggplot(data = select_result_thin_sum, aes(x = SNP_num, y = pearson)) +
  geom_point(color = "#3e90e2", pch = 1) +  # 移除数据点
  labs(x = "SNP number", y = "R") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  geom_line(color = "#3e90e2") +  # 添加拟合线，这里使用平滑曲线进行拟合
  geom_text(aes(x = (0.9 * SNP_num), y = pearson, label = paste("(", SNP_num, " ,", round(pearson, 3), " )")), 
            data = select_max_row, vjust = -0.5, hjust = 0.5, color = "black")   # 标记横坐标值



# 保存图像
ggsave(filename = paste0(trait_dir, "Improved_SNP_Cropgbm_select_accuracy.png"), plot = s, width = 8, height = 6, dpi = 600)


# 随机挑选，数量对于模型精度的影响
# 进行20次随机抽样，以去除特殊情况
random.seed <- seq(1, 20, by=1)

random_result_sum <- data.frame(SNP_num = as.numeric(), pearson = as.numeric())
snp_list <- colnames(genotype[,-1]) #去除Sample列

for ( j in random.seed){

set.seed(j)

for ( i in random_n_m){

  random_list <- sample(snp_list, i)
  geno <- genotype[,c("Sample",random_list)]

  CVlist <- list()

   for (x in 1: nfold){
    trainIdx <- cvSampleList[[x]]$trainIdx
    testIdx <- cvSampleList[[x]]$testIdx

    Sample_train <- Sample[trainIdx, ]
    Sample_train <- as.data.frame(Sample_train)
    colnames(Sample_train) <- c("Sample")

    Sample_test <- Sample[testIdx, ]
    Sample_test <- as.data.frame(Sample_test)
    colnames(Sample_test) <- c("Sample")

    #构建lightgbm所用数据

    gtrain <- as.data.frame(inner_join(geno, Sample_train, by = "Sample"))
    gdtrain <- as.matrix(gtrain[,-1])
    row.names(gdtrain) <- gtrain[[1]]

    ptrain <- as.data.frame(inner_join(phenotype, Sample_train, by = "Sample"))
    pdtrain <- as.matrix(ptrain[,-1])
    row.names(pdtrain) <- ptrain[[1]]

    dtrain <- lgb.Dataset(data = as.matrix(gdtrain), label = pdtrain[,1])

    gtest <- as.data.frame(inner_join(geno, Sample_test, by = "Sample"))
    gdtest <- as.matrix(gtest[,-1])
    row.names(gdtest) <- gtest[[1]]

    ptest <- as.data.frame(inner_join(phenotype, Sample_test, by = "Sample"))
    pdtest <- as.matrix(ptest[,-1])
    row.names(pdtest) <- ptest[[1]]

    dtest <- lgb.Dataset(data = as.matrix(gdtest), label = pdtest[,1])

    model <- lgb.train(
             params,
             data = dtrain,
             valids = list(test = dtest)
            )

    pred <- as.data.frame(predict(model, as.matrix(gdtest)))
    pred$Sample <- rownames(pred)
    colnames(pred)[1] <- "lightgbm"

    CV <- inner_join(phenotype, pred, by="Sample")

    CVlist[[x]] <- CV    

  }   
  
  CVlist <- do.call(rbind, CVlist)
  CVlist <- as.data.frame(CVlist)
  CVeval <- evaluateGS(realScores = CVlist[,2], 
                      predScores = CVlist[,3], 
                      evalMethod = c( "pearson", "kendall","spearman", "RE",
                                           "Kappa",
                                           #"AUC", "AUCpr", "NDCG", "meanNDCG",
                                           "MSE", "R2", "F1", "accuracy"), 
                      topAlpha = 1:90)
  
  random_result <- data.frame(SNP_num=i, pearson= CVeval$corMethosds["pearson",1])
  rownames(random_result) <- random_result$SNP_num

  random_result_sum <- rbind(random_result_sum, random_result)

}

}

write.table(random_result_sum, file=paste0(trait_dir, "Improved_SNP_random_select_accuracy.txt"),
            sep="\t", row.names = F, col.names=T, quote=F)

# 找到最大值的行
# 计算每个SNP_num对应的pearson均值
mean_pearson <- aggregate(pearson ~ SNP_num, data = random_result_sum, FUN = mean)

# 找到pearson均值最大的SNP_num
max_mean_row <- mean_pearson[which.max(mean_pearson$pearson), ]


# 创建绘图对象
r <- ggplot() +
  geom_boxplot(data = random_result_sum, aes(x = SNP_num, y = pearson, group = SNP_num), fill = "#367db8", color = "#367db8", alpha = 0.5) +  # 添加箱线图
  labs(x = "SNP number", y = "R") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
#  geom_smooth(method = "loess", color = "#3e90e2", se = FALSE) +  # 添加拟合线，这里使用线性模型进行拟合
  geom_line(data = mean_pearson, aes(x = SNP_num, y = pearson, color = "#367db8")) +
  geom_text(aes(x = (0.9 * SNP_num), y = pearson, label = paste("(", SNP_num, " ,", round(pearson, 3), " )")), 
            data = max_mean_row, vjust = -0.5, hjust = 0.5, color = "black")   # 标记横坐标值



# 保存图像
ggsave(filename = paste0(trait_dir, "Improved_SNP_random_select_accuracy.png"), plot = r, width = 8, height = 6, dpi = 600)


# 创建绘图对象
a <- ggplot() +
  geom_line(data = select_result_thin_sum, aes(x = SNP_num, y = pearson, color = "Cropgbm")) +  # 添加Cropgbm组的折线
  geom_point(data = select_result_thin_sum, aes(x = SNP_num, y = pearson), shape = 1, color = "#e4191c", alpha = 0.4) +  # 添加Cropgbm组的数据点，降低透明度
  geom_line(data = mean_pearson, aes(x = SNP_num, y = pearson, color = "Random")) +  # 添加Random组的平滑曲线，不显示置信区间
  geom_boxplot(data = random_result_sum, aes(x = SNP_num, y = pearson, group = SNP_num), fill = "white", color = "#367db8", alpha = 0.4) +  # 添加箱线图，降低透明度
  labs(x = "SNP number", y = "R", color = "SNP select method") +
  theme_bw() +
  theme(
    legend.position = "right",  # 显示图例并将其放置在右侧
    panel.grid = element_blank()
  ) +
  scale_color_manual(values = c("Cropgbm" = "#e4191c", "Random" = "#367db8")) +  # 设置颜色对应的标签
  geom_text(data = select_max_row, aes(x = (0.9 * SNP_num), y = pearson, label = paste("(", SNP_num, " ,", round(pearson, 3), " )")), 
            vjust = -0.5, hjust = 0.5, color = "black", size = 2.5) +  # 标记Cropgbm筛选最优值坐标值
  geom_text(data = max_mean_row, aes(x = (0.9 * SNP_num), y = pearson, label = paste("(", SNP_num, " ,", round(pearson, 3), " )")), 
            vjust = -0.5, hjust = 0.5, color = "black", size = 2.5)    # 标记随机筛选最优值坐标值


ggsave(filename = paste0(trait_dir, "Improved_SNP_cropgbm_and_random_select_accuracy.png"), plot = a, width = 8, height = 6, dpi = 600)
