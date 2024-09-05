library(getopt,quietly = T)
command=matrix(c( 
  'help', 'h', 0,'loical', 'help information',
  'genotype', 'G', 1, 'character', 'Genotype file encoded with 012, Maker in line, sample in column.',
  'phenptype', 'P', 1, 'character', 'Phenotype file, sample in line, sample name should be same as in genotype file.',
  'max','M',1,'integer','Maximum SNP number keeped in model, default 20000.',
  'selection','S',1,'character','SNP selection method,can be rrBLUP (default) or ACC (based on random forest Accuracy, slowly).',
  'BestIndividuals','B',1,'character','The position of expected phenotype in whole phenotypic data set. Can be top,bottom or middle,defult top.'
),
byrow=T,ncol=5)
args=getopt(command)

if (!is.null(args$help) || is.null(args$genotype) || is.null(args$phenptype) ) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

## 设置默认值
if ( is.null(args$max) ) {
 args$max = 20000
}

if ( is.null(args$selection) ) {
  args$selection = "rrBLUP"
}
if ( is.null(args$BestIndividuals) ) {
  args$BestIndividuals = "top"
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

# Markers<-read.delim(file = genotype_file,header = T,col.names = 1)
# Markers<-as.matrix(Markers)

genotype <- fread(file = genotype_file,header = T)
# 删除genotype中包含空值的行
genotype <- genotype[complete.cases(genotype), ]

phenotype <- fread(file = phenotype_file,header = T)
# 删除phenotype中包含空值的行
phenotype <- phenotype[complete.cases(phenotype), ]

# 找到 "ENVs_Sample" 列的交集
common_samples <- intersect(genotype$ENVs_Sample, phenotype$ENVs_Sample)
  
# 根据交集子集重新赋值
genotype <- genotype[genotype$ENVs_Sample %in% common_samples, ]
phenotype <- phenotype[phenotype$ENVs_Sample %in% common_samples, ]
phenoall <- as.data.frame(phenotype)

Markers <- as.matrix(genotype)
row.names(Markers)<-genotype[[1]]
genoall <- as.data.frame(Markers)

Pheno <- as.matrix(phenotype[,-1])
row.names(Pheno)<-phenotype[[1]]

# Markers<-round(Markers)
n_m<-ncol(Markers)-1

for(n_pheno in 1:ncol(Pheno)){
  trait_name<-colnames(Pheno)[n_pheno]
  trait_dir<-paste0("./",trait_name,"/")
  if (!dir.exists(trait_dir))  { dir.create(trait_dir, recursive= T)}

  #已知预测已知的十折
  nfold=10

  cvSampleList <- cvSampleIndex(sampleNum = nrow(phenotype), cross=nfold, seed=1)
  ENVs_Sample <- data.frame(ENVs_Sample=phenotype$ENVs_Sample)
  
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

  modellist <- c("BayesA", "BayesB", "BRR", "rrBLUP","LASSO","SPLS","SVR","RFR")

  CVlist <- NULL
  
  for (x in 1: nfold){
    trainIdx <- cvSampleList[[x]]$trainIdx
    testIdx <- cvSampleList[[x]]$testIdx

    ENVs_Sample_train <- ENVs_Sample[trainIdx, ]
    ENVs_Sample_train <- as.data.frame(ENVs_Sample_train)
    colnames(ENVs_Sample_train) <- c("ENVs_Sample")

    ENVs_Sample_test <- ENVs_Sample[testIdx, ]
    ENVs_Sample_test <- as.data.frame(ENVs_Sample_test)
    colnames(ENVs_Sample_test) <- c("ENVs_Sample")

    #构建lightgbm所用数据

    gtrain <- as.data.frame(inner_join(genotype, ENVs_Sample_train, by = "ENVs_Sample"))
    gdtrain <- as.matrix(gtrain[,-1])
    row.names(gdtrain) <- gtrain[[1]]

    ptrain <- as.data.frame(inner_join(phenotype, ENVs_Sample_train, by = "ENVs_Sample"))
    pdtrain <- as.matrix(ptrain[,-1])
    row.names(pdtrain) <- ptrain[[1]]

    dtrain <- lgb.Dataset(data = as.matrix(gdtrain), label = pdtrain[,1])

    gtest <- as.data.frame(inner_join(genotype, ENVs_Sample_test, by = "ENVs_Sample"))
    gdtest <- as.matrix(gtest[,-1])
    row.names(gdtest) <- gtest[[1]]

    ptest <- as.data.frame(inner_join(phenotype, ENVs_Sample_test, by = "ENVs_Sample"))
    pdtest <- as.matrix(ptest[,-1])
    row.names(pdtest) <- ptest[[1]]

    dtest <- lgb.Dataset(data = as.matrix(gdtest), label = pdtest[,1])

    model <- lgb.train(
             params,
             data = dtrain,
             valids = list(test = dtest)
            )

    pred <- as.data.frame(predict(model, as.matrix(gdtest)))
    pred$ENVs_Sample <- rownames(pred)
    colnames(pred)[1] <- "lightgbm"

    #构建G2P所用数据
    gtrain <- inner_join(genoall, ENVs_Sample_train, by = "ENVs_Sample")
    rownames(gtrain) <- gtrain$ENVs_Sample
    gtrain <- subset(gtrain, select = -c(ENVs_Sample))
    gtrain <- as.data.frame(apply(gtrain, 1:2, as.numeric))

    gtest <- inner_join(genoall, ENVs_Sample_test, by = "ENVs_Sample")
    rownames(gtest) <- gtest$ENVs_Sample
    gtest <- subset(gtest, select = -c(ENVs_Sample))
    gtest <- as.data.frame(apply(gtest, 1:2, as.numeric))

    ptrain <- inner_join(phenoall, ENVs_Sample_train, by = "ENVs_Sample")
    rownames(ptrain) <- ptrain$ENVs_Sample
    ptrain <- subset(ptrain, select = -c(ENVs_Sample))
    ptrain <- as.data.frame(apply(ptrain, 1:2, as.numeric))

    ptest <- inner_join(phenoall, ENVs_Sample_test, by = "ENVs_Sample")
    rownames(ptest) <- ptest$ENVs_Sample
    ptest <- subset(ptest, select = -c(ENVs_Sample))
    ptest <- as.data.frame(apply(ptest, 1:2, as.numeric))

    gtrain <- as.matrix.data.frame(gtrain)
    gtest <- as.matrix.data.frame(gtest)
    ptrain <- as.matrix.data.frame(ptrain)
    ptest <- as.matrix.data.frame(ptest)

    CV1 <- G2P(trainMarker = gtrain, 
             trainPheno = ptrain[,n_pheno], 
             testMarker = gtest, 
             testPheno = ptest, 
             modelMethods = modellist, 
             outputModel = F)
    CV1 <- as.data.frame(CV1)
    CV1$ENVs_Sample <- rownames(CV1)

    CV <- inner_join(CV1, pred, by="ENVs_Sample")
    rownames(CV) <- CV$ENVs_Sample
    CV <- subset(CV, select = -c(ENVs_Sample))

    CVlist[[x]] <- CV    

  }     

  CVlist <- do.call(rbind, CVlist)

  evalTest <- evaluateGS(realScores = CVlist[,1], 
                         predScores = CVlist[,2:ncol(CVlist)], 
                         evalMethod = c( "pearson", "kendall","spearman", "RE",
                                           "Kappa",
                                           "AUC", "AUCpr", "NDCG", "meanNDCG",
                                           "MSE", "R2", "F1", "accuracy"), 
                         topAlpha = 1:90)

  CV.data<-as.data.frame(CVlist)

  write.csv(x = CV.data,
              file = paste0(trait_dir,"Select_",n_m,"_Markers_","pred.data.csv"),
              row.names = T)
  write.csv(x = evalTest$corMethosds,
              file = paste0(trait_dir,"Select_",n_m,"_Markers_","pred.data.evalTest.csv"),
              row.names = T)

  for (i in 2:ncol(CV.data)) {
      t.x<-0.2*range(CV.data[,1])[2]+0.8*range(CV.data[,1])[1]
      t.y<-0.9*range(CV.data[,1])[2]+0.1*range(CV.data[,1])[1]
      t.t<-paste0("R = ",sprintf("%1.3f",evalTest$corMethosds["pearson",names(CV.data)[i]]))
      p<-ggplot(data = CV.data,aes(x=CV.data[,1],y = CV.data[,i]))+
        geom_point(color="#3e90e2",pch=1)+
        geom_smooth(method = "lm",se = F,formula = y~x)+
        geom_abline(slope = 1,intercept = 0)+
        geom_text(x=t.x,y=t.y,label=t.t)+
        scale_x_continuous(limits = c(range(CV.data[,1])),name = "Obs value")+
        scale_y_continuous(limits = c(range(CV.data[,1])),name = paste0("Pred value of ",names(CV.data)[i]))+
        theme_bw()+
        theme(
          legend.position = "none",
          panel.grid = element_blank()
        )
      
    ggsave(filename = paste0(trait_dir,"Select_",n_m,"_Markers_","sectorPlot_for_",names(CV.data)[i],".pdf"),
             plot = p,width = 4,height = 4)
  }

  geno_split_cols <- strsplit(as.character(genoall$ENVs_Sample), "_")
  genoall$ENVs <- sapply(geno_split_cols, function(x) x[1])
  genoall$Sample <- sapply(geno_split_cols, function(x) x[2])

  pheno_split_cols <- strsplit(as.character(phenoall$ENVs_Sample), "_")
  phenoall$ENVs <- sapply(pheno_split_cols, function(x) x[1])
  phenoall$Sample <- sapply(pheno_split_cols, function(x) x[2])

  ENVs <- as.data.frame(unique(genoall$ENVs))
  Sample <- as.data.frame(unique(genoall$Sample))

  ENVs_count <- length(ENVs[,1])
  Sample_count <- length(Sample[,1])

  #十折验证基因型
  gcvSampleList <- cvSampleIndex(sampleNum = Sample_count, cross=10, seed=1)
  
  gmodellist <- c("BayesA", "BayesB", "BRR", "rrBLUP","LASSO","SPLS","SVR","RFR")

  gCVlist <- list()

  for( x in 1:nfold ){
  
  trainIdx <- gcvSampleList[[x]]$trainIdx
  testIdx <- gcvSampleList[[x]]$testIdx

  Sample_train <- Sample[trainIdx, ]
  Sample_train <- as.data.frame(Sample_train)
  colnames(Sample_train) <- c("Sample")

  Sample_test <- Sample[testIdx, ]
  Sample_test <- as.data.frame(Sample_test)
  colnames(Sample_test) <- c("Sample")

  #构建lightgbm的模型
  gtrain <- as.data.frame(inner_join(genoall, Sample_train, by = "Sample"))
  gdtrain <- as.matrix(subset(gtrain, select = -c(Sample, ENVs_Sample, ENVs)))
  row.names(gdtrain) <- gtrain[[1]]

  ptrain <- as.data.frame(inner_join(phenoall, Sample_train, by = "Sample"))
  pdtrain <- as.matrix(subset(ptrain, select = -c(Sample, ENVs_Sample, ENVs)))
  row.names(pdtrain) <- ptrain[[1]]

  dtrain <- lgb.Dataset(data = as.matrix(gdtrain), label = pdtrain[,1])

  gtest <- as.data.frame(inner_join(genoall, Sample_test, by = "Sample"))
  gdtest <- as.matrix(subset(gtest, select = -c(Sample, ENVs_Sample, ENVs)))
  row.names(gdtest) <- gtest[[1]]

  ptest <- as.data.frame(inner_join(phenoall, Sample_test, by = "Sample"))
  pdtest <- as.matrix(subset(ptest, select = -c(Sample, ENVs_Sample, ENVs)))
  row.names(pdtest) <- ptest[[1]]

  dtest <- lgb.Dataset(data = as.matrix(gdtest), label = pdtest[,1])

  model <- lgb.train(
             params,
             data = dtrain,
             valids = list(test = dtest)
            )

  pred <- as.data.frame(predict(model, as.matrix(gdtest)))
  pred$ENVs_Sample <- rownames(pred)
  colnames(pred)[1] <- "lightgbm"

  #构建G2P的相关模型

  gtrain <- inner_join(genoall, Sample_train, by = "Sample")
  rownames(gtrain) <- gtrain$ENVs_Sample
  gtrain <- subset(gtrain, select = -c(Sample, ENVs_Sample, ENVs))
  gtrain <- as.data.frame(apply(gtrain, 1:2, as.numeric))

  gtest <- inner_join(genoall, Sample_test, by = "Sample")
  rownames(gtest) <- gtest$ENVs_Sample
  gtest <- subset(gtest, select = -c(Sample, ENVs_Sample, ENVs))
  gtest <- as.data.frame(apply(gtest, 1:2, as.numeric))

  ptrain <- inner_join(phenoall, Sample_train, by = "Sample")
  rownames(ptrain) <- ptrain$ENVs_Sample
  ptrain <- subset(ptrain, select = -c(Sample, ENVs_Sample, ENVs))
  ptrain <- as.data.frame(apply(ptrain, 1:2, as.numeric))

  ptest <- inner_join(phenoall, Sample_test, by = "Sample")
  rownames(ptest) <- ptest$ENVs_Sample
  ptest <- subset(ptest, select = -c(Sample, ENVs_Sample, ENVs))
  ptest <- as.data.frame(apply(ptest, 1:2, as.numeric))

  gtrain <- as.matrix.data.frame(gtrain)
  gtest <- as.matrix.data.frame(gtest)
  ptrain <- as.matrix.data.frame(ptrain)
  ptest <- as.matrix.data.frame(ptest)

  #正式操作替换成以下模型
  # "BayesA", "BayesB", "BRR", "rrBLUP","LASSO","SPLS","SVR","RFR"

  gCV <- G2P(trainMarker = gtrain, 
             trainPheno = ptrain[,n_pheno], 
             testMarker = gtest, 
             testPheno = ptest, 
             modelMethods = gmodellist, 
             outputModel = F)
    
    gCV <- as.data.frame(gCV)
    gCV$ENVs_Sample <- rownames(gCV)

    gCV <- inner_join(gCV, pred, by="ENVs_Sample")
    rownames(gCV) <- gCV$ENVs_Sample
    gCV <- subset(gCV, select = -c(ENVs_Sample))

  gCVlist[[x]] <- gCV
  
}

# 将列表中的结果合并成一个数据框
gCVlist <- do.call(rbind, gCVlist)

gCVeval <- evaluateGS(realScores = gCVlist[,1], 
                      predScores = gCVlist[,2:ncol(gCVlist)], 
                      evalMethod = c( "pearson", "kendall","spearman", "RE",
                                           "Kappa",
                                           #"AUC", "AUCpr", "NDCG", "meanNDCG",
                                           "MSE", "R2", "F1", "accuracy"), 
                      topAlpha = 1:90)


gCV.data <- as.data.frame(gCVlist)

#调整trait_Dir
write.csv(x = gCV.data,
         file = paste0(trait_dir,"gCV_Select_",n_m,"_Markers_","pred.data.csv"),
         row.names = T)
write.csv(x = gCVeval$corMethosds,
          file = paste0(trait_dir,"gCV_Select_",n_m,"_Markers_","pred.data.evalTest.csv"),
          row.names = T)

gCV.data$ENVs_Sample <- rownames(gCV.data)
gCV_split_cols <- strsplit(as.character(gCV.data$ENVs_Sample), "_")
gCV.data$ENVs <- sapply(gCV_split_cols, function(x) x[1])
gCV.data <- gCV.data[, -which(names(gCV.data) == "ENVs_Sample")]

gCV.env <- data.frame(ENVs = unique(gCV.data$ENVs))

for (i in 2:(ncol(gCV.data) - 1)) {
  gCV.temp <- gCV.data[,c(1,i,ncol(gCV.data))]
  
  pearson_cor <- gCV.temp %>%
  group_by(ENVs) %>%
  summarize(pearson_cor = cor(realPhenScore, !!as.name(colnames(gCV.temp)[2])))
  
  pearson_cor <- as.data.frame(pearson_cor)

  gCV.temp <- inner_join(gCV.temp, pearson_cor, by="ENVs")
  gCV.temp$ENVs_R <- paste(gCV.temp$ENVs, "(", round(gCV.temp$pearson_cor, 3), ")", sep = "")

  colnames(pearson_cor)[2] <- names(gCV.temp)[2]
  gCV.env <- inner_join(gCV.env, pearson_cor, by="ENVs")

  t.x<-0.4*range(gCV.temp[,1])[2]+0.6*range(gCV.temp[,1])[1]
  t.y<-0.99*range(gCV.temp[,1])[2]
  t.t <- paste0("R = ", sprintf("%1.3f", gCVeval$corMethosds["pearson", names(gCV.temp)[2]]))
  
  p <- ggplot(data = gCV.temp, aes(x = gCV.temp[, 1], y = gCV.temp[, 2], color = ENVs_R)) +
    geom_point(pch = 1) +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    geom_text(x = t.x, y = t.y, label = t.t, hjust = 0, vjust = 1, color = "black", size = 5) +  # 调整文本位置
  
    scale_x_continuous(limits = c(range(gCV.temp[, 1])), name = "Obs value") +
    scale_y_continuous(limits = c(range(gCV.temp[, 2])), name =paste0("Pred value of ", names(gCV.temp)[2])) +
  
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = c(0.95, 0.05),  # 将图例放到右下角
      legend.justification = c(1, 0),  # 将图例靠右对齐并底部对齐
      legend.margin = margin(0, 0, 0, 0),  # 调整图例边距
      legend.key.size = unit(1, "lines"),  # 调整图例键的大小
      legend.title = element_blank()  # 删除图例的标题
    )

  ggsave(
    filename = paste0(trait_dir, "gCV_Select_", n_m, "_Markers_", "sectorPlot_for_", names(gCV.temp)[2], ".pdf"),
    plot = p, width = 6, height = 6  # 根据需要调整宽度和高度
  )
}

write.csv(x = gCV.env,
          file = paste0(trait_dir,"gCV_perENV_Select_",n_m,"_Markers_","pred.data.evalTest.csv"),
          row.names = F)

 #留一法交叉验证验证环境影响

  ecvSampleList <- cvSampleIndex(sampleNum = ENVs_count, cross=ENVs_count, seed=1)

  emodellist <- c("BayesA", "BayesB", "BRR", "rrBLUP","LASSO","SPLS","SVR","RFR")


  eCVlist <- list()
  for( x in 1:ENVs_count ){
  etrainIdx <- ecvSampleList[[x]]$trainIdx
  etestIdx <-  ecvSampleList[[x]]$testIdx

  ENVs_train <- ENVs[etrainIdx, ]
  ENVs_train <- as.data.frame(ENVs_train)
  colnames(ENVs_train) <- c("ENVs")

  ENVs_test <- ENVs[etestIdx, ]
  ENVs_test <- as.data.frame(ENVs_test)
  colnames(ENVs_test) <- c("ENVs")

 #构建lightgbm的模型
  gtrain <- as.data.frame(inner_join(genoall, ENVs_train, by = "ENVs"))
  gdtrain <- as.matrix(subset(gtrain, select = -c(Sample, ENVs_Sample, ENVs)))
  row.names(gdtrain) <- gtrain[[1]]

  ptrain <- as.data.frame(inner_join(phenoall, ENVs_train, by = "ENVs"))
  pdtrain <- as.matrix(subset(ptrain, select = -c(Sample, ENVs_Sample, ENVs)))
  row.names(pdtrain) <- ptrain[[1]]

  dtrain <- lgb.Dataset(data = as.matrix(gdtrain), label = pdtrain[,1])

  gtest <- as.data.frame(inner_join(genoall, ENVs_test, by = "ENVs"))
  gdtest <- as.matrix(subset(gtest, select = -c(Sample, ENVs_Sample, ENVs)))
  row.names(gdtest) <- gtest[[1]]

  ptest <- as.data.frame(inner_join(phenoall, ENVs_test, by = "ENVs"))
  pdtest <- as.matrix(subset(ptest, select = -c(Sample, ENVs_Sample, ENVs)))
  row.names(pdtest) <- ptest[[1]]

  dtest <- lgb.Dataset(data = as.matrix(gdtest), label = pdtest[,1])

  model <- lgb.train(
             params,
             data = dtrain,
             valids = list(test = dtest)
            )

  pred <- as.data.frame(predict(model, as.matrix(gdtest)))
  pred$ENVs_Sample <- rownames(pred)
  colnames(pred)[1] <- "lightgbm"

  #构建G2P的相关模型

  etrain <- inner_join(genoall, ENVs_train, by = "ENVs")
  rownames(etrain) <- etrain$ENVs_Sample
  etrain <- subset(etrain, select = -c(Sample, ENVs_Sample, ENVs))
  etrain <- as.data.frame(apply(etrain, 1:2, as.numeric))

  etest <- inner_join(genoall, ENVs_test, by = "ENVs")
  rownames(etest) <- etest$ENVs_Sample
  etest <- subset(etest, select = -c(Sample, ENVs_Sample, ENVs))
  etest <- as.data.frame(apply(etest, 1:2, as.numeric))

  eptrain <- inner_join(phenoall, ENVs_train, by = "ENVs")
  rownames(eptrain) <- eptrain$ENVs_Sample
  eptrain <- subset(eptrain, select = -c(Sample, ENVs_Sample, ENVs))
  eptrain <- as.data.frame(apply(eptrain, 1:2, as.numeric))

  eptest <- inner_join(phenoall, ENVs_test, by = "ENVs")
  rownames(eptest) <- eptest$ENVs_Sample
  eptest <- subset(eptest, select = -c(Sample, ENVs_Sample, ENVs))
  eptest <- as.data.frame(apply(eptest, 1:2, as.numeric))

  etrain <- as.matrix.data.frame(etrain)
  etest <- as.matrix.data.frame(etest)
  eptrain <- as.matrix.data.frame(eptrain)
  eptest <- as.matrix.data.frame(eptest)


  #正式操作替换成以下模型
  # "BayesA", "BayesB", "BRR", "rrBLUP","LASSO","SPLS","SVR","RFR"
  
  eCV <- G2P(trainMarker = etrain, 
                    trainPheno = eptrain[,n_pheno], 
                    testMarker = etest, 
                    testPheno = eptest, 
                    modelMethods = emodellist, 
                    outputModel = F)
  eCV <- as.data.frame(eCV)
  eCV$ENVs_Sample <- rownames(eCV)

  eCV <- inner_join(eCV, pred, by="ENVs_Sample")
  rownames(eCV) <- eCV$ENVs_Sample
  eCV <- subset(eCV, select = -c(ENVs_Sample))

  eCVlist[[x]] <- eCV
  
}

# 将列表中的结果合并成一个数据框
eCVlist <- do.call(rbind, eCVlist)

eCVeval <- evaluateGS(realScores = eCVlist[,1], 
                      predScores = eCVlist[,2:ncol(eCVlist)], 
                      evalMethod = c( "pearson", "kendall","spearman", "RE",
                                           "Kappa",
                                           #"AUC", "AUCpr", "NDCG", "meanNDCG",
                                           "MSE", "R2", "F1", "accuracy"), 
                      topAlpha = 1:90)

eCV.data <- as.data.frame(eCVlist)

#调整trait_Dir
write.csv(x = eCV.data,
         file = paste0(trait_dir,"eCV_Select_",n_m,"_Markers_","pred.data.csv"),
         row.names = T)
write.csv(x = eCVeval$corMethosds,
          file = paste0(trait_dir,"eCV_Select_",n_m,"_Markers_","pred.data.evalTest.csv"),
          row.names = T)

eCV.data$ENVs_Sample <- rownames(eCV.data)
eCV_split_cols <- strsplit(as.character(eCV.data$ENVs_Sample), "_")
eCV.data$ENVs <- sapply(eCV_split_cols, function(x) x[1])
eCV.data <- eCV.data[, -which(names(eCV.data) == "ENVs_Sample")]

eCV.env <- data.frame(ENVs = unique(eCV.data$ENVs))

for (i in 2:(ncol(eCV.data) - 1)) {
  eCV.temp <- eCV.data[,c(1,i,ncol(eCV.data))]
  
  pearson_cor <- eCV.temp %>%
  group_by(ENVs) %>%
  summarize(pearson_cor = cor(realPhenScore, !!as.name(colnames(eCV.temp)[2])))
  
  pearson_cor <- as.data.frame(pearson_cor)

  eCV.temp <- inner_join(eCV.temp, pearson_cor, by="ENVs")
  eCV.temp$ENVs_R <- paste(eCV.temp$ENVs, "(", round(eCV.temp$pearson_cor, 3), ")", sep = "")

  colnames(pearson_cor)[2] <- names(eCV.temp)[2]
  eCV.env <- inner_join(eCV.env, pearson_cor, by="ENVs")

  t.x<-0.4*range(eCV.temp[,1])[2]+0.6*range(eCV.temp[,1])[1]
  t.y<-0.99*range(eCV.temp[,1])[2]
  t.t <- paste0("R = ", sprintf("%1.3f", eCVeval$corMethosds["pearson", names(eCV.temp)[2]]))
  
  p <- ggplot(data = eCV.temp, aes(x = eCV.temp[, 1], y = eCV.temp[, 2], color = ENVs_R)) +
    geom_point(pch = 1) +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    geom_text(x = t.x, y = t.y, label = t.t, hjust = 0, vjust = 1, color = "black", size = 5) +  # 调整文本位置
  
    scale_x_continuous(limits = c(range(eCV.temp[, 1])), name = "Obs value") +
    scale_y_continuous(limits = c(range(eCV.temp[, 2])), name =paste0("Pred value of ", names(eCV.temp)[2])) +
  
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = c(0.95, 0.05),  # 将图例放到右下角
      legend.justification = c(1, 0),  # 将图例靠右对齐并底部对齐
      legend.margin = margin(0, 0, 0, 0),  # 调整图例边距
      legend.key.size = unit(1, "lines"),  # 调整图例键的大小
      legend.title = element_blank()  # 删除图例的标题
    )

  ggsave(
    filename = paste0(trait_dir, "eCV_Select_", n_m, "_Markers_", "sectorPlot_for_", names(eCV.temp)[2], ".pdf"),
    plot = p, width = 6, height = 6  # 根据需要调整宽度和高度
  )
}

write.csv(x = eCV.env,
          file = paste0(trait_dir,"eCV_perENV_Select_",n_m,"_Markers_","pred.data.evalTest.csv"),
          row.names = F)

 #十折交叉验证基因型和环境
  
  num=0
  allCVlist <- list() 

  # 定义所有模型列表
  allmodellist <- c("BayesA", "BayesB", "BRR", "rrBLUP","LASSO","SPLS","SVR","RFR")

  # 十折交叉验证基因型和环境
  for( x in 1:nfold ){
    for( y in 1:ENVs_count ){
      
      num <- num+1

    gtrainIdx <- gcvSampleList[[x]]$trainIdx
    gtestIdx <-  gcvSampleList[[x]]$testIdx
    etrainIdx <- ecvSampleList[[y]]$trainIdx
    etestIdx <- ecvSampleList[[y]]$testIdx

    Sample_train <- Sample[gtrainIdx, ]
    Sample_train <- as.data.frame(Sample_train)
    colnames(Sample_train) <- c("Sample")

    Sample_test <- Sample[gtestIdx, ]
    Sample_test <- as.data.frame(Sample_test)
    colnames(Sample_test) <- c("Sample")

    ENVs_train <- ENVs[etrainIdx, ]
    ENVs_train <- as.data.frame(ENVs_train)
    colnames(ENVs_train) <- c("ENVs")

    ENVs_test <- ENVs[etestIdx, ]
    ENVs_test <- as.data.frame(ENVs_test)
    colnames(ENVs_test) <- c("ENVs")

    combin_train <- expand.grid(ENVs = ENVs_train$ENVs, Sample = Sample_train$Sample)
    combin_test <- expand.grid(ENVs = ENVs_test$ENVs, Sample = Sample_test$Sample)

    all_train <- data.frame(ENVs_Sample = paste(combin_train$ENVs, combin_train$Sample, sep = "_"))
    all_test <- data.frame(ENVs_Sample = paste(combin_test$ENVs, combin_test$Sample, sep = "_"))

     #构建lightgbm的模型
  gtrain <- as.data.frame(inner_join(genoall, all_train, by = "ENVs_Sample"))
  gdtrain <- as.matrix(subset(gtrain, select = -c(Sample, ENVs_Sample, ENVs)))
  row.names(gdtrain) <- gtrain[[1]]

  ptrain <- as.data.frame(inner_join(phenoall, all_train, by = "ENVs_Sample"))
  pdtrain <- as.matrix(subset(ptrain, select = -c(Sample, ENVs_Sample, ENVs)))
  row.names(pdtrain) <- ptrain[[1]]

  dtrain <- lgb.Dataset(data = as.matrix(gdtrain), label = pdtrain[,1])

  gtest <- as.data.frame(inner_join(genoall, all_test, by = "ENVs_Sample"))
  gdtest <- as.matrix(subset(gtest, select = -c(Sample, ENVs_Sample, ENVs)))
  row.names(gdtest) <- gtest[[1]]

  ptest <- as.data.frame(inner_join(phenoall, all_test, by = "ENVs_Sample"))
  pdtest <- as.matrix(subset(ptest, select = -c(Sample, ENVs_Sample, ENVs)))
  row.names(pdtest) <- ptest[[1]]

  dtest <- lgb.Dataset(data = as.matrix(gdtest), label = pdtest[,1])

  model <- lgb.train(
             params,
             data = dtrain,
             valids = list(test = dtest)
            )

  pred <- as.data.frame(predict(model, as.matrix(gdtest)))
  pred$ENVs_Sample <- rownames(pred)
  colnames(pred)[1] <- "lightgbm"

  #构建G2P的相关模型

    alltrain <- inner_join(genoall, all_train, by = "ENVs_Sample")
    rownames(alltrain) <- alltrain$ENVs_Sample
    alltrain <- subset(alltrain, select = -c(Sample, ENVs_Sample, ENVs))
    alltrain <- as.data.frame(apply(alltrain, 1:2, as.numeric))

    alltest <- inner_join(genoall, all_test, by = "ENVs_Sample")
    rownames(alltest) <- alltest$ENVs_Sample
    alltest <- subset(alltest, select = -c(Sample, ENVs_Sample, ENVs))
    alltest <- as.data.frame(apply(alltest, 1:2, as.numeric))

    allptrain <- inner_join(phenoall, all_train, by = "ENVs_Sample")
    rownames(allptrain) <- allptrain$ENVs_Sample
    allptrain <- subset(allptrain, select = -c(Sample, ENVs_Sample, ENVs))
    allptrain <- as.data.frame(apply(allptrain, 1:2, as.numeric))

    allptest <- inner_join(phenoall, all_test, by = "ENVs_Sample")
    rownames(allptest) <- allptest$ENVs_Sample
    allptest <- subset(allptest, select = -c(Sample, ENVs_Sample, ENVs))
    allptest <- as.data.frame(apply(allptest, 1:2, as.numeric))

    alltrain <- as.matrix.data.frame(alltrain)
    alltest <- as.matrix.data.frame(alltest)
    allptrain <- as.matrix.data.frame(allptrain)
    allptest <- as.matrix.data.frame(allptest)

    allCV <- G2P(trainMarker = alltrain, 
                  trainPheno = allptrain[,n_pheno], 
                  testMarker = alltest, 
                  testPheno = allptest, 
                  modelMethods = allmodellist, 
                  outputModel = F)
  allCV <- as.data.frame(allCV)
  allCV$ENVs_Sample <- rownames(allCV)

  allCV <- inner_join(allCV, pred, by="ENVs_Sample")
  rownames(allCV) <- allCV$ENVs_Sample
  allCV <- subset(allCV, select = -c(ENVs_Sample))

  allCVlist[[num]] <- allCV  # 将每次循环的结果存储在列表中
  
 } 
}

# 停止并行处理

allCVlist <- do.call(rbind, allCVlist)

allCVeval <- evaluateGS(realScores = allCVlist[,1], 
                      predScores = allCVlist[,2:ncol(allCVlist)], 
                      evalMethod = c( "pearson", "kendall","spearman", "RE",
                                           "Kappa",
                                           #"AUC", "AUCpr", "NDCG", "meanNDCG",
                                           "MSE", "R2", "F1", "accuracy"), 
                      topAlpha = 1:90)


allCV.data <- as.data.frame(allCVlist)

#调整trait_Dir
write.csv(x = allCV.data,
         file = paste0(trait_dir,"allCV_Select_",n_m,"_Markers_","pred.data.csv"),
         row.names = T)
write.csv(x = allCVeval$corMethosds,
          file = paste0(trait_dir,"allCV_Select_",n_m,"_Markers_","pred.data.evalTest.csv"),
          row.names = T)

allCV.data$ENVs_Sample <- rownames(allCV.data)
allCV_split_cols <- strsplit(as.character(allCV.data$ENVs_Sample), "_")
allCV.data$ENVs <- sapply(allCV_split_cols, function(x) x[1])
allCV.data <- allCV.data[, -which(names(allCV.data) == "ENVs_Sample")]
allCV.env <- data.frame(ENVs = unique(allCV.data$ENVs))

for (i in 2:(ncol(allCV.data) - 1)) {
  allCV.temp <- allCV.data[,c(1,i,ncol(allCV.data))]
  
  pearson_cor <- allCV.temp %>%
  group_by(ENVs) %>%
  summarize(pearson_cor = cor(realPhenScore, !!as.name(colnames(allCV.temp)[2])))
  
  pearson_cor <- as.data.frame(pearson_cor)

  allCV.temp <- inner_join(allCV.temp, pearson_cor, by="ENVs")
  allCV.temp$ENVs_R <- paste(allCV.temp$ENVs, "(", round(allCV.temp$pearson_cor, 3), ")", sep = "")

  colnames(pearson_cor)[2] <- names(allCV.temp)[2]
  allCV.env <- inner_join(allCV.env, pearson_cor, by="ENVs")

  t.x<-0.4*range(allCV.temp[,1])[2]+0.6*range(allCV.temp[,1])[1]
  t.y<-0.99*range(allCV.temp[,1])[2]
  t.t <- paste0("R = ", sprintf("%1.3f", allCVeval$corMethosds["pearson", names(allCV.temp)[2]]))
  
  p <- ggplot(data = allCV.temp, aes(x = allCV.temp[, 1], y = allCV.temp[, 2], color = ENVs_R)) +
    geom_point(pch = 1) +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    geom_text(x = t.x, y = t.y, label = t.t, hjust = 0, vjust = 1, color = "black", size = 5) +  # 调整文本位置
  
    scale_x_continuous(limits = c(range(allCV.temp[, 1])), name = "Obs value") +
    scale_y_continuous(limits = c(range(allCV.temp[, 2])), name =paste0("Pred value of ", names(allCV.temp)[2])) +
  
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = c(0.95, 0.05),  # 将图例放到右下角
      legend.justification = c(1, 0),  # 将图例靠右对齐并底部对齐
      legend.margin = margin(0, 0, 0, 0),  # 调整图例边距
      legend.key.size = unit(1, "lines"),  # 调整图例键的大小
      legend.title = element_blank()  # 删除图例的标题
    )

  ggsave(
    filename = paste0(trait_dir, "allCV_Select_", n_m, "_Markers_", "sectorPlot_for_", names(allCV.temp)[2], ".pdf"),
    plot = p, width = 6, height = 6  # 根据需要调整宽度和高度
  )
}

write.csv(x = allCV.env,
          file = paste0(trait_dir,"allCV_perENV_Select_",n_m,"_Markers_","pred.data.evalTest.csv"),
          row.names = F)

}