# Sample_plasticity ----------------------------------------------------------------------

library(getopt,quietly = T)

command=matrix(c( 
  'help', 'h', 0,'loical', 'Compute sample plasticity of each trait.',
  'pheno', 'p', 1,'character', 'Phenotype file, text Use Tab separator, col name: ENVs Sample Trait1 [Trait2] ...'
),
byrow=T,ncol=5)
args=getopt(command)

if (!is.null(args$help) || is.null(args$pheno)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}
pheno_f<-args$pheno
library(data.table)
library(dplyr)
library("lme4")
library(FW) 
# library(ggplot2)
# pheno_f="Pheno_ENVs.txt"
trait.data<-read.delim(file = pheno_f,header = T)
vars<-names(trait.data)[c(-1,-2)]
sample.lm.all<-NULL
sample.lm1.all<-NULL
# trait="BW"
for (trait in vars) {
  print(paste0("Doing whith ",trait))
  lInd <- which(colnames(trait.data) == 'Sample')
  eInd <- which(colnames(trait.data) == 'ENVs')
  tInd <- which(colnames(trait.data) == trait)
  exp_trait <- trait.data[,c(lInd, eInd, tInd)]
  colnames(exp_trait) <- c("line_code","env_code","Yobs");
  exp_trait$Yobs <- as.numeric(exp_trait$Yobs)
  exp_trait <- aggregate(Yobs ~  line_code + env_code, data = exp_trait, mean,na.rm=T) ## To make sure only one phenotype record per each line in each environment
  exp_trait <- exp_trait[!is.na(exp_trait$Yobs),]
 
  line_codes <- unique(exp_trait$line_code)
  
  exp_trait_m <- exp_trait
  exp_trait_m$env_code <- as.factor(exp_trait_m$env_code)
  exp_trait_m$line_code <- as.factor(exp_trait_m$line_code)
  exp_trait_m$Yobs <- as.numeric(exp_trait_m$Yobs)
  lm_ <- FW(exp_trait_m$Yobs, exp_trait_m$line_code, exp_trait_m$env_code, method = "OLS")
  print(summary(lm_))
  #blue <- data.frame(ENVs = exp_trait_m$env_code, Sample = exp_trait_m$line_code, Yobs=lm_$y, meanY=lm_$yhat)
  blue <- data.frame(ENVs = exp_trait_m$env_code, Sample = exp_trait_m$line_code, Yobs=lm_$y)
  env <- data.frame(ENVs = lm_$ENVlevels, meanY=lm_$h)
  blue <- merge(blue, env, by.x = "ENVs", by.y="ENVs")

  sample.lm<-NULL
  sample.lm1<-NULL
  for (sample1 in unique(blue$Sample)) {
    trait_sample<-blue[blue$Sample==sample1,]
    lm(Yobs~meanY,data = trait_sample)->fit1
    summary(fit1)->fit.sum
    
    p.value<-1-pf(fit.sum[["fstatistic"]][[1]],
                  fit.sum[["fstatistic"]][[2]],
                  fit.sum[["fstatistic"]][[3]])
    
    Estimate<-fit.sum[["coefficients"]][,1]
    # pr<-round(-log(fit.sum[["coefficients"]][,4],10), digits = 4)
    RMSE<-sqrt(sum(residuals(fit1)^2)/fit1$df.residual)
    R2<-fit.sum[["r.squared"]]
    
    trait_sample$Intercept<-Estimate[1]
    trait_sample$slope<-Estimate[2]
    trait_sample$RMSE<-RMSE
    trait_sample$R2<-R2
    sample.lm1<-rbind(trait_sample,sample.lm1)
    
    fit.sum1<-data.frame(sample1,matrix(c(Estimate,RMSE,R2),nrow = 1,ncol = 4))
    names(fit.sum1)<-c("Sample","Intercept","slope","RMSE","R2")
    sample.lm<-rbind(fit.sum1,sample.lm)
  }
  sample.lm.all<-rbind(sample.lm.all,data.frame(Trait=trait,sample.lm))
  sample.lm1.all<-rbind(sample.lm1.all,data.frame(Trait=trait,sample.lm1))
  write.csv(x = sample.lm,file = paste0(trait,"_Sample_plasticity.csv"),row.names = F)
  write.csv(x = sample.lm1,file = paste0(trait,"_Sample_plasticity_with_data.csv"),row.names = F)
  
}

sample.lm.dcast<-sample.lm.all[,c(1:4)] %>% setDT() %>% dcast(Sample~Trait,value.var = c("Intercept", "slope"))
#write.csv(x = sample.lm.dcast,file = "Sample_plasticity_Intercept_slope.csv",row.names = F)
write.table(x = sample.lm.dcast,file = "Sample_plasticity_Intercept_slope.txt", quote = F, sep ="\t", row.names = F, col.names = TRUE)