library(getopt)

command=matrix(c( 
  'help', 'h', 0,'loical', 'Please keep the internet connection.',
  'phenodata','p',1,'character','TXT file (Tab separator) whith colnames: ENVs Sample Pheno1 Pheno2 ...'
),
byrow=T,ncol=5)
args=getopt(command)

if (!is.null(args$help) || is.null(args$phenodata)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}


phenofile<-args$phenodata

if(is.null(phenofile)){
  phenofile<-"pheno.txt"
}

if (file.exists(phenofile)) {
   pheno.data<-read.delim(file = phenofile,header = T)
}else{
  print("Phenotype file not exit!")
  q(status = 1)
}

library(readxl)
library("lme4")
library(ggplot2)
library("dplyr")
library(data.table)
library(VIM)
library(mice)
library(ComplexHeatmap)
library(dendsort)
library(stringr)
# library(ggtree)

# phenotype variation analysis  -------------------------------------------
unique(pheno.data$ENVs)
head(pheno.data)

# pheno.data<-read_excel(path = "Phnotype_GWS_rep.xlsx")

vars<-c("BW","LP","SI","FL","FS","FM")
all.sum<-NULL
all.aov<-NULL
result <- data.frame(line=unique(pheno.data$Sample))
for (var in vars) {
  lInd <- which(colnames(pheno.data) == 'Sample')
  eInd <- which(colnames(pheno.data) == 'ENVs')
  tInd <- which(colnames(pheno.data) == var)
  exp_trait <- pheno.data[,c(lInd, eInd, tInd)]
  colnames(exp_trait)<-c("lInd", "eInd", "tInd")
  exp_trait <- na.omit(exp_trait)
  aov<-aov(tInd~lInd*eInd, data = exp_trait)
  sumaov<-summary(aov)
  # as.character(row.names(sumaov[[1]]))
  res.aov <- data.frame(variable=var,
                        grp=row.names(sumaov[[1]]),
                        sumaov[[1]])
  # names(res.aov)
  res.aov$vcov_p<-res.aov$Sum.Sq/sum(res.aov$Sum.Sq)*100
  
  all.aov <- rbind(all.aov,res.aov)
  ## 建模
  blup <- lmer(tInd~(1|lInd)+(1|eInd), data = exp_trait)
  ss<-summary(blup)
  #test内容
  blp <- ranef(blup)
  LINE=blp$lInd+blup@beta
  res=data.frame(id=rownames(LINE),blup=LINE)
  colnames(res)<-c("line",var)
  result <- inner_join(result, res, by = "line")
  #test结束

  ss1<-as.data.frame(ss[["varcor"]])[,c(1,4,5)]
  ss1$vcov_p<-ss1$vcov/sum(ss1$vcov)*100
  
  ss2<-data.frame(
    variable=var,
    ss1
  )
  all.sum<-rbind(all.sum,ss2)
}
write.table(result,file="pheno_blup.txt",row.names = F,quote = F,sep="\t")
write.csv(x = all.aov,file = "ANOVA.csv",row.names = F)
write.csv(x = all.sum,file = "lmer_Variance.csv",row.names = F)
all.sum$variable<-factor(x = all.sum$variable,levels = c("LYA","RYA","RYP","LYP","BW","LP","SI","LS","FL","FN","FS","FM","FE"))
all.sum$grp<-factor(x = all.sum$grp,
                    levels = c("lInd","eInd","Residual"),
                    labels = c("lInd","ENVs","Residual"))
p<-ggplot(data = all.sum,aes(x = variable,y = vcov_p,group=grp))+
  geom_bar(aes(fill=grp),
           position = "stack",stat = "identity",
           color="black",alpha=.9)+
  scale_x_discrete(name="")+
  scale_y_continuous(name = "Variance precentge")+
  scale_fill_manual(values = c("#ef6c60","#00a0e9","#efbc60"),name="Effects")+
  theme_bw()+
  theme(axis.text.y =element_text(size= 10,color = "black"),
        axis.text.x =element_text(face= "bold",size= 10,color = "black"),
        axis.title = element_text(face= "bold",size= 10,color = "black"),
        panel.grid=element_blank(),
        legend.position = "right")
ggsave(plot = p,filename = "Variance_Line_ENVs.png",dpi=600,width=6,height=4)

all.aov1<-read.csv(file = "ANOVA.csv")

all.aov1$variable<-factor(x = all.sum$variable,levels = c("LYA","RYA","RYP","LYP","BW","LP","SI","LS","FL","FN","FS","FM","FE"))
# all.aov$grp <- row.names(all.aov)

all.aov1$grp<-factor(x = all.aov1$grp,
                    levels = c("lInd","eInd","Residuals"),
                    labels = c("lInd","ENVs","Residual"))

p<-ggplot(data = all.aov1,aes(x = variable,y = vcov_p,group=grp))+
  geom_bar(aes(fill=grp),
           position = "stack",stat = "identity",
           color="black",alpha=.9)+
  scale_x_discrete(name="")+
  scale_y_continuous(name = "Variance precentge")+
  scale_fill_manual(values = c("#ef6c60","#00a0e9","#efbc60"),name="Effects")+
  theme_bw()+
  theme(axis.text.y =element_text(size= 10,color = "black"),
        axis.text.x =element_text(face= "bold",size= 10,color = "black"),
        axis.title = element_text(face= "bold",size= 10,color = "black"),
        panel.grid=element_blank(),
        legend.position = "right")
ggsave(plot = p,filename = "ANOVA_Line_ENVs.png",dpi=600,width=6,height=4)
