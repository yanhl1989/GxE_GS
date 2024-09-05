library(getopt)

command=matrix(c( 
  'help', 'h', 0,'loical', 'Please keep the internet connection.',
  'pheno', 'p', 1,'character', 'Phenotype file, text Use Tab separator, col name: ENVs Sample Trait1 [Trait2] ...',
  'siteinfo','s',1,'character','Excel (xlsx) file whith colnames: Site Year Envs LON LAT SeedingDate',
  'metdata', 'm', 1,'character', 'Met data from prepareMet4search.R',
  'metvar' , 'v',0,'character', 'Met Variables used in ENV search, character can be: "All"(default) or specified one or more variable name in Metdata file like: "radn,dh,tmean,tmax,tmin"'
),
byrow=T,ncol=5)
args=getopt(command)

if (!is.null(args$help) || is.null(args$pheno) || is.null(args$siteinfo) || is.null(args$metdata)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

library(corrgram,quietly =T)
library(colorspace,quietly =T)
library(readxl,quietly =T)
library(lme4,quietly =T)
library(dplyr,quietly =T)
library(data.table,quietly =T)
library(ggpubr,quietly =T)
library(stringr,quietly =T)
# colors
col_wdw <- 25;
col_palette <- diverge_hcl(col_wdw + 1, h = c(260, 0), c = 100, l = c(50, 90), power = 1)
col_palette1 <-sequential_hcl(n = 11, h = 10, c = c(65, 100, NA), l = c(30, 90), power = 0.8 )

# pheno_f<-"pheno.txt"
# siteinfo<-"site_info.xlsx"
# metdatafile<-"MetData_allENVs.csv"
# metvar<-NULL


pheno_f<-args$pheno
siteinfo<-args$siteinfo
metdatafile<-args$metdata
metvar<-args$metvar


if (is.null(metvar)){
  metvar <- "All"
}

# pheno_f<-"pheno_orign_gjw.txt"
# siteinfo<-"site_info.xlsx"
# metdatafile<-"MetData_allENVs.csv"
# metvar <- "Dim.1,Dim.2,Dim.3,Dim.4"

trait.data<-read.delim(file = pheno_f,header = T)
#trait.data <- na.omit(trait.data)
env_meta_info1<-read_excel(siteinfo)
env_meta_info<-env_meta_info1[,c("Envs","LAT","LON","SeedingDate","Year","Site")]
names(env_meta_info)<-c("env_code","lat","lon","PlantingDate","TrialYear","env_note")
head(env_meta_info)
print(paste0("Met file content: ",length(unique(env_meta_info$env_code))," Envs."))

## those move to Rscript 'prepareMet4search.R'
# met.data<-read.csv(file = metdatafile,header = T,encoding = "UTF-8")
# met.data$date<-as.Date(paste(met.data$year,met.data$month,met.data$day,sep = "-"))
# met.data$Envs<-paste0(substr(met.data$year,start = 3,stop = 4),met.data$Site)
# met.data1<-merge(met.data,env_meta_info[,c(3,5)])
# met.data1$date<-as.Date(met.data1$date)
# met.data1$SeedingDate<-as.Date(met.data1$SeedingDate)
# met.data1<-met.data1[met.data1$date>=met.data1$SeedingDate,]
# met.data2<-met.data1[,c("Envs","date","radn","dh","tmean","tmax","tmin","pcp","rhmean","rhmin","wsmean","wsmax")]
# met.data2$GDD<-apply(select(met.data2,tmean), 1, FUN = function(x){max(0,x-10)})
# met.data2$DTR<-abs(met.data2$tmax-met.data2$tmin)
# met.data2$PRDTR<-met.data2$pcp/met.data2$DTR
# met.data2$PTT<-met.data2$GDD*met.data2$dh
# # met.data2$PTQ<-met.data2$radn/met.data2$GDD
# met.data2[met.data2$GDD!=0,"PTQ"]<-met.data2[met.data2$GDD!=0,"radn"]/met.data2[met.data2$GDD!=0,"GDD"]
# met.data2[is.na(met.data2$PTQ),"PTQ"]<-met.data2[is.na(met.data2$PTQ),"radn"]/0.1
# met.data2<-met.data2 %>% setorder(Envs,date)
# write.csv(file="MetData_allENVs.csv", x=met.data2,row.names = F)

met.data2 <- read.csv(file = metdatafile,header = T)

vars <- names(trait.data)[c(-1,-2)]

if (metvar=="All") {
  Paras <- names(met.data2)[c(-1,-2)]
}else{
  metvars<-unlist(str_split(string = metvar,pattern = ","))
  # as.character(str_split(string = metvar,pattern = ","))
  Paras <- names(met.data2)[c(-1,-2)]
  Paras <- intersect(metvars, Paras)
}

# trait<-"BW"
source("CERIS-JGRA_functions.r")
# source("Z:/yanhl/MyScript/R/functions/CERIS-JGRA_functions.r")
searching_daps<-150
trait<-"BW"
print(">>>Now doing Exhaustive search ...")
for (trait in vars) {
  print(paste0("Exhaustive search for ", trait," ..."))
  lInd <- which(colnames(trait.data) == 'Sample')
  eInd <- which(colnames(trait.data) == 'ENVs')
  tInd <- which(colnames(trait.data) == trait)
  
  exp_trait_dir <- paste("./", trait,  '/',  sep = '')
  
  if (!dir.exists(exp_trait_dir))  { dir.create(exp_trait_dir, recursive= T)}
  
  exp_trait <- trait.data[,c(lInd, eInd, tInd)]
  colnames(exp_trait) <- c("line_code","env_code","Yobs")
  exp_trait$Yobs <- as.numeric(exp_trait$Yobs)
  exp_trait <- aggregate(Yobs ~  line_code + env_code, data = exp_trait, mean,na.rm=T) 
  ## To make sure only one phenotype record per each line in each environment
  exp_trait <- exp_trait[!is.na(exp_trait$Yobs),]
  
  all_env_codes <- unique(exp_trait$env_code)
  
  env_cols <- rainbow_hcl(length(all_env_codes), c = 80, l = 60, start = 0, end = 300, fixup = TRUE, alpha = 0.75)
  
  line_codes <- unique(exp_trait$line_code)
  exp_trait_m <- exp_trait
  exp_trait_m$env_code <- as.factor(exp_trait_m$env_code)
  exp_trait_m$line_code <- as.factor(exp_trait_m$line_code)
  
  lm_ <- lmer(Yobs ~ env_code + (1|line_code), data = exp_trait_m)
  mlm_col <- c(fixef(lm_)[1], fixef(lm_)[1] + fixef(lm_)[-1]) ## BLUE for environnment
  
  # ranef(lm_)$env_code+fixef(lm_)[1]
  
  env_mean_trait_0 <- data.frame(
    env_code = as.vector(unique(exp_trait_m$env_code)), 
    meanY = mlm_col)
  
  env_mean_trait <- env_mean_trait_0[order(env_mean_trait_0$meanY),]
  
  ### pairwise correlations among environments; trait distribution across environments;
  ### two figures and the correspondent output files will be saved in the trait directory;
  try(Pairwise_trait_env_distribution_plot(exp_trait, exp_trait_dir, trait, all_env_codes, env_meta_info))
  
  ##### searching the critical window having the strongest correlation with environmental mean
  ##### the window can be adjusted based on biological meaning
  ##### 'FT_9Envs_PTTPTR_0LOO_cor.txt' stores all correlations from all the tested windows and environmental parameters;
  ##### 'MaxR_FTgdd_9Envs_0LOO.png' is the visualization 
  pop_cor_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Envs_PTTPTR_', 0, 'LOO_cor.txt', sep = '')
  # if (file.exists(pop_cor_file)) {
  #   file.remove(pop_cor_file)
  # }
  met.data2<-as.data.frame(met.data2)
  Exhaustive_search2(env_mean_trait, met.data2, searching_daps, exp_trait_dir, trait.data[[trait]], trait, 1, searching_daps, searching_daps, 0, Paras, pop_cor_file)#; searching_daps, searching_daps);
}

#

# CERIS_line_plot ---------------------------------------------------------

print(">>>Now ploting CERIS_line_plot ...")

if (!dir.exists("./CERIS_line_plot"))  { dir.create("./CERIS_line_plot", recursive= T)}

# trait<-"LP"
for (trait in vars){
  print(paste0("ploting for ", trait," ..."))
  exp_trait_dir <- paste("./", trait,  '/',  sep = ''); if (!dir.exists(exp_trait_dir))  { dir.create(exp_trait_dir, recursive= T)};
  pop_cor_file <- list.files(path = exp_trait_dir,
                             pattern = 'LOO_cor.txt',
                             full.names = T)
  pop_cor<-fread(file = pop_cor_file,header = T)
  windows_num <- nrow(pop_cor)

  # grep(pattern = "R_|P_",x = names(pop_cor),value = T)  
  
  r.data<-pop_cor%>%select(c("Day_x","Day_y",grep(pattern = "R_",x = names(pop_cor),value = T),"R2")) %>% setDT()
  p.data<-pop_cor%>%select(c("Day_x","Day_y",grep(pattern = "P_",x = names(pop_cor),value = T),"P_Value")) %>% setDT()
 
  r.data.melt<-melt(data = r.data,
                    id.vars = c("Day_x","Day_y"),
                    variable.name = "Traits",value.name = "R")
# unique(r.data.melt$Traits)

  r.data.melt$Traits<-str_split_fixed(string = r.data.melt$Traits,pattern = "_",n = 2)[,2]
  r.data.melt[r.data.melt$Traits=="","Traits"]<-"Mix"
  
  p.data.melt<-melt(data = p.data,
                    id.vars = c("Day_x","Day_y"),
                    variable.name = "Traits",value.name = "P")
  p.data.melt$Traits<-str_split_fixed(string = p.data.melt$Traits,pattern = "_",n = 2)[,2]
  p.data.melt[p.data.melt$Traits=="Value","Traits"]<-"Mix"
  
  rp.data.melt1<-merge(r.data.melt,p.data.melt)
  # unique(rp.data.melt$Traits)
  #rp.data.melt<-rp.data.melt1[rp.data.melt1$Traits %in% c("DL","PTR","tmax","tmean","tmin","radn","dh","pcp","rhmean","wsmean","DTR","GDD","PRDTR","PTQ","PTT","Dim.1","Dim.2","Dim.3","Dim.4","Dim.5","Mix"),]
  rp.data.melt<-rp.data.melt1[rp.data.melt1$Traits %in% c("DL","DTR","GDD","PRDTR","PTR","PTQ","PTT","dh","pcp","radn","rhmean","tmax","tmean","tmin","wsmean"),]
  #rp.data.melt$Traits<-factor(x = rp.data.melt$Traits,levels = c("DL","PTR","tmax","tmean","tmin","radn","dh","pcp","rhmean","wsmean","DTR","GDD","PRDTR","PTQ","PTT","Dim.1","Dim.2","Dim.3","Dim.4","Dim.5","Mix"))
  rp.data.melt$Traits<-factor(x = rp.data.melt$Traits,levels = c("DL","DTR","GDD","PRDTR","PTR","PTQ","PTT","dh","pcp","radn","rhmean","tmax","tmean","tmin","wsmean"))
  # rp.data.melt<-melt(data = rp.data,
  #                    id.vars = c("Day_x","Day_y"),
  #                    variable.name = "Traits")
  # 
  # rp.data.melt$Y<-str_split_fixed(string = rp.data.melt$Traits,pattern = "_",n = 2)[,1]
  # rp.data.melt$Traits<-str_split_fixed(string = rp.data.melt$Traits,pattern = "_",n = 2)[,2]
  # rp.data.cast<-dcast(data = rp.data.melt,formula = Day_x+Day_y+Traits~Y)
  # 
  p.select<-rp.data.melt[,.SD[which.max(P)],by=c("Day_x","Traits")]
  r.select<-rp.data.melt[,.SD[which.max(abs(R))],by=c("Day_x","Traits")]
  
  p.select1<-rp.data.melt[,.SD[which.max(P)],by=c("Traits")]
  write.csv(x = p.select1,file = paste0("./CERIS_line_plot/","select_p_min_" ,trait,".csv"))

  # Bonferroni correction at the a = 0.05 level with 10440 tests
  p_bjust1<--log10(0.05/windows_num)
  p_bjust<--log10(0.05)
  p.select1<-p.select1[p.select1$P>p_bjust,]
  p.select2<-p.select1[p.select1$P>p_bjust1,]
  if(nrow(p.select1)==0){
    p.select1[1,]<-NA
  }else{
    write.csv(x = p.select1,
              file = paste0("./CERIS_line_plot/","select_p_bjust_" ,trait,".csv"))
    
  };

  if(nrow(p.select2)==0){
    p.select2[1,]<-NA
  }else{
    write.csv(x = p.select2,
              file = paste0("./CERIS_line_plot/","select_p_bjust2_" ,trait,".csv"))
  };
  
  p.p<-ggplot(data = p.select,mapping = aes(x=Day_x,y=P))+
    facet_grid(~Traits)+
    scale_x_continuous(name = "Days after seeding",breaks = seq(0,150,40))+
    scale_y_continuous(name = "-log(P)")+
    geom_line(color="#3e90e2")+
    geom_hline(yintercept = p_bjust,color="darkgrey")+
    geom_hline(yintercept = p_bjust1,color="darkgrey",linetype=2)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.spacing.x = unit(0,"cm"),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "grey"),
      # axis.line.x = element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      # axis.line.y = element_line(color = "black")
    )
  p.r<-ggplot(data = p.select,mapping = aes(x=Day_x,y=R))+
    facet_grid(~Traits)+
    scale_x_continuous(name = "Days after seeding",breaks = seq(0,150,40))+
    scale_y_continuous(name = "r",breaks = seq(-1,1,0.5))+
    geom_line(color="#3e90e2")+
    geom_hline(yintercept = 0,color="black")+
        # geom_text(data = p.select1,
    #           aes(x=Day_x,y=P,label=paste0(Day_x,"-",Day_y)))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.spacing.x = unit(0,"cm"),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "grey"),
      strip.text = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  if(nrow(p.select1)!=0){
    p.p<-p.p+geom_segment(data = p.select1,
                          aes(x=Day_x,y=0,xend=Day_x,yend=P),
                          color="red",linetype="dashed")+
      geom_text(data = p.select1,
                aes(x=Day_x,y=P,label=paste0(Day_x,"-",Day_y)))
    
    p.r<-p.r+geom_segment(data = p.select1,
                          aes(x=Day_x,y=0,xend=Day_x,yend=R),
                          color="red",linetype="dashed")
  }
  
  p<-ggarrange(p.p,p.r,ncol = 1,nrow = 2,align="v")
  ggsave(filename = paste0("./CERIS_line_plot/", trait,"_log_p.png"),
         plot = p,width = 14,height = 4)
  write.csv(x = p.select,file = paste0("./CERIS_line_plot/", trait,"_p_select.csv"),
            row.names = F)
  
  r.p<-ggplot(data = r.select,mapping = aes(x=Day_x,y=P))+
    facet_wrap(facets = vars(Traits),nrow = 1)+
    scale_x_continuous(name = "Days after seeding",breaks = seq(0,150,40))+
    scale_y_continuous(name = "-log(P)")+
    geom_line(color="#3e90e2")+
    geom_hline(yintercept = p_bjust,color="darkgrey")+
    geom_hline(yintercept = p_bjust1,color="darkgrey",linetype=2)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.spacing.x = unit(0,"cm"),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "grey"),
      strip.text = element_blank(),
      axis.line = element_line(color = "black")
    )
  r.r<-ggplot(data = r.select,mapping = aes(x=Day_x,y=R))+
    facet_wrap(facets = vars(Traits),nrow = 1)+
    # scale_x_continuous(name = "Days after seeding",breaks = seq(0,150,40))+
    scale_y_continuous(name = "r",breaks = seq(-1,1,0.5))+
    geom_line(color="#3e90e2")+
    geom_hline(yintercept = 0,color="black")+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.spacing.x = unit(0,"cm"),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "grey"),
      # axis.line.x = element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      # axis.line.y = element_line(color = "black")
    )
  p<-ggarrange(r.r,r.p,ncol = 1,nrow = 2,align="v")
  ggsave(filename = paste0("./CERIS_line_plot/", trait,"_r_select.png"),
         plot = p,width = 14,height = 4)
  write.csv(x = r.select,file = paste0("./CERIS_line_plot/", trait,"_r_select.csv"),row.names = F)
}


# Sample predict ---------------------------------------------------------------
trait="BW"

print(">>>Now doing Sample predict ...")
for (trait in vars) {
  print(paste0("Sample predict for ", trait," ..."))
  lInd <- which(colnames(trait.data) == 'Sample')
  eInd <- which(colnames(trait.data) == 'ENVs')
  tInd <- which(colnames(trait.data) == trait)
  exp_trait_dir <- paste("./", trait,  '/',  sep = '')
  if (!dir.exists(exp_trait_dir))  { dir.create(exp_trait_dir, recursive= T)}
  exp_trait <- trait.data%>% as.data.frame() %>% select(all_of(c(lInd, eInd, tInd)))
  colnames(exp_trait) <- c("line_code","env_code","Yobs")
  exp_trait <- aggregate(Yobs ~  line_code + env_code, data = exp_trait, mean,na.rm=T) ## To make sure only one phenotype record per each line in each environment
  exp_trait <- exp_trait[!is.na(exp_trait$Yobs),]
  
   all_env_codes <- unique(exp_trait$env_code)
  
  env_cols <- rainbow_hcl(length(all_env_codes), c = 80, l = 60, start = 0, end = 300, fixup = TRUE, alpha = 0.75)
  ### remove outlier environments, such as one with high missing rate
  # if (trait == 'FT') {
  #   env_outliers <- c('03NY06', '08MO07');
  #   exp_trait <- exp_trait[!(exp_trait$env_code %in% env_outliers),];
  # }
  
  line_codes <- unique(exp_trait$line_code)
  
  env_mean_method <- 'mlm'; ### 'mlm', 'emm'. Modifying this based on the method to calculate environmental mean
  exp_trait_m <- exp_trait
  exp_trait_m$env_code <- as.factor(exp_trait_m$env_code)
  exp_trait_m$line_code <- as.factor(exp_trait_m$line_code)
  
  if (env_mean_method == 'ari') {
    env_mean_trait_0 <- na.omit(aggregate(x = exp_trait$Yobs, by = list(env_code = exp_trait$env_code), mean, na.rm = T));
    colnames(env_mean_trait_0)[2] <- 'meanY'    
  } else if (env_mean_method == 'mlm') {
    lm_ <- lmer(Yobs ~ env_code + (1|line_code), data = exp_trait_m)
    mlm_col <- c(fixef(lm_)[1], fixef(lm_)[1] + fixef(lm_)[-1]) ## BLUE for environnment
    env_mean_trait_0 <- data.frame(env_code = as.vector(unique(exp_trait_m$env_code)), meanY = mlm_col)
  } else if (env_mean_method == 'emm') {
    env_n <- length(as.vector(unique(exp_trait_m$env_code)))
    trait_ori_lm1 <- lm(Yobs ~ env_code + line_code, data = exp_trait_m)
    trait_ori.pred1 <- matrix(predict(ref_grid(trait_ori_lm1)), nrow = env_n)
    emm_col <- apply(trait_ori.pred1, 1, mean) ### marginal mean for environments
    env_mean_trait_0 <- data.frame(env_code = as.vector(unique(exp_trait_m$env_code)), meanY = emm_col)
  }
  
  env_mean_trait <- env_mean_trait_0[order(env_mean_trait_0$meanY),]
  
  ### Change the following three parameters for the window and environmental parameter with the strongest correlation
  selet.p<-read.csv(file = paste0("./CERIS_line_plot/","select_p_min_" ,trait,".csv"),
                    row.names = 1)
  # vv="tmin"
  vvs<-intersect(unique(selet.p$Traits),c("DL","DTR","GDD","PRDTR","PTR","PTQ","PTT","dh","pcp","radn","rhmean","tmax","tmean","tmin","wsmean","Dim.1","Dim.2","Dim.3","Dim.4","Dim.5"))
  for (vv in vvs) {
    
    maxR_dap1 <- selet.p[selet.p$Traits==vv,"Day_x"]
    maxR_dap2 <- selet.p[selet.p$Traits==vv,"Day_y"]
    kPara_Name <- vv
    
    PTT_PTR<-met.data2
    names(PTT_PTR)[1]<-"env_code"
    PTT_PTR_ind <-  which(colnames(PTT_PTR) == kPara_Name); 
    
    #### Visualization of the relationships between environmental mean and environmental parameters from the selected window.  
    Plot_Trait_mean_envParas(env_mean_trait, PTT_PTR, maxR_dap1, maxR_dap2, trait, exp_trait_dir, env_cols, Paras)
    
    #### Output intercept and slope estimation for each line based on environmental mean and environmental parameter
    #Slope_Intercept(maxR_dap1, maxR_dap2, env_mean_trait, PTT_PTR, exp_trait, line_codes, exp_trait_dir);
    
    Slope_Intercept_para(maxR_dap1, maxR_dap2, env_mean_trait, PTT_PTR, PTT_PTR_ind, exp_trait, line_codes, exp_trait_dir);
    #### LOOCV function for 1 -> 2
    obs_prd_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Env_LOO_by_Lines_', kPara_Name, 'D', maxR_dap1, '_', maxR_dap2, '.txt', sep = '');
    LOO_pdf_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Env_LOO_by_Lines_', kPara_Name, 'D', maxR_dap1, '_', maxR_dap2, '.png', sep = '');
    p=1
    # if (!file.exists(obs_prd_file)) { 
      prdM <- LOOCV(maxR_dap1, maxR_dap2, env_mean_trait, PTT_PTR, PTT_PTR_ind, exp_trait, obs_prd_file, p)
    # }
    Plot_prediction_result(obs_prd_file, all_env_codes, prdM, kPara_Name, LOO_pdf_file,env_cols);
    
  }  
}

print(">>>All have done <<<")