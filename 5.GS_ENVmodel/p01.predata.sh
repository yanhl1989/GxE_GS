#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH -p normal
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24


cul=$(pwd)

if [ -e $cul/predata.list ]; then
    rm $cul/predata.list
fi

if [ ! -e $cul/1.data/1.genotype ]; then
     mkdir -p $cul/1.data/1.genotype
fi

if [ ! -e $cul/1.data/2.phenotype ]; then
     mkdir -p $cul/1.data/2.phenotype
fi

if [ ! -e $cul/1.data/3.ENV ]; then
     mkdir -p $cul/1.data/3.ENV
fi

if [ ! -e $cul/1.data/4.geno_ENV ]; then
     mkdir -p $cul/1.data/4.geno_ENV
fi

cd $cul/1.data/2.phenotype

ln -s $cul/../1.data/phenotype/0-153_phenotype.txt

Rscript $cul/predata.R "0-153_phenotype.txt"

cd $cul/1.data/3.ENV

cd $cul
# BW LP SI FL FS FM
for trait in BW LP SI FL FS FM;do
      
      if [ ! -e $cul/1.data/3.ENV/$trait ] ; then
          mkdir -p $cul/1.data/3.ENV/$trait
      fi
    
      cd $cul/1.data/3.ENV/$trait

        if [ -e env.list ] ; then
             date1=`date +%m_%d_%H:%M` 
           rm env.list
        fi

   awk -v trait="$trait" '$1 == trait {print $2 "\t" $3}' $cul/all_env.list > env.list

   Rscript $cul/extract_ENV.R env.list 
    

cd $cul

cd $cul/1.data/1.genotype

ln -s $cul/../1.data/genotype/Geno_map.xlsx

cd $cul 

   if [ ! -e $cul/1.data/1.genotype/$trait ]; then
     mkdir $cul/1.data/1.genotype/$trait
   fi
   cd $cul/1.data/1.genotype/$trait

   # 提取表型选择所得的snp的012信息
  if [ ! -e ${trait}_plasticity_ENV.raw ] ;then 

   if [ ! -e $cul/1.data/1.genotype/$trait/BLUP ]; then
      mkdir $cul/1.data/1.genotype/$trait/BLUP
   fi
   cd $cul/1.data/1.genotype/$trait/BLUP
   awk -F',' '{gsub(/featureid/,"IID",$1); print $1 "," $2}' $cul/../4.cropgbm/result/$trait/gbm_result/engine/genotype_data_origin.feature > ${trait}.list
   Rscript $cul/manhattan_blue.R -M $cul/1.data/1.genotype/Geno_map.xlsx -S ${trait}.list
   Rscript $cul/select_compare_random.R -G $cul/../1.data/genotype/genotype_data_origin.raw -P $cul/../1.data/phenotype/${trait}.trait -C ${trait}.list


   # 提取表型可塑性所得的snp的012信息
   if [ ! -e $cul/1.data/1.genotype/$trait/plasticity/Slope ]; then
      mkdir -p $cul/1.data/1.genotype/$trait/plasticity/Slope
   fi
   if [ ! -e $cul/1.data/1.genotype/$trait/plasticity/Intercept ]; then
      mkdir -p $cul/1.data/1.genotype/$trait/plasticity/Intercept
   fi
   cd $cul/1.data/1.genotype/$trait/plasticity/Slope 
   Rscript $cul/select_compare_random.R -G $cul/../1.data/genotype/genotype_data_origin.raw -P $cul/../2.plasticity_FWR/slope_${trait}.trait -C $cul/../4.cropgbm/result/slope_${trait}/gbm_result/engine/genotype_data_origin.feature
   cd $cul/1.data/1.genotype/$trait/plasticity/Intercept 
   Rscript $cul/select_compare_random.R -G $cul/../1.data/genotype/genotype_data_origin.raw -P $cul/../2.plasticity_FWR/Intercept_${trait}.trait -C $cul/../4.cropgbm/result/Intercept_${trait}/gbm_result/engine/genotype_data_origin.feature
   cd $cul/1.data/1.genotype/$trait/plasticity   
   awk -F',' 'NR>1 { print $1 "," $2 ",Slope" }' $cul/../4.cropgbm/result/slope_${trait}/gbm_result/engine/genotype_data_origin.feature > slope_${trait}.list
   awk -F',' 'NR>1 { print $1 "," $2 ",Intercept" }' $cul/../4.cropgbm/result/Intercept_${trait}/gbm_result/engine/genotype_data_origin.feature > Intercept_${trait}.list
   cat slope_${trait}.list Intercept_${trait}.list > plasticity.list   #含表型snp选择所得的snp
   Rscript $cul/manhattan_plasticity.R -M $cul/1.data/1.genotype/Geno_map.xlsx -S plasticity.list
   

   # 提取环境搜索所得的snp的012信息
   find "$cul/../4.cropgbm/result/" -type d -name "${trait}_Intcp_Slope_*" | while read -r line
   do
     line1=$(basename "$line") #FTdap_Intcp_Slope_dh_D64_70_7envs
     line2=${line1#*Slope_}    #dh_D64_70_7envs
     env=${line2%%_D*}         #dh
     # 提取最优数量的snp， slope
     if [ ! -e $cul/1.data/1.genotype/$trait/ENV/$env/Slope ]; then
        mkdir -p $cul/1.data/1.genotype/$trait/ENV/$env/Slope
     fi  
     cd $cul/1.data/1.genotype/$trait/ENV/$env/Slope
     awk -F'\t' '{ print $1 "," $3 }' $cul/../3.ENV_search/$trait/Intcp_Slope_${line2}.txt > modified.trait
     Rscript $cul/select_compare_random.R -G $cul/../1.data/genotype/genotype_data_origin.raw -P modified.trait -C $cul/../4.cropgbm/result/$line1/Slope_${line1}/gbm_result/engine/genotype_data_origin.feature
     # 提取最优数量的snp， Intercept
     if [ ! -e $cul/1.data/1.genotype/$trait/ENV/$env/Intercept ]; then
        mkdir -p $cul/1.data/1.genotype/$trait/ENV/$env/Intercept
     fi  
     cd $cul/1.data/1.genotype/$trait/ENV/$env/Intercept
     awk -F'\t' ' { print $1 "," $2}' $cul/../3.ENV_search/$trait/Intcp_Slope_${line2}.txt > modified.trait
     Rscript $cul/select_compare_random.R -G $cul/../1.data/genotype/genotype_data_origin.raw -P modified.trait -C $cul/../4.cropgbm/result/$line1/Intcp_${line1}/gbm_result/engine/genotype_data_origin.feature
     awk -F',' 'NR>1 { print $1 "," $2 "," "'${env}_Slope'" }' "$cul/../4.cropgbm/result/$line1/Slope_${line1}/gbm_result/engine/genotype_data_origin.feature" >> "$cul/1.data/1.genotype/$trait/ENV/ENV_${trait}.list"        
     awk -F',' 'NR>1 { print $1 "," $2 "," "'${env}_Intcp'" }' "$cul/../4.cropgbm/result/$line1/Intcp_${line1}/gbm_result/engine/genotype_data_origin.feature" >> "$cul/1.data/1.genotype/$trait/ENV/ENV_${trait}.list"

   done
   cd $cul/1.data/1.genotype/$trait/ENV
   Rscript $cul/manhattan_env.R -M $cul/1.data/1.genotype/Geno_map.xlsx -S ENV_${trait}.list
   cd $cul/1.data/1.genotype/$trait
   cat ./BLUP/BLUP_select_snp.list ./plasticity/Plasticity_select_snp.list ./ENV/ENV_search_select_snp.list > ratio_select.list
   sort -u ratio_select.list -o ratio_select.list
   awk -F',' 'NR==FNR{a[$1]; next} FNR==1 {for (i=1; i<=NF; i++) if ($i in a) col[i]=1} {row_data=$1; for (i=2; i<=NF; i++) if (i in col) row_data=row_data","$i; print row_data}' ratio_select.list $cul/../1.data/genotype/genotype_data_origin.raw | sed 's/^,/,NA,/' > ratio_select.raw
   fi

   if [ ! -e $cul/1.data/4.geno_ENV/ratio_select/$trait ]; then
     mkdir -p $cul/1.data/4.geno_ENV/ratio_select/$trait
   fi

    cd $cul/1.data/4.geno_ENV/ratio_select/$trait

    Rscript $cul/pregeno.R $cul/1.data/2.phenotype/ENVs_Sample.list $cul/1.data/1.genotype/$trait/ratio_select.raw $cul/1.data/3.ENV/$trait/output.txt
    
    cd $cul/1.data/1.genotype/$trait

    echo "IID" > best_num_select.tmp

    find "$cul"/1.data/1.genotype/"$trait" -name "genotype.list" | while read file
    do
       cat best_num_select.tmp "$file" > temp_output.tmp
       mv temp_output.tmp best_num_select.tmp
    done

    mv best_num_select.tmp best_num_select.list

    awk -F',' 'NR==FNR{a[$1]; next} FNR==1 {for (i=1; i<=NF; i++) if ($i in a) col[i]=1} {row_data=$1; for (i=2; i<=NF; i++) if (i in col) row_data=row_data","$i; print row_data}' best_num_select.list $cul/../1.data/genotype/genotype_data_origin.raw | sed 's/^,/,NA,/' > best_num_select.raw

    if [ ! -e $cul/1.data/4.geno_ENV/best_num_select/$trait ]; then
     mkdir -p $cul/1.data/4.geno_ENV/best_num_select/$trait
    fi

    cd $cul/1.data/4.geno_ENV/best_num_select/$trait

    Rscript $cul/pregeno.R $cul/1.data/2.phenotype/ENVs_Sample.list $cul/1.data/1.genotype/$trait/best_num_select.raw $cul/1.data/3.ENV/$trait/output.txt

    cd $cul
  done



