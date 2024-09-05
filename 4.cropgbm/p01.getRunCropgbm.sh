#!/bin/bash
#$1的文件格式参考/public/home/zzuaga04/hanxv/yancaoNP600/5.GSmodel/1.data/1.genotype/genotype_data_modified.raw
#$2参考/public/home/zzuaga04/hanxv/yancaoNP600/5.GSmodel/1.data/2.pheotype/trait.list
genoRaw=$1
traitFile=$2
cul=`pwd`

if [ ! "$traitFile" ] ; then
    echo "not find trait.lib file!"
    exit
fi
if [ ! -e ./result ] ; then
    mkdir ./result
fi

if [ -e ./cropgbm_run.list ] ; then
    date1=`date +%m_%d_%H:%M`
    #mv cropgbm_run.list cropgbm_run.list_${date1} 
    rm cropgbm_run.list
fi

cat $traitFile | while read line
do 

    line_basename=$(basename "$line")
    if [[ $line_basename == *".trait" ]]; then
        trait="${line_basename%%.*}"
    else
        trait=$(basename $(dirname $(readlink -f "$line")))_${line_basename%%.txt}
    fi
    if [ ! -e ./result/$trait ] ; then
        mkdir ./result/$trait
    fi
    if [ -e $cul/result/$trait/${trait}_snpselect.sh ] ; then
        rm $cul/result/$trait/${trait}_snpselect.sh
    fi
    if [[ $line_basename == *".trait" ]]; then
        echo  "
cd $cul/result/$trait
if [ ! -e $cul/result/$trait/gbm_result ] ; then
   mkdir $cul/result/$trait/gbm_result
fi
if [ ! -e $cul/result/$trait/gbm_result/engine ] ; then
   mkdir $cul/result/$trait/gbm_result/engine
fi
if [ ! -e ./gbm_result/engine/genotype_data_cleaned.feature ] ; then
    sed 's/\t/,/g' "$line" > modified.trait
    cropgbm -o gbm_result/ -e -t -sf --bygain-boxplot --traingeno "$genoRaw" --trainphe modified.trait --min-gain 0.05 --max-colorbar 0.6 --cv-times 6
    sync
fi
">>$cul/result/$trait/${trait}_snpselect.sh
    else
       echo  "
cd $cul/result/$trait
if [ ! -e $cul/result/$trait/Intcp_${trait} ] ; then
   mkdir $cul/result/$trait/Intcp_${trait}
fi 
if [ ! -e $cul/result/$trait/Intcp_${trait}/gbm_result ] ; then
   mkdir $cul/result/$trait/Intcp_${trait}/gbm_result
fi 
if [ ! -e $cul/result/$trait/Intcp_${trait}/gbm_result/engine ] ; then
   mkdir $cul/result/$trait/Intcp_${trait}/gbm_result/engine
fi 
cd $cul/result/$trait/Intcp_${trait}
awk -v OFS='\t' '{print \$1, \$2}' "$line" > "Intcp_${trait}.trait"
sed 's/\t/,/g' "Intcp_${trait}.trait" > Intcp_modified.trait
cropgbm -o gbm_result/ -e -t -sf --bygain-boxplot --traingeno "$genoRaw" --trainphe Intcp_modified.trait --min-gain 0.05 --max-colorbar 0.6 --cv-times 6
sync
cd $cul/result/$trait
if [ ! -e $cul/result/$trait/Slope_${trait} ] ; then
   mkdir $cul/result/$trait/Slope_${trait}
fi 
if [ ! -e $cul/result/$trait/Slope_${trait}/gbm_result ] ; then
   mkdir $cul/result/$trait/Slope_${trait}/gbm_result
fi 
if [ ! -e $cul/result/$trait/Slope_${trait}/gbm_result/engine ] ; then
   mkdir $cul/result/$trait/Slope_${trait}/gbm_result/engine
fi 
cd $cul/result/$trait/Slope_${trait}
awk -v OFS='\t' '{print \$1, \$3}' "$line" > "Slope_${trait}.trait"
sed 's/\t/,/g' "Slope_${trait}.trait" > Slope_modified.trait
cropgbm -o gbm_result/ -e -t -sf --bygain-boxplot --traingeno "$genoRaw" --trainphe Slope_modified.trait --min-gain 0.05 --max-colorbar 0.6 --cv-times 6
sync

">>$cul/result/$trait/${trait}_snpselect.sh 
fi
 echo $cul/result/$trait/${trait}_snpselect.sh >> cropgbm_run.list
done

