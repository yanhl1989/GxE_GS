# G × E integrated genomic selection model Framework
This project is to construct genomic selection models that incorporate genotype-environment interactions.
Please refer to our paper for details.

## Pheno analysis
####./1.data/phenotype/phenoAnalysis.R
To calculate BLUP and ANOVA, phenotype files are required, with the phenotype format as follows: 
ENVs	Sample	Trait1	Trait2	Trait3	Trait4	Trait5	Trait6
ENV1	Sample1	5.77	31.02	5.41	29.55	38.71	10.9
ENV1	Sample2	4.8	31.03	4.24	34.5	37.79	11.6
...
ENVn	Samplen	5.01	30.3	4.17	32.2	43.5	10.6


## Plasticity analysis
####./2.plasticity_FWR/sample_plasticity.r
Phenotypic plasticity for analyzing sample traits requires phenotype files, with the phenotype format as follows: 
ENVs	Sample	Trait1	Trait2	Trait3	Trait4	Trait5	Trait6
ENV1	Sample1	5.77	31.02	5.41	29.55	38.71	10.9
ENV1	Sample2	4.8	31.03	4.24	34.5	37.79	11.6
...
ENVn	Samplen	5.01	30.3	4.17	32.2	43.5	10.6


## Key environmental factors search
####./3.ENV_search/ENV_search.r
Used for searching key environmental factors and their window periods. We need to prepare phenotype data, meteorological factor data, and station data.
The phenotype data is as follows: 
ENVs	Sample	Trait1	Trait2	Trait3	Trait4	Trait5	Trait6...
ENV1	Sample1	5.77	31.02	5.41	29.55	38.71	10.9...
ENV1	Sample2	4.8	31.03	4.24	34.5	37.79	11.6...
...
ENVn	Samplen	5.01	30.3	4.17	32.2	43.5	10.6...

The meteorological factor data is as follows: 
Envs,date,factor1,factor2...
ENV1,date1,18.85,24...
ENV1,date2,18.5,18.1...
...
ENVn,daten,20.64,15.7...

The station data is as follows: 
Year	Site	Envs	Area	SeedingDate	LON	LAT
Year1	Site1	ENV1	Northwest	SeedingDate1	80.27	41.17
Year2	Site1	ENV2	Northwest	SeedingDate2	81.29	40.54
...
Yearn	Siten	ENVn	Northwest	SeedingDate2	81.29	40.54


## Core markers select
####./4.cropgbm/p01.getRunCropgbm.sh
Using Cropgbm to screen core markers. We need to prepare phenotype data and genotype data.
The phenotype data is as follows: 
Sample,Trait
Sample1,5.20318117341928
Sample2,4.8
...
Samplen,5.01

The genotype data is as follows: 
IID,Marker1,Marker2...
Sample1,0,2...
Sample2,2,1...
...
Samplen,0,0...


## Model construction and validation
####./5.GS_ENVmodel/p01.predata.sh
Used for preparing G × E prediction mixture matrix. Key environmental factor information needs to be prepared. 
The content of all_dev.list is as follows: 
Trait1	../3.ENV_search/Trait1/Trait1_envMeanPara_day1_day2.txt	factor1
Trait1	../3.ENV_search/Trait1/Trait1_envMeanPara_day1_day2.txt	factor2
...
Traitn	../3.ENV_search/Traitn/Traitn_envMeanPara_day1_day2.txt	factorn

####./5.GS_ENVmodel/extract_ENV.R
Used to extract the mean of environmental factors and assist in the use of p01.pretata.sh.

####./5.GS_ENVmodel/predata.R
Adjust the phenotype file format and assist in the use of p01.pretata.sh.

####./5.GS_ENVmodel/manhattan_blue.R OR manhattan_plasticity.R OR manhattan_env.R
Used for drawing Manhattan maps and assist in the use of p01.pretata.sh.

####./5.GS_ENVmodel/select_compare_random.R
Obtain the optimal core markers and assist in the use of p01.pretata.sh.

####./5.GS_ENVmodel/pregeno.R
Used for preparing G × E prediction mixture matrix and assist in the use of p01.pretata.sh.

####./5.GS_ENVmodel/CV_Serial.R
Model architecture and cross validation.
