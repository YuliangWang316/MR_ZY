library(readxl)
library(tidyverse)
library(data.table)
library(TwoSampleMR)
#library(colochelpR)
library(tidyverse)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(MRPRESSO)
outcome = fs::dir_ls("/mnt/ndisk1/Student/zhouyi/data/otherdat/NAFLD",recurse = TRUE,glob="*.tsv")
outcome1 = read.table(outcome[1],header = T)
outcome2 = read.table(outcome[2],header = T)

## Step 1 定义函数
Newharharmonise_data = function(exp_dat,out_dat){
    harmonise_data(exposure_dat = exp_dat,outcome_dat = out_dat,action=2)
}
computmr_presso = function(x){
  mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = x, NbDistribution = 1000,  SignifThreshold = 0.05)

}

search_outcome = function(exp,out){
  extract_outcome_data(
  snps=exp$SNP,
  outcomes=out,
  proxies = FALSE,
  maf_threshold = 0.01,
  access_token = NULL)
}
## step 2 循环输出结果
setwd("/mnt/ndisk1/Student/zhouyi/data/otherdat/NAFLD")
Eur_guts <- read_excel("Eur_guts.xlsx")
for(i in 1:5){
    print(i)
    exp_datas = fs::dir_ls(paste0(paste0("/mnt/ndisk1/Student/zhouyi/data/otherdat/Germangut/",i),"_testdeal"),recurse = TRUE,glob = "*.csv")
    sumres <-c()
    for(t in 1:length(exp_datas)){
        print(exp_datas[t])
        exp_dat <- fread(exp_datas[t],header = T)
  #outcome_1 <- read.table(exposure[4],header = T)
  
        Gutname = str_extract(exp_datas[t],"(?<=testdeal/).*(?=.csv)")
        Gutrow = Eur_guts %>% filter(`Study accession` == Gutname)
        print(Gutname)
        mydata <- harmonise_data(
                exposure_dat=exp_dat,
                outcome_dat=out_dat_1,
                action= 2
            )
        # 绘图
        if(nrow(mydata)==0){
                print("byebye")
            }else{
    
            pleio <- mr_pleiotropy_test(mydata)#多效性检验，检测mr-egger的截距项与0有无统计学差异
            sumpleio <- rbind(sumpleio,pleio)
            sumpleio$exposure = Gutrow$`Reported trait`
            sumpleio$outcome = NAFLDname
            het <- mr_heterogeneity(mydata)#异质性检验，如果Q_pval的值小于0.05，则使用随机效应模型估计mr
    
            res <- mr(mydata)
            res$exposure = Gutrow$`Reported trait`
            res$outcome = NAFLDname
            res <- rbind(res,res[5,]) 
            res$method[6] = "MR-PRESSO"
            mp = computmr_presso(mydata)
            if(length(na.omit(mp$`Main MR results`[,3]))==2){
                    res$b[6] = mp$`Main MR results`$'Causal Estimate'[2]
                    res$pval[6] = mp$`Main MR results`$'P-value'[2]
  
                }else{
                    res$b[6] = mp$`Main MR results`$'Causal Estimate'[1]
                    res$pval[6] = mp$`Main MR results`$'P-value'[1]
            }
            sumres = rbind(sumres,res)
        }
    
        
    }
    setwd("/mnt/ndisk1/Student/zhouyi/data/otherdat/NAFLD")
    filename = paste0(i,"_NAFLD.csv")
    print(filename)
    fwrite(sumres,filename)


}


### NAFLD芬兰数据库
### finn-b-NAFLD
id = c("finn-b-NAFLD")


for(i in 1:5){
    print(i)
    exp_datas = fs::dir_ls(paste0(paste0("/mnt/ndisk1/Student/zhouyi/data/otherdat/Germangut/",i),"_testdeal"),recurse = TRUE,glob = "*.csv")
    sumres <-c()
    for(t in 1:length(exp_datas)){
        print(exp_datas[t])
        exp_dat <- fread(exp_datas[t],header = T)
  #outcome_1 <- read.table(exposure[4],header = T)
        Gutname = str_extract(exp_datas[t],"(?<=testdeal/).*(?=.csv)")
        Gutrow = Eur_guts %>% filter(`Study accession` == Gutname)
        print(Gutname)
        mydata <- Newharharmonise_data(exp_dat,search_outcome(exp_dat,"finn-b-NAFLD"))             
        # 绘图
        if(nrow(mydata)==0){
                print("byebye")
            }else{
    
            #pleio <- mr_pleiotropy_test(mydata)#多效性检验，检测mr-egger的截距项与0有无统计学差异
            #sumpleio <- rbind(sumpleio,pleio)
            #sumpleio$exposure = Gutrow$`Reported trait`
            #sumpleio$outcome = NAFLDname
            #het <- mr_heterogeneity(mydata)#异质性检验，如果Q_pval的值小于0.05，则使用随机效应模型估计mr
    
            res <- mr(mydata)
            res$exposure = Gutrow$`Reported trait`
            
            res <- rbind(res,res[5,]) 
            res$method[6] = "MR-PRESSO"
            mp = computmr_presso(mydata)
            if(length(na.omit(mp$`Main MR results`[,3]))==2){
                    res$b[6] = mp$`Main MR results`$'Causal Estimate'[2]
                    res$pval[6] = mp$`Main MR results`$'P-value'[2]
  
                }else{
                    res$b[6] = mp$`Main MR results`$'Causal Estimate'[1]
                    res$pval[6] = mp$`Main MR results`$'P-value'[1]
            }
            sumres = rbind(sumres,res)
        }
    
        
    }
    setwd("/mnt/ndisk1/Student/zhouyi/data/otherdat/NAFLD")
    filename = paste0(paste0(i,"_NAFLD_"),"Finn.csv")
    print(filename)
    fwrite(sumres,filename)


}


## 芬兰数据库

