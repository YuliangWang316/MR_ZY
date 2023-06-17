library(readxl)
library(tidyverse)
library(data.table)
library(TwoSampleMR)
#library(colochelpR)
library(tidyverse)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

sumpleio <- c()
sumres <-c()

## Step 1 文件导入
setwd("/mnt/ndisk1/Student/zhouyi/data/otherdat/NAFLD")
Eur_guts <- read_excel("Eur_guts.xlsx")
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
#available.SNPs()
exposure  = fs::dir_ls("/mnt/ndisk1/Student/zhouyi/data/otherdat/Germangut/2_test",recurse = TRUE,glob="*.tsv")
outcome = fs::dir_ls("/mnt/ndisk1/Student/zhouyi/data/otherdat/NAFLD",recurse = TRUE,glob="*.tsv")
#outcome5 = read.table("/mnt/ndisk1/Student/zhouyi/data/otherdat/NAFLD/NamjouB_31311600_NAFLD.txt",sep = " ",header =T)
#outcome3 = read.table(outcome[3],header = T)
outcome1 = read.table(outcome[1],header = T)
#outcome2 = read.table(outcome[2],header = T)
#outcome4 = read.table(outcome[4],header = T)

#BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
#BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
n = length(exposure)
n
## 处理outcome数据
#NAFLDname 
#outcome1$phenotype <- NAFLDname
#outcome_2 = outcome_1 %>% filter(p_value<1e-5)

for(i in 1:5){
    print(i)
      for(t in 1:n){
        ## Step 2 暴露与结局数据处理
        print(t)
        exposure  = fs::dir_ls(paste0(paste0("/mnt/ndisk1/Student/zhouyi/data/otherdat/Germangut/",i),"_test"),recurse = TRUE,glob="*.tsv")
        exposure_1 <- read.table(exposure[t],header = T)
  #outcome_1 <- read.table(exposure[4],header = T)
  
        Gutname = str_extract(exposure[t],"(?<=test/).*(?=_)")
        print(Gutname)
        NAFLDname = str_extract(outcome[1],"(?<=LD/).*(?=_)")
  
        Gutrow = Eur_guts %>% filter(`Study accession` == Gutname)
  
        colnames(exposure_1)[1:2] = c("CHR","BP")
  
        data = exposure_1 %>% filter(p_value<1e-5)
        data = data[complete.cases(data[, 2]), ]
        for (i in unique(data$CHR)){
              my_pos <- data$BP[data$CHR == i]
              chr_snps <- snpsBySeqname(snps, as.character(i))
              idx <- match(my_pos,pos(chr_snps))
              rsids <- mcols(chr_snps)$RefSNP_id[idx]
              data$rsid[data$CHR == i] <- rsids
    
        }
  
  # fwrite(data,".csv")
  
            exp_dat <- format_data(
    data,
    type = "exposure",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col ="effect_allele",
    other_allele_col = "other_allele",
    pval_col = "p_value"
            )
            print("yes")
  #a = convert_loc_to_rs(data[,1:2], dbSNP)
  
            exp_dat$id.exposure = Gutrow$`Study accession`
            exp_dat$exposurename = Gutrow$`Reported trait`
  ## 去除连锁不平衡
            exp_dat <-clump_data(exp_dat,clump_r2=0.001,clump_kb=10000)
            setwd(paste0(paste0("/mnt/ndisk1/Student/zhouyi/data/otherdat/Germangut/",i),"_testdeal")) 
            fwrite(exp_dat,paste0(Gutname,".csv"))


  }

}  