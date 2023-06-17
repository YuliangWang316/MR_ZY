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
exposure  = fs::dir_ls("/mnt/ndisk1/Student/zhouyi/data/otherdat/Germangut/1_test",recurse = TRUE,glob="*.tsv")
outcome = fs::dir_ls("/mnt/ndisk1/Student/zhouyi/data/otherdat/NAFLD",recurse = TRUE,glob="*.tsv")
#outcome5 = read.table("/mnt/ndisk1/Student/zhouyi/data/otherdat/NAFLD/NamjouB_31311600_NAFLD.txt",sep = " ",header =T)
#outcome3 = read.table(outcome[3],header = T)
outcome1 = read.table(outcome[1],header = T)
#outcome2 = read.table(outcome[2],header = T)
#outcome4 = read.table(outcome[4],header = T)

#BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
#BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
n = length(exposure)
## 处理outcome数据
#NAFLDname 
#outcome1$phenotype <- NAFLDname
#outcome_2 = outcome_1 %>% filter(p_value<1e-5)

out_dat_1 <- format_data(
  dat=outcome1,
  type = "outcome",
  # snps = exp_dat$SNP,
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "variant_id",
  beta_col = "lnOR",
  se_col = "standard_error",
  effect_allele_col ="effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  chr_col = "chromosome"
)
for(t in 1:n){
  ## Step 2 暴露与结局数据处理
  print(t)
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
  
  #a = convert_loc_to_rs(data[,1:2], dbSNP)
  
  exp_dat$id.exposure = Gutrow$`Study accession`
  exp_dat$exposurename = Gutrow$`Reported trait`
  ## 去除连锁不平衡
  exp_dat <-clump_data(exp_dat,clump_r2=0.01,clump_kb=5000)
  
  

  
  # 合并
  mydata <- harmonise_data(
    exposure_dat=exp_dat,
    outcome_dat=out_dat_1,
    action= 2
  )
  # 绘图
  if(nrow(mydata)==0){
    print("byebye")
  } else {
    
    pleio <- mr_pleiotropy_test(mydata)#多效性检验，检测mr-egger的截距项与0有无统计学差异
    sumpleio <- rbind(sumpleio,pleio)
    sumpleio$exposure = Gutrow$`Reported trait`
    sumpleio$outcome = NAFLDname
    het <- mr_heterogeneity(mydata)#异质性检验，如果Q_pval的值小于0.05，则使用随机效应模型估计mr
    
    res <- mr(mydata)
    res$exposure = Gutrow$`Reported trait`
    res$outcome = NAFLDname
    sumres <- rbind(sumres,res)
    
    #散点图
    #png('Scatter_Sleep_duration.png',units="in",width=7,height=7,res=600)
    #mr_scatter_plot(res,mydata)#查看各模型估计值与原始估计值的位置
    #dev.off()
    
    #res_single <- mr_singlesnp(mydata)
    #森林图
    #png('Forest_Sleep_duration.png',units="in",width=7,height=7,res=600)
    #mr_forest_plot(res_single)#对比单个模型和单个SNP估计值之间的差异
    #dev.off()
    #漏斗图
    #png('Funnel_Sleep_duration.png',units="in",width=7,height=7,res=600)
    #mr_funnel_plot(res_single)#漏斗图,各模型之间的评估结果
    #dev.off()
  }
}
colnames(sumpleio)[1] = c("Study accession")
t = inner_join(sumpleio,Eur_guts,by="Study accession")
t$exposure = 0
fwrite(t,"1_1_res.csv")
fwrite(sumres,"1_1_res.csv")
