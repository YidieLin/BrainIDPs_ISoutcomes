################################################################################
##############################MR-Forward##################################
##################################By:Yidie######################################
################################2025年3月18日##################################

library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(mr.raps)
library(MendelianRandomization)
library(RadialMR)

setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\相关材料")

idpinfo <- fread("IDPsinfo.csv", header = T, stringsAsFactors = F)
idpinfo <- idpinfo %>% mutate(exposure = gsub("IDP.","IDP",`IDP ID`))

########################Part1 主分析##########################
###########（1）寻找并剔除混杂SNP############
###匹配暴露SNP的rsID
library(tidyverse)
library(data.table)
setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\data\\analysisData\\proxy")
mydat <- fread("data_harmonised_proxy.tsv", sep = "\t", header = T, stringsAsFactors = F)

snplist <- mydat %>% distinct(SNP, .keep_all = TRUE) %>%
  select(chr.exposure, pos.exposure, SNP, exposure)
expnames <- unique(snplist$exposure)
data_fill <- data.frame()
for(expname in expnames){
  tryCatch({
    print(paste0("-----开始处理",expname,"-------"))
    
    data <- snplist %>% filter(exposure == expname)
    print(paste0("-----待处理SNP数量",nrow(data),"-------"))
    
    fill_c <- get_rsID_from_chrpos(data,
                                   col_chr = "chr.exposure",
                                   col_pos = "pos.exposure",
                                   build = "37",
                                   pop = "EUR",
                                   database = "dbsnp"
    )
    
    data_fill <- rbind(data_fill, fill_c)
    
    print(paste("已合并:", Sys.time ()))
  }, error = function(e) {
    message(paste0("处理 ", expname, " 时发生错误：", e$message))
  })
}

fwrite(data_fill, paste0("G:/奕蝶/03文章/BrainPoststroke/分析空间/results/proxy/","data_fill.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

#####处理暴露SNP信息
library(tidyverse)
library(data.table)
setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\data\\analysisData\\proxy")
mydat <- fread("G:/奕蝶/03文章/BrainPoststroke/分析空间/results/proxy/data_fill.tsv", sep = "\t", header = T, stringsAsFactors = F)

library(MungeSumstats)
colnames(mydat) <- c("CHR","BP", "SNP", "exposure")

mydat2 <- liftover(mydat, 
                   chain_source = "ucsc",
                   ref_genome = "hg19", convert_ref_genome = "hg38")
colnames(mydat2) <- c("SNP", "CHR_hg38", "POS_hg38", "exposure", "Imputed_Gen_build") 

mydat3 <- mydat %>% 
  left_join(mydat2, by = c("SNP","exposure")) %>%
  mutate(ID_hg19 = paste0(CHR, ":", BP),
         ID_hg38 = paste0(CHR_hg38, ":", POS_hg38))

mydat3$ID_hg38[mydat3$SNP == "rs7544145"] <- "1:150166501"
mydat3$ID_hg38[mydat3$SNP == "rs62010099"] <- "15:82327998"
mydat3$ID_hg38[mydat3$SNP == "rs77515140"] <- "17:36899951"

#####处理混杂因素SNP数据
setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\data\\confounderSNPs")
files <- list.files(pattern = "_export.tsv")

data <- data.frame()
for(file in files){
  print(paste0("----开始合并",file,"的数据-------"))
  
  dat <- fread(file, header = T, stringsAsFactors = F)
  
  data <- rbind(data, dat)
  
}

data2 <- data %>%
  select("riskAllele","pValue","traitName","efoTraits","accessionId","locations") %>%
  distinct(locations, .keep_all = TRUE)

######剔除暴露中的混杂SNP
filterd <- data2 %>% filter(locations %in% mydat3$ID_hg38)
colnames(filterd) <- c("riskAllele", "pValue", "traitName", "efoTraits", "accessionId", 
                       "ID_hg38")
SNP_confoundlist <- filterd %>%
  left_join(mydat3, by = "ID_hg38")
fwrite(SNP_confoundlist, 
       paste0("G:/奕蝶/03文章/BrainPoststroke/分析空间/results/","SNP_confoundlists.tsv"), 
       sep = "\t", row.names = FALSE, quote = FALSE)

setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\data\\analysisData\\proxy")
e_o_dat <- fread("data_harmonised_proxy_radial.tsv", sep = "\t", header = T, stringsAsFactors = F)

e_o_dat2 <- e_o_dat %>% 
  filter(!SNP %in% SNP_confoundlist$ID_hg19)
#check <- e_o_dat %>% 
#  filter(SNP %in% SNP_confoundlist$ID_hg19)
fwrite(e_o_dat2, "data_harmonised_proxy_radial_noconfounder.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

###########（2）剔除混杂SNP后进行分析############
data_mr_res <- data.frame()
data_het <- data.frame()
data_pleio <- data.frame()
data_Fstat_I2 <- data.frame()
##
setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\data\\analysisData\\proxy")
mydat <- fread("data_harmonised_proxy_radial_noconfounder.tsv", sep = "\t", header = T, stringsAsFactors = F)
expnames <- unique(mydat$exposure)
#outnames <- unique(mydat$outcome)
outnames <- c("mRS_012vs3456_adj", "mRS_012vs3456_noadj", "nihss")
for (expname in expnames) {
  for (outname in outnames) {
    print(paste0("-----------BEGIN: ", expname," + ",outname, "-------------"))
    ##读取数据
    e_o_dat <- mydat %>% 
      filter(exposure == expname) %>%
      filter(outcome == outname)
    
    if(nrow(e_o_dat %>% filter(mr_keep == TRUE)) == 0){
      print("无可用SNPs")
    } else if(nrow(e_o_dat %>% filter(mr_keep == TRUE)) == 1) {
      print("SNP只有1个，采用Wald Ratio")
      mr_main <- mr(e_o_dat, method_list = "mr_wald_ratio")
      mr_res_main <- generate_odds_ratios(mr_main) %>% 
        mutate(`OR(95%CI)` = sprintf("%.2f (%.2f to %.2f)", or, or_lci95, or_uci95),
               P = ifelse(pval < 0.001, sprintf("%.3e", pval), sprintf("%.3f", pval))) %>%
        select(outcome, exposure, method, nsnp, b, se, or, or_lci95, or_uci95, `OR(95%CI)`, P)
      mr_res <-mr_res_main
      data_mr_res <- rbind(data_mr_res, mr_res)
    } else {
      ###MR分析
      mr_main <- mr(e_o_dat, method_list = c('mr_ivw', 'mr_weighted_median', 'mr_egger_regression'))
      mr_res_main <- generate_odds_ratios(mr_main) %>% 
        mutate(`OR(95%CI)` = sprintf("%.2f (%.2f to %.2f)", or, or_lci95, or_uci95),
               P = ifelse(pval < 0.001, sprintf("%.3e", pval), sprintf("%.3f", pval))) %>%
        select(outcome, exposure, method, nsnp, b, se, or, or_lci95, or_uci95, `OR(95%CI)`, P)
      
      ###MR_raps
      nsnp <- nrow(e_o_dat %>% filter(mr_keep == TRUE))
      if(nsnp < 3){
        print("SNPs数量小于3，不进行MR-Raps分析")
        mr_res <- mr_res_main
        print("MR分析结果汇总完毕")
      } else{
        mr_sen_raps <- data.frame(mr.raps(e_o_dat$beta.exposure, e_o_dat$beta.outcome,
                                          e_o_dat$se.exposure, e_o_dat$se.outcome))
        mr_sen_raps <- mr_sen_raps %>%
          mutate(exposure = expname,
                 outcome = outname,
                 method = "Robust adjusted profile score",
                 nsnp = nsnp) %>%
          dplyr::select("outcome","exposure","method","nsnp","beta.hat", "beta.se", "beta.p.value")
        
        colnames(mr_sen_raps) <- c("outcome", "exposure", "method", "nsnp", "b", "se", "pval")
        mr_res_raps <- generate_odds_ratios(mr_sen_raps) %>%
          mutate(`OR(95%CI)` = sprintf("%.2f (%.2f to %.2f)", or, or_lci95, or_uci95),
                 P = ifelse(pval < 0.001, sprintf("%.3e", pval), sprintf("%.3f", pval))) %>%
          select(outcome, exposure, method, nsnp, b, se, or, or_lci95, or_uci95, `OR(95%CI)`, P)
        
        print("MR-Raps分析完毕")
        
        mr_res <- rbind(mr_res_main, mr_res_raps)
        print("MR分析结果汇总完毕")
      }
      
      ###Calculate heterogeneity ###snps数量< 2的不会计算
      if(nsnp < 3){
        print("SNP < 3, 不计算Heterogeneity")
        het <- data.frame()
      } else {
        het <- mr_heterogeneity(e_o_dat) %>%  
          filter(method == "Inverse variance weighted") %>%   
          select("outcome", "exposure", "method", "Q", "Q_df", "Q_pval")
      }
      ###Calculate Egger intercept and Pvalue ###snps数量< 3的不会计算
      if(nsnp < 3){
        print("SNP < 3, 不计算Egger intercept")
        pleio <- data.frame()
      }else{
        pleio <- mr_pleiotropy_test(e_o_dat) %>%
          select("outcome", "exposure", "egger_intercept", "se", "pval")
      }
      
      ##Calculate Fstats and I2
      e_o_dat2 <- e_o_dat %>% filter(mr_keep == TRUE)
      object_in <- mr_input(bx = e_o_dat2$beta.exposure, bxse = e_o_dat2$se.exposure, 
                            by = e_o_dat2$beta.outcome, byse = e_o_dat2$se.outcome)
      if(nsnp < 3){
        print("SNPs < 3, 不计算I2")
        I2 <- NA
      } else {
        I2 <- mr_egger(object_in)@I.sq
      }
      nsnp <- mr_divw(object_in)@SNPs
      Fstat <- round(mr_divw(object_in)@Condition,2)
      Fstat_I2 <- data.frame(outname, expname, nsnp, Fstat, I2)
      colnames(Fstat_I2) <- c("outcome", "exposure", "nsnp", "F-statistics", "I2")
      
      ####Summarize Results
      data_mr_res <- rbind(data_mr_res, mr_res)
      data_het <- rbind(data_het, het)
      data_pleio <- rbind(data_pleio, pleio)
      data_Fstat_I2 <- rbind(data_Fstat_I2, Fstat_I2)
    }
    
    print(paste0("---------------------END--------------------"))
  }
}   
fwrite(data_mr_res, file = "../../../results/proxy/BaseRadialMR/Results_MRmain_noconfounder.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_het, file = "../../../results/proxy/BaseRadialMR/Results_heterogeneity_noconfounder.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_pleio, file = "../../../results/proxy/BaseRadialMR/Results_pleiotropy_noconfounder.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_Fstat_I2, file = "../../../results/proxy/BaseRadialMR/Results_Fstat_I2_noconfounder.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


##################（3）补充：匹配IDP info##################
library(data.table)
library(tidyverse)
setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\相关材料")
idpinfo <- fread("IDPsinfo.csv", header = T, stringsAsFactors = F)


setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\results\\proxy\\baseRadialMR")
data <- fread("Results_MRmain_noconfounder.tsv", sep = "\t", header = T)

data2 <- data %>%
  left_join(idpinfo, by = "IDP ID") %>%
  mutate(b_lci95 = b-1.96*se,
         b_uci95 = b+1.96*se) %>%
  mutate(`b(95%CI)` = sprintf("%.3f(%.3f, %.3f)", b, b_lci95, b_uci95)) %>%
  select("outcome", "exposure", "method", "nsnp", "b(95%CI)", "OR(95%CI)", "P", "IDP ID", "IDP name.x",
         "Tissue/Region.x") %>%
  mutate(Estimate = ifelse(outcome == "nihss",`b(95%CI)`, `OR(95%CI)`)) %>%
  filter(outcome != "mRS_012vs3456_noadj")

colnames(data2) <- c("outcome", "exposure", "method", "nsnp", "b(95%CI)", "OR(95%CI)", "P", "IDP ID", "IDP name",
                      "Tissue/Region","Estimate")
res_ivw <- data2 %>%
  filter(method == "Inverse variance weighted") %>%
  select("outcome", "IDP ID", "IDP name","Tissue/Region","nsnp","Estimate", "P")

res_median <- data2 %>%
  filter(method == "Weighted median")%>%
  select("outcome", "IDP ID","nsnp", "Estimate", "P")

res_egger <- data2 %>%
  filter(method == "MR Egger")%>%
  select("outcome", "IDP ID","nsnp", "Estimate", "P")

res_raps <- data2 %>%
  filter(method == "Robust adjusted profile score")%>%
  select("outcome", "IDP ID","nsnp", "Estimate", "P")

res_wald <- data2 %>%
  filter(method == "Wald ratio")%>%
  select("outcome", "IDP ID","nsnp", "Estimate", "P")


res_all <- res_ivw %>%
  left_join(res_median, by = c("outcome", "IDP ID")) %>%
  left_join(res_egger, by = c("outcome", "IDP ID")) %>%
  left_join(res_raps, by = c("outcome", "IDP ID"))

colnames(res_all) <- c("outcome", "IDP ID", "IDP name", "Tissue/Region", "nsnp", 
                       "IVW", "PIVW", "nsnp.y", "Median", "Pmedian", "nsnp.x.x", 
                       "egger", "Pegger", "nsnp.y.y", "raps", "Praps")


setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\results")
fwrite(res_all, paste0("Results_MRmain_forward_tidy.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)




###########################Part2. 敏感性分析####################################
#########################01 不剔除混杂SNP########################
data_mr_res <- data.frame()
data_het <- data.frame()
data_pleio <- data.frame()
data_Fstat_I2 <- data.frame()
##
setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\data\\analysisData\\proxy")
mydat <- fread("data_harmonised_proxy_radial.tsv", sep = "\t", header = T, stringsAsFactors = F)
expnames <- unique(mydat$exposure)
#outnames <- unique(mydat$outcome)
outnames <- c("mRS_012vs3456_adj", "mRS_012vs3456_noadj", "nihss")
for (expname in expnames) {
  for (outname in outnames) {
    print(paste0("-----------BEGIN: ", expname," + ",outname, "-------------"))
    ##读取数据
    e_o_dat <- mydat %>% 
      filter(exposure == expname) %>%
      filter(outcome == outname)
    
    if(nrow(e_o_dat %>% filter(mr_keep == TRUE)) == 0){
      print("无可用SNPs")
    } else if(nrow(e_o_dat %>% filter(mr_keep == TRUE)) == 1) {
      print("SNP只有1个，采用Wald Ratio")
      mr_main <- mr(e_o_dat, method_list = "mr_wald_ratio")
      mr_res_main <- generate_odds_ratios(mr_main) %>% 
        mutate(`OR(95%CI)` = sprintf("%.2f (%.2f to %.2f)", or, or_lci95, or_uci95),
               P = ifelse(pval < 0.001, sprintf("%.3e", pval), sprintf("%.3f", pval))) %>%
        select(outcome, exposure, method, nsnp, b, se, or, or_lci95, or_uci95, `OR(95%CI)`, P)
      mr_res <-mr_res_main
      data_mr_res <- rbind(data_mr_res, mr_res)
    } else {
      ###MR分析
      mr_main <- mr(e_o_dat, method_list = c('mr_ivw', 'mr_weighted_median', 'mr_egger_regression'))
      mr_res_main <- generate_odds_ratios(mr_main) %>% 
        mutate(`OR(95%CI)` = sprintf("%.2f (%.2f to %.2f)", or, or_lci95, or_uci95),
               P = ifelse(pval < 0.001, sprintf("%.3e", pval), sprintf("%.3f", pval))) %>%
        select(outcome, exposure, method, nsnp, b, se, or, or_lci95, or_uci95, `OR(95%CI)`, P)
      
      ###MR_raps
      nsnp <- nrow(e_o_dat %>% filter(mr_keep == TRUE))
      if(nsnp < 3){
        print("SNPs数量小于3，不进行MR-Raps分析")
        mr_res <- mr_res_main
        print("MR分析结果汇总完毕")
      } else{
        mr_sen_raps <- data.frame(mr.raps(e_o_dat$beta.exposure, e_o_dat$beta.outcome,
                                          e_o_dat$se.exposure, e_o_dat$se.outcome))
        mr_sen_raps <- mr_sen_raps %>%
          mutate(exposure = expname,
                 outcome = outname,
                 method = "Robust adjusted profile score",
                 nsnp = nsnp) %>%
          dplyr::select("outcome","exposure","method","nsnp","beta.hat", "beta.se", "beta.p.value")
        
        colnames(mr_sen_raps) <- c("outcome", "exposure", "method", "nsnp", "b", "se", "pval")
        mr_res_raps <- generate_odds_ratios(mr_sen_raps) %>%
          mutate(`OR(95%CI)` = sprintf("%.2f (%.2f to %.2f)", or, or_lci95, or_uci95),
                 P = ifelse(pval < 0.001, sprintf("%.3e", pval), sprintf("%.3f", pval))) %>%
          select(outcome, exposure, method, nsnp, b, se, or, or_lci95, or_uci95, `OR(95%CI)`, P)
        
        print("MR-Raps分析完毕")
        
        mr_res <- rbind(mr_res_main, mr_res_raps)
        print("MR分析结果汇总完毕")
      }
      
      ###Calculate heterogeneity ###snps数量< 2的不会计算
      if(nsnp < 3){
        print("SNP < 3, 不计算Heterogeneity")
        het <- data.frame()
      } else {
        het <- mr_heterogeneity(e_o_dat) %>%  
          filter(method == "Inverse variance weighted") %>%   
          select("outcome", "exposure", "method", "Q", "Q_df", "Q_pval")
      }
      ###Calculate Egger intercept and Pvalue ###snps数量< 3的不会计算
      if(nsnp < 3){
        print("SNP < 3, 不计算Egger intercept")
        pleio <- data.frame()
      }else{
        pleio <- mr_pleiotropy_test(e_o_dat) %>%
          select("outcome", "exposure", "egger_intercept", "se", "pval")
      }
      
      ##Calculate Fstats and I2
      e_o_dat2 <- e_o_dat %>% filter(mr_keep == TRUE)
      object_in <- mr_input(bx = e_o_dat2$beta.exposure, bxse = e_o_dat2$se.exposure, 
                            by = e_o_dat2$beta.outcome, byse = e_o_dat2$se.outcome)
      if(nsnp < 3){
        print("SNPs < 3, 不计算I2")
        I2 <- NA
      } else {
        I2 <- mr_egger(object_in)@I.sq
      }
      nsnp <- mr_divw(object_in)@SNPs
      Fstat <- round(mr_divw(object_in)@Condition,2)
      Fstat_I2 <- data.frame(outname, expname, nsnp, Fstat, I2)
      colnames(Fstat_I2) <- c("outcome", "exposure", "nsnp", "F-statistics", "I2")
      
      ####Summarize Results
      data_mr_res <- rbind(data_mr_res, mr_res)
      data_het <- rbind(data_het, het)
      data_pleio <- rbind(data_pleio, pleio)
      data_Fstat_I2 <- rbind(data_Fstat_I2, Fstat_I2)
    }
    
    print(paste0("---------------------END--------------------"))
  }
}

data_mr_res2 <- data_mr_res %>%
  left_join(idpinfo, by = "exposure") %>%
  select("outcome", "exposure", "method", "nsnp", "b", "se", "or", "or_lci95", 
         "or_uci95", "OR(95%CI)", "P", "IDP name", "IDP ID", "Tissue/Region")

data_het2 <- data_het %>%
  left_join(idpinfo, by = "exposure") %>%
  select("outcome", "exposure", "method", "Q", "Q_df", "Q_pval", "IDP name", 
         "IDP ID", "Tissue/Region")

data_pleio2 <- data_pleio %>%
  left_join(idpinfo, by = "exposure") %>%
  select("outcome", "exposure", "egger_intercept", "se", "pval", "IDP name", 
         "IDP ID", "Tissue/Region")

data_Fstat_I22 <- data_Fstat_I2 %>%
  left_join(idpinfo, by = "exposure") %>%
  select("outcome", "exposure", "nsnp", "F-statistics", "I2", "IDP name", 
         "IDP ID", "Tissue/Region")

fwrite(data_mr_res2, file = "../../../results/proxy/BaseRadialMR/Results_MRmain.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_het2, file = "../../../results/proxy/BaseRadialMR/Results_heterogeneity.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_pleio2, file = "../../../results/proxy/BaseRadialMR/Results_pleiotropy.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_Fstat_I22, file = "../../../results/proxy/BaseRadialMR/Results_Fstat_I2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



############补充：整合剔除混杂SNP前后的结果################
library(tidyverse)
library(data.table)
setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\results")
data1 <- fread("./proxy/BaseRadialMR/Results_MRmain.tsv", sep = "\t", header = T)
data2 <- fread("Results_MRmain_forward_tidy.tsv", sep = "\t", header = T)
head(data1)

data3 <- data1 %>%
  filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>%
  mutate(b_lci95 = b-1.96*se,
         b_uci95 = b+1.96*se) %>%
  mutate(`b(95%CI)` = sprintf("%.3f(%.3f, %.3f)", b, b_lci95, b_uci95)) %>%
  select("outcome", "exposure", "method", "nsnp", "b(95%CI)", "OR(95%CI)", "P", "IDP ID", "IDP name",
         "Tissue/Region") %>%
  mutate(Estimate = ifelse(outcome == "nihss",`b(95%CI)`, `OR(95%CI)`)) %>%
  filter(outcome != "mRS_012vs3456_noadj") %>%
  select("outcome", "IDP ID","nsnp", "Estimate", "P")
  
data4 <- data2 %>%
  left_join(data3, by = c("outcome", "IDP ID"))


setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\results")
fwrite(data4, paste0("Results_MRmain_forward_tidy_V2.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)



###########################02 不剔除Radial检测的outliers########################
data_mr_res <- data.frame()
data_het <- data.frame()
data_pleio <- data.frame()
data_Fstat_I2 <- data.frame()

#SNP_confoundlist <- fread(paste0("G:/奕蝶/03文章/BrainPoststroke/分析空间/results/","SNP_confoundlists.tsv"), 
#                           sep = "\t", header = T, stringsAsFactors = F)
#setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\data\\analysisData\\proxy")
#mydat <- fread("data_harmonised_proxy.tsv", sep = "\t", header = T, stringsAsFactors = F)

#e_o_dat2 <- mydat %>% 
#  filter(!SNP %in% SNP_confoundlist$ID_hg19)
#check <- e_o_dat %>% 
#  filter(SNP %in% SNP_confoundlist$ID_hg19)
#fwrite(e_o_dat2, "data_harmonised_proxy_noconfounder_noradial.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\data\\analysisData\\proxy")
mydat <- fread("data_harmonised_proxy_noconfounder_noradial.tsv", sep = "\t", header = T, stringsAsFactors = F)

expnames <- unique(mydat$exposure)
#outnames <- unique(mydat$outcome)
outnames <- c("mRS_012vs3456_adj", "mRS_012vs3456_noadj", "nihss")
for (expname in expnames) {
  for (outname in outnames) {
    print(paste0("-----------BEGIN: ", expname," + ",outname, "-------------"))
    ##读取数据
    e_o_dat <- mydat %>% 
      filter(exposure == expname) %>%
      filter(outcome == outname)
    
    if(nrow(e_o_dat %>% filter(mr_keep == TRUE)) == 0){
      print("无可用SNPs")
    } else if(nrow(e_o_dat %>% filter(mr_keep == TRUE)) == 1) {
      mr_main <- mr(e_o_dat, method_list = "mr_wald_ratio")
      mr_res_main <- generate_odds_ratios(mr_main) %>% 
        mutate(`OR(95%CI)` = sprintf("%.2f (%.2f to %.2f)", or, or_lci95, or_uci95),
               P = ifelse(pval < 0.001, sprintf("%.3e", pval), sprintf("%.3f", pval))) %>%
        select(outcome, exposure, method, nsnp, b, se, or, or_lci95, or_uci95, `OR(95%CI)`, P)
    } else {
      ###MR分析
      mr_main <- mr(e_o_dat, method_list = c('mr_ivw', 'mr_weighted_median', 'mr_egger_regression'))
      mr_res_main <- generate_odds_ratios(mr_main) %>% 
        mutate(`OR(95%CI)` = sprintf("%.2f (%.2f to %.2f)", or, or_lci95, or_uci95),
               P = ifelse(pval < 0.001, sprintf("%.3e", pval), sprintf("%.3f", pval))) %>%
        select(outcome, exposure, method, nsnp, b, se, or, or_lci95, or_uci95, `OR(95%CI)`, P)
      
      ###MR_raps
      nsnp <- nrow(e_o_dat %>% filter(mr_keep == TRUE))
      if(nsnp < 3){
        print("SNPs数量小于3，不进行MR-Raps分析")
        mr_res <- mr_res_main
        print("MR分析结果汇总完毕")
      } else{
        mr_sen_raps <- data.frame(mr.raps(e_o_dat$beta.exposure, e_o_dat$beta.outcome,
                                          e_o_dat$se.exposure, e_o_dat$se.outcome))
        mr_sen_raps <- mr_sen_raps %>%
          mutate(exposure = expname,
                 outcome = outname,
                 method = "Robust adjusted profile score",
                 nsnp = nsnp) %>%
          dplyr::select("outcome","exposure","method","nsnp","beta.hat", "beta.se", "beta.p.value")
        
        colnames(mr_sen_raps) <- c("outcome", "exposure", "method", "nsnp", "b", "se", "pval")
        mr_res_raps <- generate_odds_ratios(mr_sen_raps) %>%
          mutate(`OR(95%CI)` = sprintf("%.2f (%.2f to %.2f)", or, or_lci95, or_uci95),
                 P = ifelse(pval < 0.001, sprintf("%.3e", pval), sprintf("%.3f", pval))) %>%
          select(outcome, exposure, method, nsnp, b, se, or, or_lci95, or_uci95, `OR(95%CI)`, P)
        
        print("MR-Raps分析完毕")
        
        mr_res <- rbind(mr_res_main, mr_res_raps)
        print("MR分析结果汇总完毕")
      }
      
      ###Calculate heterogeneity ###snps数量< 2的不会计算
      if(nsnp < 3){
        print("SNP < 3, 不计算Heterogeneity")
        het <- data.frame()
      } else {
        het <- mr_heterogeneity(e_o_dat) %>%  
          filter(method == "Inverse variance weighted") %>%   
          select("outcome", "exposure", "method", "Q", "Q_df", "Q_pval")
      }
      ###Calculate Egger intercept and Pvalue ###snps数量< 3的不会计算
      if(nsnp < 3){
        print("SNP < 3, 不计算Egger intercept")
        pleio <- data.frame()
      }else{
        pleio <- mr_pleiotropy_test(e_o_dat) %>%
          select("outcome", "exposure", "egger_intercept", "se", "pval")
      }
      
      ##Calculate Fstats and I2
      e_o_dat2 <- e_o_dat %>% filter(mr_keep == TRUE)
      object_in <- mr_input(bx = e_o_dat2$beta.exposure, bxse = e_o_dat2$se.exposure, 
                            by = e_o_dat2$beta.outcome, byse = e_o_dat2$se.outcome)
      if(nsnp < 3){
        print("SNPs < 3, 不计算I2")
        I2 <- NA
      } else {
        I2 <- mr_egger(object_in)@I.sq
      }
      nsnp <- mr_divw(object_in)@SNPs
      Fstat <- round(mr_divw(object_in)@Condition,2)
      Fstat_I2 <- data.frame(outname, expname, nsnp, Fstat, I2)
      colnames(Fstat_I2) <- c("outcome", "exposure", "nsnp", "F-statistics", "I2")
      
      ####Summarize Results
      data_mr_res <- rbind(data_mr_res, mr_res)
      data_het <- rbind(data_het, het)
      data_pleio <- rbind(data_pleio, pleio)
      data_Fstat_I2 <- rbind(data_Fstat_I2, Fstat_I2)
    }
    
    print(paste0("---------------------END--------------------"))
  }
}   
fwrite(data_mr_res, file = "../../../results/proxy/sensitivity/noBaseRadialMR/Results_MRmain.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_het, file = "../../../results/proxy/sensitivity/noBaseRadialMR/Results_heterogeneity.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_pleio, file = "../../../results/proxy/sensitivity/noBaseRadialMR/Results_pleiotropy.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_Fstat_I2, file = "../../../results/proxy/sensitivity/noBaseRadialMR/Results_Fstat_I2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



########################03 MR PRESSO分析############################
library(MRPRESSO)
data_mrpresso <- data.frame()
data_distort <- data.frame()
data_global<- data.frame()
#暴露列表和结局列表
setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\data\\analysisData\\proxy")
mydat <- fread("data_harmonised_proxy_radial_noconfounder.tsv", sep = "\t", header = T, stringsAsFactors = F)

expnames <- unique(mydat$exposure)
#outnames <- unique(mydat$outcome)
outnames <- c("mRS_012vs3456_adj", "nihss")
for (expname in expnames) {
  for (outname in outnames) {
    print(paste0("-----------BEGIN: ", expname," + ",outname, "-------------"))
    #step 1: prepare the data
    e_o_dat <- mydat %>% 
      filter(exposure == expname) %>%
      filter(outcome == outname) %>%
      filter(mr_keep == TRUE)
    
    ##(1)采用PRESSO检测
    if (nrow(e_o_dat) > 5){
      
      #initial
      outlier_index <- "NA"
      dat_presso <- as.data.frame(e_o_dat)
      
      #while loop
      mrpresso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                            SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                            OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                            data = dat_presso, 
                            NbDistribution = 2000, 
                            SignifThreshold = 0.05, 
                            seed = 20130813)
      
      if (is.null(mrpresso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)) {
        print(" No outlier were identified")
      } else {
        outlier_index <- row.names(mrpresso$`MR-PRESSO results`$`Outlier Test`[mrpresso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, ])
        out_distort <- as.data.frame(rbind(mrpresso$`MR-PRESSO results`$`Distortion Test`))
        outlier_index
        data_distort <- rbind(data_distort, out_distort)
      }
      
      out_global <- as.data.frame(rbind(mrpresso$`MR-PRESSO results`$`Global Test`)) %>%
        mutate(outcome = outname,
               exposure = expname) %>%
        select("outcome","exposure","RSSobs","Pvalue")
      
      out_mrpressomain <- as.data.frame(mrpresso$`Main MR results`) %>%
        mutate(outcome = outname,
               exposure = expname) %>%
        select("outcome", "exposure", "MR Analysis", "Causal Estimate", "Sd", "T-stat", "P-value")
      
      #汇总数据
      data_mrpresso <- rbind(data_mrpresso, out_mrpressomain)
      data_global <- rbind(data_global, out_global)
    }
    print(paste("-------------------------------END：分割线--------------------------------------"))
    print(paste("END时间为：", Sys.time ()))
    ##
    gc()
  }
}

fwrite(data_mrpresso, file = "../../../results/proxy/sensitivity/Results_mrpresso_main.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_global, file = "../../../results/proxy/sensitivity/Results_mrpresso_global.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_distort, file = "../../../results/proxy/sensitivity/Results_mrpresso_distort.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


