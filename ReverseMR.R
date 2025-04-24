################################################################################
##############################MR-Reverse##################################
##################################By:Yidie######################################
################################2025.03.20##################################

library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(mr.raps)
library(MendelianRandomization)
library(RadialMR)


##Note: Reverse分析中暴露数据里没有混杂SNP
#########################Part1. 主分析：radial########################
###
data_mr_res <- data.frame()
data_het <- data.frame()
data_pleio <- data.frame()
data_Fstat_I2 <- data.frame()
##
setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\reverse\\data\\analysisData")
mydat <- fread("data_harmonised_radial.tsv", sep = "\t", header = T, stringsAsFactors = F)
expnames <- unique(mydat$exposure)
outnames <- unique(mydat$outcome)

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

fwrite(data_mr_res, file = "../../results/Results_MRmain_reverse.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_het, file = "../../results/Results_heterogeneity_reverse.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_pleio, file = "../../results//Results_pleiotropy_reverse.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_Fstat_I2, file = "../../results/Results_Fstat_I2_reverse.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

###########################Part2. 敏感性分析####################################
###########################01 不剔除Radial检测的outliers########################
data_mr_res <- data.frame()
data_het <- data.frame()
data_pleio <- data.frame()
data_Fstat_I2 <- data.frame()
##
setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\reverse\\data\\analysisData")
mydat <- fread("data_harmonized.tsv", sep = "\t", header = T, stringsAsFactors = F)
expnames <- unique(mydat$exposure)
outnames <- unique(mydat$outcome)

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

fwrite(data_mr_res, file = "../../results/sensitivity/noBaseRadial/Results_MRmain_reverse.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_het, file = "../../results/sensitivity/noBaseRadial/Results_heterogeneity_reverse.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_pleio, file = "../../results/sensitivity/noBaseRadial/Results_pleiotropy_reverse.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_Fstat_I2, file = "../../results/sensitivity/noBaseRadial/Results_Fstat_I2_reverse.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



########################02 Steiger filtering and 方向检验##########################
setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\reverse\\data\\analysisData")
mydat <- fread("data_harmonised_radial.tsv", sep = "\t", header = T, stringsAsFactors = F)

#######Steiger direction tests
data_direc <- directionality_test(mydat) %>%
  select("outcome", "exposure", "snp_r2.exposure", "snp_r2.outcome", "correct_causal_direction", "steiger_pval")
fwrite(data_direc, file = "../../results/sensitivity/Results_steigerDir.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

########Steiger filtering and MR analysis
data_steiger <- steiger_filtering(mydat)

if(nrow(data_steiger %>% filter(steiger_dir == FALSE)) > 0 ){
  filters_dat <- data_steiger %>% filter(steiger_dir == FALSE)
  filters <- filters_dat$rsID
  snp_string <- paste(filters, collapse = ", ")
  
  mydat2 <- data_steiger %>% filter(steiger_dir == TRUE)
  res_before <- generate_odds_ratios(mr(mydat, method_list = "mr_ivw")) %>%
    mutate(`Beta (95%CI)` = sprintf("%.3f (%.3f to %.3f)", b, lo_ci, up_ci),
           P = ifelse(pval < 0.001, sprintf("%.3e", pval), sprintf("%.3f", pval))) %>%
    select("outcome", "exposure","nsnp", "Beta (95%CI)", "P")
  
  res_after <- generate_odds_ratios(mr(mydat2, method_list = "mr_ivw")) %>%
    mutate(`Beta (95%CI)` = sprintf("%.3f (%.3f to %.3f)", b, lo_ci, up_ci),
           P = ifelse(pval < 0.001, sprintf("%.3e", pval), sprintf("%.3f", pval))) %>%
    select("outcome", "exposure","nsnp", "Beta (95%CI)", "P")
  
  steiger_res <- res_before %>%
    left_join(res_after, by = c("outcome", "exposure"))
} else {
  print("No SNPs needed to be filtered.")
}

##Note：No SNPs needed to be filtered.


########################02 MR PRESSO分析############################
library(MRPRESSO)
data_mrpresso <- data.frame()
data_distort <- data.frame()
data_global<- data.frame()
#暴露列表和结局列表
setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\reverse\\data\\analysisData")
mydat <- fread("data_harmonised_radial.tsv", sep = "\t", header = T, stringsAsFactors = F)
expnames <- unique(mydat$exposure)
outnames <- unique(mydat$outcome)

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

fwrite(data_mrpresso, file = "../../results/sensitivity/Results_mrpresso_main.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_global, file = "../../results/sensitivity/Results_mrpresso_global.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(data_distort, file = "../../results/sensitivity/Results_mrpresso_distort.tsv", sep = "\t", row.names = FALSE, quote = FALSE)








##################补充：匹配IDP info##################
library(data.table)
library(tidyverse)
setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\相关材料")
idpinfo <- fread("IDPsinfo.csv", header = T, stringsAsFactors = F)


setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\reverse\\results")
data <- fread("Results_MRmain_reverse.tsv", sep = "\t", header = T)

data2 <- data %>%
  mutate(`IDP ID` = outcome) %>%
  left_join(idpinfo, by = "IDP ID") %>%
  mutate(b_lci95 = b-1.96*se,
         b_uci95 = b+1.96*se) %>%
  mutate(`b(95%CI)` = sprintf("%.3f(%.3f, %.3f)", b, b_lci95, b_uci95)) %>%
  select("outcome", "exposure", "method", "nsnp", "b(95%CI)", "OR(95%CI)", "P", "IDP ID", "IDP name",
         "Tissue/Region") %>%
  mutate(Estimate = `b(95%CI)`) %>%
  filter(exposure != "mRS_012vs3456_noadj")


res_ivw <- data2 %>%
  filter(method == "Inverse variance weighted") %>%
  select("exposure", "IDP ID", "IDP name","Tissue/Region","nsnp","Estimate", "P")

res_median <- data2 %>%
  filter(method == "Weighted median")%>%
  select("exposure", "IDP ID","nsnp", "Estimate", "P")

res_egger <- data2 %>%
  filter(method == "MR Egger")%>%
  select("exposure", "IDP ID","nsnp", "Estimate", "P")

res_raps <- data2 %>%
  filter(method == "Robust adjusted profile score")%>%
  select("exposure", "IDP ID","nsnp", "Estimate", "P")

res_wald <- data2 %>%
  filter(method == "Wald ratio")%>%
  select("exposure", "IDP ID","nsnp", "Estimate", "P")


res_all <- res_ivw %>%
  left_join(res_median, by = c("exposure", "IDP ID")) %>%
  left_join(res_egger, by = c("exposure", "IDP ID")) %>%
  left_join(res_raps, by = c("exposure", "IDP ID"))

colnames(res_all) <- c("exposure", "IDP ID", "IDP name", "Tissue/Region", "nsnp", 
                       "IVW", "PIVW", "nsnp.y", "Median", "Pmedian", "nsnp.x.x", 
                       "egger", "Pegger", "nsnp.y.y", "raps", "Praps")


setwd("G:\\奕蝶\\03文章\\BrainPoststroke\\分析空间\\reverse\\results")
fwrite(res_all, paste0("Results_MRmain_reverse_tidy.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
