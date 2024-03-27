library(data.table)
library(ggplot2)
library(plyr)
library(tidyr)
library(haven)
library(twopartm)

source("scripts/prs_dict.R")
source("scripts/immune_dict.R")

load("data/merged_all.Rdata")

##---------------------------------------------------------
## CIBERSORTx
##---------------------------------------------------------

for (p_thres in c("0.05","none")){
  for (method in c("abs_cibersortx","rel_cibersortx")){
    for (subtype in c("all","erp","ern","LumA")){
      for(stage in c("all","stage12","stage34")){
        
        prs_base_res_list <- list()
        prs_full_res_list <- list()
        i <- 1
        j <- 1
        
        for (tissue in c("tumor","normal")){
          
          if (p_thres=="0.05" & method=="abs_cibersortx"){
            temp.dat <- abs_cibersort_all_dat[abs_cibersort_all_dat$TissueType==tissue,]
            temp.dat <- temp.dat[temp.dat$P.value < 0.05,]
          } else if (p_thres=="none" & method=="abs_cibersortx"){
            temp.dat <- abs_cibersort_all_dat[abs_cibersort_all_dat$TissueType==tissue,]
          } else if (p_thres=="0.05" & method=="rel_cibersortx"){
            temp.dat <- rel_cibersort_all_dat[rel_cibersort_all_dat$TissueType==tissue,]
            temp.dat <- temp.dat[temp.dat$P.value < 0.05,]
          } else if (p_thres=="none" & method=="rel_cibersortx"){
            temp.dat <- rel_cibersort_all_dat[rel_cibersort_all_dat$TissueType==tissue,]
          }
          
          if (stage=="all"){
            temp.dat <- temp.dat
          } else if (stage=="stage12"){
            temp.dat <- temp.dat[temp.dat$stage<=2,]
          } else if (stage=="stage34") {
            temp.dat <- temp.dat[temp.dat$stage>=3,]
          } 
          
          if (subtype=="all"){
            temp.dat <- temp.dat
          } else if (subtype=="erp"){
            temp.dat <- temp.dat[temp.dat$er==1,]
          } else if (subtype=="ern"){
            temp.dat <- temp.dat[temp.dat$er==0,]
          } else if (subtype %in% c("LumA")){
            temp.dat <- temp.dat[temp.dat$subtype_PAM50_parker==subtype,]
          }
          
          for (pgs_id in pgs_ids){
            
            scaled_prs <- temp.dat[[pgs_id]]
            
            for (cell.type in lm22.selected){
              
              infiltrates <- qnorm((rank(temp.dat[[cell.type]],na.last="keep")-0.5)/sum(!is.na(temp.dat[[cell.type]])))
              mod.family <- gaussian()
              temp.dat$CombinedBatch <- as.factor(temp.dat$CombinedBatch)
              
              if (tissue == "tumor"){
                
                temp.dat$stage_factor <- ifelse(temp.dat$stage %in% c(1,2), "Stage I/II", "Stage III/IV")
                temp.dat$grade <- as.factor(temp.dat$grade)
                temp.dat$er_factor <- ifelse(temp.dat$er==1,"ER+","ER-")
                
                if (pgs_id!="PGS000841"){
                  temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                            onebef_bmi,platform,stage_factor,grade,er_factor,CombinedBatch))
                } else {
                  temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                            platform,stage_factor,grade,er_factor,CombinedBatch))
                }
                
              } else if (tissue == "normal"){
                if (pgs_id!="PGS000841"){
                  temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                            onebef_bmi,platform,CombinedBatch))
                } else{
                  temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                            platform,CombinedBatch))
                }
                
              }
              
              if (tissue=="tumor"){
                if (subtype %in% c("erp","ern")){temp.dat_sub <- subset(temp.dat_sub, select=-er_factor)}
                if (stage!="all") {temp.dat_sub <- subset(temp.dat_sub, select=-stage_factor)}
              }
              
              temp.mod_base <- glm(infiltrates ~ scaled_prs + agedx + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family=mod.family, data=temp.dat)
              temp.mod_full <- glm(infiltrates ~ scaled_prs + ., family=mod.family, data=temp.dat_sub)
              
              prs_base_res_list[[i]] <- c(tissue,pgs_id,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_traits,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_fullname,method,cell.type,p_thres,subtype,stage,dim(temp.dat)[1],"baseline",
                                          as.character(mod.family[1]),formula(temp.mod_base),
                                          summary(temp.mod_base)$coefficients["scaled_prs",],
                                          rep(NA,8),
                                          summary(temp.mod_base)$coefficients["agedx",][c(1,4)])
              prs_full_res_list[[j]] <- c(tissue,pgs_id,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_traits,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_fullname,method,cell.type,p_thres,subtype,stage,dim(temp.dat_sub)[1],"full",
                                          as.character(mod.family[1]),formula(temp.mod_full),
                                          summary(temp.mod_full)$coefficients["scaled_prs",],
                                          rep(NA,8),
                                          summary(temp.mod_full)$coefficients["agedx",][c(1,4)])
              i <- i + 1
              j <- j + 1
              
              if (cell.type %in% lm22.selected_twopart & !(subtype == "ern" & stage == "stage34") & !(subtype == "LumA" & stage == "stage34")){
                
                tpm.dat <- temp.dat
                tpm.dat$infiltrates <- temp.dat[[cell.type]]
                tpm.dat$scaled_prs <- scaled_prs
                
                temp.mod_base_tpm <- tpm(infiltrates ~ scaled_prs + agedx + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=tpm.dat, link_part1="logit", family_part2=gaussian())
                temp.ame <- AME(temp.mod_base_tpm)
                
                prs_base_res_list[[i]] <- c(tissue,pgs_id,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_traits,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_fullname,method,cell.type,p_thres,subtype,stage,dim(temp.dat)[1],"baseline",
                                            "twopartm","linear_logistic",
                                            temp.ame[temp.ame$Variable=="scaled_prs",][2:5],
                                            summary(temp.mod_base_tpm)$Secondpart.model$coefficients["scaled_prs",],
                                            summary(temp.mod_base_tpm)$Firstpart.model$coefficients["scaled_prs",],
                                            rep(NA,2))
                
                i <- i + 1
              }
              
            }
          }
        }
        
        prs_base_res_dat <- as.data.frame(do.call("rbind",prs_base_res_list))
        prs_full_res_dat <- as.data.frame(do.call("rbind",prs_full_res_list))
        
        tmp_colnames <- c("tissue","pgs_id","pgs_trait","pgs_trait_full","method","cell_type","p_threshold","subtype","stage","sample_size","model",
                          "model_family","model_formula",
                          "linear_beta_prs","linear_se_prs","linear_t_prs","linear_p_prs",
                          "tpm_linear_beta_prs","tpm_linear_se_prs","tpm_linear_t_prs","tpm_linear_p_prs",
                          "tpm_logistic_beta_prs","tpm_logistic_se_prs","tpm_logistic_t_prs","tpm_logistic_p_prs",
                          "linear_beta_age","linear_p_age")
        
        colnames(prs_base_res_dat) <- tmp_colnames
        colnames(prs_full_res_dat) <- tmp_colnames
        
        prs_base_res_dat <- prs_base_res_dat %>% mutate(model_formula = sapply(model_formula, toString))
        prs_full_res_dat <- prs_full_res_dat %>% mutate(model_formula = sapply(model_formula, toString))
        
        fwrite(prs_base_res_dat, paste0("results/model/",method,"_prs_basemod_p.",p_thres,"_subtype.",subtype,"_stage.",stage,".csv"), row.names=F)
        fwrite(prs_full_res_dat, paste0("results/model/",method,"_prs_fullmod_p.",p_thres,"_subtype.",subtype,"_stage.",stage,".csv"), row.names=F)
      }
    }
  }
}

##---------------------------------------------------------
## TIMER 2.0, XCell, EPIC, quanTIseq, MCPcounter
##---------------------------------------------------------

other_methods <- c("timer2","xcell","epic","quantiseq","mcpcounter")
other_data_list <- list(timer2_all_dat,
                        xcell_all_dat,
                        epic_all_dat,
                        quantiseq_all_dat,
                        mcpcounter_all_dat)
names(other_data_list) <- other_methods
other_selected_list <- list(timer2.selected,
                            xcell.selected,
                            epic.selected,
                            quantiseq.selected,
                            mcpcounter.selected) 
names(other_selected_list) <- other_methods

for (method in other_methods){
  
  all.dat <- other_data_list[[method]]
  
  for (subtype in c("all","erp","ern","LumA")){
    for(stage in c("all","stage12","stage34")){
      
      prs_base_res_list <- list()
      prs_full_res_list <- list()
      i <- 1
      j <- 1
      
      for (tissue in c("tumor","normal")){
        
        temp.dat <- all.dat[all.dat$TissueType==tissue,]
        
        if (stage=="all"){
          temp.dat <- temp.dat
        } else if (stage=="stage12"){
          temp.dat <- temp.dat[temp.dat$stage<=2,]
        } else if (stage=="stage34") {
          temp.dat <- temp.dat[temp.dat$stage>=3,]
        } 
        
        if (subtype=="all"){
          temp.dat <- temp.dat
        } else if (subtype=="erp"){
          temp.dat <- temp.dat[temp.dat$er==1,]
        } else if (subtype=="ern"){
          temp.dat <- temp.dat[temp.dat$er==0,]
        } else if (subtype %in% c("LumA")){
          temp.dat <- temp.dat[temp.dat$subtype_PAM50_parker==subtype,]
        }
        
        for (pgs_id in pgs_ids){
          
          scaled_prs <- temp.dat[[pgs_id]]
          
          for (cell.type in other_selected_list[[method]]){
            
            infiltrates <- qnorm((rank(temp.dat[[cell.type]],na.last="keep")-0.5)/sum(!is.na(temp.dat[[cell.type]])))
            mod.family <- gaussian()
            temp.dat$CombinedBatch <- as.factor(temp.dat$CombinedBatch)
            
            if (tissue == "tumor"){
              
              temp.dat$stage_factor <- ifelse(temp.dat$stage %in% c(1,2), "Stage I/II", "Stage III/IV")
              temp.dat$grade <- as.factor(temp.dat$grade)
              temp.dat$er_factor <- ifelse(temp.dat$er==1,"ER+","ER-")
              
              if (pgs_id!="PGS000841"){
                temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                          onebef_bmi,platform,stage_factor,grade,er_factor,CombinedBatch))
              } else {
                temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                          platform,stage_factor,grade,er_factor,CombinedBatch))
              }
              
            } else if (tissue == "normal"){
              if (pgs_id!="PGS000841"){
                temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                          onebef_bmi,platform,CombinedBatch))
              } else{
                temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                          platform,CombinedBatch))
              }
              
            }
            
            if (tissue=="tumor"){
              if (subtype %in% c("erp","ern")){temp.dat_sub <- subset(temp.dat_sub, select=-er_factor)}
              if (stage!="all") {temp.dat_sub <- subset(temp.dat_sub, select=-stage_factor)}
            }
            
            temp.mod_base <- glm(infiltrates ~ scaled_prs + agedx + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family=mod.family, data=temp.dat)
            temp.mod_full <- glm(infiltrates ~ scaled_prs + ., family=mod.family, data=temp.dat_sub)
            
            prs_base_res_list[[i]] <- c(tissue,pgs_id,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_traits,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_fullname,method,cell.type,NA,subtype,stage,dim(temp.dat)[1],"baseline",
                                        as.character(mod.family[1]),formula(temp.mod_base),
                                        summary(temp.mod_base)$coefficients["scaled_prs",],
                                        rep(NA,8),
                                        summary(temp.mod_base)$coefficients["agedx",][c(1,4)])
            prs_full_res_list[[j]] <- c(tissue,pgs_id,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_traits,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_fullname,method,cell.type,NA,subtype,stage,dim(temp.dat_sub)[1],"full",
                                        as.character(mod.family[1]),formula(temp.mod_full),
                                        summary(temp.mod_full)$coefficients["scaled_prs",],
                                        rep(NA,8),
                                        summary(temp.mod_full)$coefficients["agedx",][c(1,4)])
            i <- i + 1
            j <- j + 1
            
            if (method=="xcell" & cell.type %in% xcell.selected_twopart & !(subtype == "ern" & stage == "stage34") & !(subtype == "LumA" & stage == "stage34")){
              
              tpm.dat <- temp.dat
              tpm.dat$infiltrates <- temp.dat[[cell.type]]
              tpm.dat$scaled_prs <- scaled_prs
              
              temp.mod_base_tpm <- tpm(infiltrates ~ scaled_prs + agedx + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=tpm.dat, link_part1="logit", family_part2=gaussian())
              temp.ame <- AME(temp.mod_base_tpm)
              
              prs_base_res_list[[i]] <- c(tissue,pgs_id,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_traits,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_fullname,method,cell.type,NA,subtype,stage,dim(temp.dat)[1],"baseline",
                                          "twopartm","linear_logistic",
                                          temp.ame[temp.ame$Variable=="scaled_prs",][2:5],
                                          summary(temp.mod_base_tpm)$Secondpart.model$coefficients["scaled_prs",],
                                          summary(temp.mod_base_tpm)$Firstpart.model$coefficients["scaled_prs",],
                                          rep(NA,2))
              
              i <- i + 1
            }
            
          }
        }
      }
      
      prs_base_res_dat <- as.data.frame(do.call("rbind",prs_base_res_list))
      prs_full_res_dat <- as.data.frame(do.call("rbind",prs_full_res_list))
      
      tmp_colnames <- c("tissue","pgs_id","pgs_trait","pgs_trait_full","method","cell_type","p_threshold","subtype","stage","sample_size","model",
                        "model_family","model_formula",
                        "linear_beta_prs","linear_se_prs","linear_t_prs","linear_p_prs",
                        "tpm_linear_beta_prs","tpm_linear_se_prs","tpm_linear_t_prs","tpm_linear_p_prs",
                        "tpm_logistic_beta_prs","tpm_logistic_se_prs","tpm_logistic_t_prs","tpm_logistic_p_prs",
                        "linear_beta_age","linear_p_age")
      
      colnames(prs_base_res_dat) <- tmp_colnames
      colnames(prs_full_res_dat) <- tmp_colnames
      
      prs_base_res_dat <- prs_base_res_dat %>% mutate(model_formula = sapply(model_formula, toString))
      prs_full_res_dat <- prs_full_res_dat %>% mutate(model_formula = sapply(model_formula, toString))
      
      fwrite(prs_base_res_dat, paste0("results/model/",method,"_prs_basemod_subtype.",subtype,"_stage.",stage,".csv"), row.names=F)
      fwrite(prs_full_res_dat, paste0("results/model/",method,"_prs_fullmod_subtype.",subtype,"_stage.",stage,".csv"), row.names=F)
    }
  }
}

##---------------------------------------------------------
## Wolf Signatures
##---------------------------------------------------------

for (subtype in c("all","erp","ern","LumA")){
  for(stage in c("all","stage12","stage34")){
    
    prs_base_res_list <- list()
    prs_full_res_list <- list()
    i <- 1
    
    for (tissue in c("tumor","normal")){
      
      temp.dat <- wolf_all_dat[wolf_all_dat$TissueType==tissue,]
      
      if (stage=="all"){
        temp.dat <- temp.dat
      } else if (stage=="stage12"){
        temp.dat <- temp.dat[temp.dat$stage<=2,]
      } else if (stage=="stage34") {
        temp.dat <- temp.dat[temp.dat$stage>=3,]
      } 
      
      if (subtype=="all"){
        temp.dat <- temp.dat
      } else if (subtype=="erp"){
        temp.dat <- temp.dat[temp.dat$er==1,]
      } else if (subtype=="ern"){
        temp.dat <- temp.dat[temp.dat$er==0,]
      } else if (subtype %in% c("LumA")){
        temp.dat <- temp.dat[temp.dat$subtype_PAM50_parker==subtype,]
      }
      
      for (pgs_id in pgs_ids){
        
        scaled_prs <- temp.dat[[pgs_id]]
        
        for (cell.type in wolf.selected){
          
          infiltrates <- qnorm((rank(temp.dat[[cell.type]],na.last="keep")-0.5)/sum(!is.na(temp.dat[[cell.type]])))
          mod.family <- gaussian()
          temp.dat$CombinedBatch <- as.factor(temp.dat$CombinedBatch)
          
          if (tissue == "tumor"){
            
            temp.dat$stage_factor <- ifelse(temp.dat$stage %in% c(1,2), "Stage I/II", "Stage III/IV")
            temp.dat$grade <- as.factor(temp.dat$grade)
            temp.dat$er_factor <- ifelse(temp.dat$er==1,"ER+","ER-")
            
            if (pgs_id!="PGS000841"){
              temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                        onebef_bmi,platform,stage_factor,grade,er_factor,CombinedBatch))
            } else {
              temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                        platform,stage_factor,grade,er_factor,CombinedBatch))
            }
            
          } else if (tissue == "normal"){
            if (pgs_id!="PGS000841"){
              temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                        onebef_bmi,platform,CombinedBatch))
            } else{
              temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                        platform,CombinedBatch))
            }
            
          }
          
          if (tissue=="tumor"){
            if (subtype %in% c("erp","ern")){temp.dat_sub <- subset(temp.dat_sub, select=-er_factor)}
            if (stage!="all") {temp.dat_sub <- subset(temp.dat_sub, select=-stage_factor)}
          }
          
          temp.mod_base <- glm(infiltrates ~ scaled_prs + agedx + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family=mod.family, data=temp.dat)
          temp.mod_full <- glm(infiltrates ~ scaled_prs + ., family=mod.family, data=temp.dat_sub)
          
          prs_base_res_list[[i]] <- c(tissue,pgs_id,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_traits,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_fullname,"wolf",cell.type,NA,subtype,stage,dim(temp.dat)[1],"baseline",
                                      as.character(mod.family[1]),formula(temp.mod_base),
                                      summary(temp.mod_base)$coefficients["scaled_prs",],
                                      rep(NA,8),
                                      summary(temp.mod_base)$coefficients["agedx",][c(1,4)])
          prs_full_res_list[[i]] <- c(tissue,pgs_id,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_traits,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_fullname,"wolf",cell.type,NA,subtype,stage,dim(temp.dat_sub)[1],"full",
                                      as.character(mod.family[1]),formula(temp.mod_full),
                                      summary(temp.mod_full)$coefficients["scaled_prs",],
                                      rep(NA,8),
                                      summary(temp.mod_full)$coefficients["agedx",][c(1,4)])
          i <- i + 1
          
        }
      }
    }
    
    prs_base_res_dat <- as.data.frame(do.call("rbind",prs_base_res_list))
    prs_full_res_dat <- as.data.frame(do.call("rbind",prs_full_res_list))
    
    tmp_colnames <- c("tissue","pgs_id","pgs_trait","pgs_trait_full","method","cell_type","p_threshold","subtype","stage","sample_size","model",
                      "model_family","model_formula",
                      "linear_beta_prs","linear_se_prs","linear_t_prs","linear_p_prs",
                      "tpm_linear_beta_prs","tpm_linear_se_prs","tpm_linear_t_prs","tpm_linear_p_prs",
                      "tpm_logistic_beta_prs","tpm_logistic_se_prs","tpm_logistic_t_prs","tpm_logistic_p_prs",
                      "linear_beta_age","linear_p_age")
    
    colnames(prs_base_res_dat) <- tmp_colnames
    colnames(prs_full_res_dat) <- tmp_colnames
    
    prs_base_res_dat <- prs_base_res_dat %>% mutate(model_formula = sapply(model_formula, toString))
    prs_full_res_dat <- prs_full_res_dat %>% mutate(model_formula = sapply(model_formula, toString))
    
    fwrite(prs_base_res_dat, paste0("results/model/wolf_prs_basemod_subtype.",subtype,"_stage.",stage,".csv"), row.names=F)
    fwrite(prs_full_res_dat, paste0("results/model/wolf_prs_fullmod_subtype.",subtype,"_stage.",stage,".csv"), row.names=F)
  }
}

##---------------------------------------------------------
## IHC 
##---------------------------------------------------------

for (subtype in c("all","erp","ern","LumA")){
  for(stage in c("all","stage12")){
    
    prs_base_res_list <- list()
    prs_full_res_list <- list()
    i <- 1
    
    temp.dat <- ihc_all_dat 
    
    if (stage=="all"){
      temp.dat <- temp.dat
    } else if (stage=="stage12"){
      temp.dat <- temp.dat[temp.dat$stage<=2,]
    } else if (stage=="stage34") {
      temp.dat <- temp.dat[temp.dat$stage>=3,]
    } 
    
    if (subtype=="all"){
      temp.dat <- temp.dat
    } else if (subtype=="erp"){
      temp.dat <- temp.dat[temp.dat$er==1,]
    } else if (subtype=="ern"){
      temp.dat <- temp.dat[temp.dat$er==0,]
    } else if (subtype %in% c("LumA")){
      temp.dat <- temp.dat[temp.dat$subtype_PAM50_parker==subtype,]
    }
    
    for (pgs_id in pgs_ids){
      
      scaled_prs <- temp.dat[[pgs_id]]
      
      for (cell.type in ihc.selected){
        
        infiltrates <- qnorm((rank(temp.dat[[cell.type]],na.last="keep")-0.5)/sum(!is.na(temp.dat[[cell.type]])))
        mod.family <- gaussian()
        temp.dat$CombinedBatch <- as.factor(temp.dat$CombinedBatch)
        temp.dat$stage_factor <- ifelse(temp.dat$stage %in% c(1,2), "Stage I/II", "Stage III/IV")
        temp.dat$grade <- as.factor(temp.dat$grade)
        temp.dat$er_factor <- ifelse(temp.dat$er==1,"ER+","ER-")
        
        if (pgs_id!="PGS000841"){
          temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                    onebef_bmi,platform,stage_factor,grade,er_factor,CombinedBatch))
        } else {
          temp.dat_sub <- subset(temp.dat, select=c(agedx,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                                    platform,stage_factor,grade,er_factor,CombinedBatch))
        }
        
        if (subtype %in% c("erp","ern")){temp.dat_sub <- subset(temp.dat_sub, select=-er_factor)}
        if (stage!="all") {temp.dat_sub <- subset(temp.dat_sub, select=-stage_factor)}
        
        temp.mod_base <- glm(infiltrates ~ scaled_prs + agedx + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family=mod.family, data=temp.dat)
        temp.mod_full <- glm(infiltrates ~ scaled_prs + ., family=mod.family, data=temp.dat_sub)
        
        prs_base_res_list[[i]] <- c("tumor",pgs_id,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_traits,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_fullname,"ihc",cell.type,NA,subtype,stage,dim(temp.dat)[1],"baseline",
                                    as.character(mod.family[1]),formula(temp.mod_base),
                                    summary(temp.mod_base)$coefficients["scaled_prs",],
                                    rep(NA,8),
                                    summary(temp.mod_base)$coefficients["agedx",][c(1,4)])
        prs_full_res_list[[i]] <- c("tumor",pgs_id,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_traits,pgs_dict[pgs_dict$pgs_ids==pgs_id,]$pgs_fullname,"ihc",cell.type,NA,subtype,stage,dim(temp.dat_sub)[1],"full",
                                    as.character(mod.family[1]),formula(temp.mod_full),
                                    summary(temp.mod_full)$coefficients["scaled_prs",],
                                    rep(NA,8),
                                    summary(temp.mod_full)$coefficients["agedx",][c(1,4)])
        i <- i + 1
        
      }
    }
    
    prs_base_res_dat <- as.data.frame(do.call("rbind",prs_base_res_list))
    prs_full_res_dat <- as.data.frame(do.call("rbind",prs_full_res_list))
    
    tmp_colnames <- c("tissue","pgs_id","pgs_trait","pgs_trait_full","method","cell_type","p_threshold","subtype","stage","sample_size","model",
                      "model_family","model_formula",
                      "linear_beta_prs","linear_se_prs","linear_t_prs","linear_p_prs",
                      "tpm_linear_beta_prs","tpm_linear_se_prs","tpm_linear_t_prs","tpm_linear_p_prs",
                      "tpm_logistic_beta_prs","tpm_logistic_se_prs","tpm_logistic_t_prs","tpm_logistic_p_prs",
                      "linear_beta_age","linear_p_age")
    
    colnames(prs_base_res_dat) <- tmp_colnames
    colnames(prs_full_res_dat) <- tmp_colnames
    
    prs_base_res_dat <- prs_base_res_dat %>% mutate(model_formula = sapply(model_formula, toString))
    prs_full_res_dat <- prs_full_res_dat %>% mutate(model_formula = sapply(model_formula, toString))
    
    fwrite(prs_base_res_dat, paste0("results/model/ihc_prs_basemod_subtype.",subtype,"_stage.",stage,".csv"), row.names=F)
    fwrite(prs_full_res_dat, paste0("results/model/ihc_prs_fullmod_subtype.",subtype,"_stage.",stage,".csv"), row.names=F)
  }
}
