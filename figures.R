library(data.table)
library(ggplot2)
library(stringr)
library(dplyr)

##--------------------------------------------
## Data for figures
##--------------------------------------------

source("scripts/display_dict.R")

### Primary results

prs_immune_res_all <- read.csv("results/model/anno_immune_prs_results.csv")
primary_res <- read.csv("results/model/anno_fullmod_immune_prs_results.csv")
primary_res <- merge(primary_res, tissue_dict, by="tissue")
primary_res <- merge(primary_res, pgs_dict, by="pgs_trait")
primary_res <- merge(primary_res, method_dict, by="method")
primary_res <- merge(primary_res, cell_dict, by="cell_type")

primary_res$tissue_full <- factor(primary_res$tissue_full, 
                                  levels=c("Tumor","Tumor-adjacent normal tissue"))

primary_res$pgs_full <- factor(primary_res$pgs_full, 
                               levels=c("Breast cancer (overall)","Breast Cancer (ER+)","Breast Cancer (ER-)",
                                        "IBD","UC","CD","RA","Allergic disease","Asthma",
                                        "T2D",
                                        "BMI","BMI-adjusted WHR",
                                        "Age at menarche","Age at menopause",
                                        "Cigarettes per day","Alcohol consumption"))

primary_res$pgs_group <- factor(primary_res$pgs_group, 
                                levels=c("Breast cancer",
                                         "Immune-related diseases",
                                         "BMI-related diseases/traits",
                                         "Hormone-related traits",
                                         "Lifestyle factors"))

primary_res$method_full <- factor(primary_res$method_full, 
                                  levels=c("Expression signature","IHC","CIBERSORTx",
                                           "TIMER 2.0","xCell","EPIC","MCPcounter","quanTIseq"))


unique(primary_res$pgs_full)
unique(primary_res$pgs_group)
unique(primary_res$method_full)
unique(primary_res$cell_type_full)

### Replication and meta-analysis

pgs_trait_selected <- c("ibd","uc","cd")
cell_type_seleced <- c("GP11_Immune_IFN","IFN_21978456","IFNG_score_21050467","Immune_NSCLC_score","Interferon_19272155","Interferon_Cluster_21214954","MHC1_21978456","Minterferon_Cluster_21214954","Module3_IFN_score","NHI_5gene_score","STAT1_19272155","STAT1_score")

nhs_nbcs <- read.csv("results/model/validation/combined_results_nhs_nbcs.csv") %>% 
  subset(select = c(study, cell_type, pgs_trait, linear_beta_prs, linear_se_prs, linear_p_prs)) %>%
  filter(pgs_trait %in% pgs_trait_selected,
         cell_type %in% cell_type_seleced)
colnames(nhs_nbcs) <- c("study", "cell_type", "pgs_trait", "beta", "se", "p")

meta_res_dat <- read.csv("results/model/validation/meta_nhs_nbcs.csv") %>% 
  subset(select = c(cell_type, pgs_trait, beta_fixed, se_fixed, p_fixed)) %>%
  mutate(study = "meta") %>%
  filter(pgs_trait %in% pgs_trait_selected,
         cell_type %in% cell_type_seleced)
colnames(meta_res_dat) <- c("cell_type", "pgs_trait", "beta", "se", "p", "study")

combined_res <- as.data.frame(rbind(nhs_nbcs, meta_res_dat))
combined_res <- merge(combined_res, pgs_dict, by="pgs_trait")
combined_res <- merge(combined_res, cell_dict, by="cell_type")

combined_res$pgs_full <- factor(combined_res$pgs_full, 
                                levels=c("IBD","UC","CD"))
combined_res <- combined_res %>%
  mutate(study_full = case_when(study == "meta" ~ "Meta-analysis",
                                study == "nhs" ~ "NHS/NHSII",
                                study == "nbcs" ~ "NBCS"))
combined_res$study_full <- factor(combined_res$study_full, 
                                  levels=c("NHS/NHSII","NBCS","Meta-analysis"))

##--------------------------------------------
## Figure 2
##--------------------------------------------

### Heatmap (tumor and normal)

temp <- primary_res
temp$sig <- ifelse(temp$fdr_byprs < 0.25, "*", "")
temp$sig <- ifelse(temp$fdr_byprs < 0.1, "**", temp$sig)
temp$sig <- ifelse(temp$fdr_byprs < 0.05, "***", temp$sig)
sig_features <- unique(temp[temp$fdr_byprs < 0.25,]$method_celltype)
temp <- temp[temp$method_celltype %in% sig_features,]
sum(temp$sig!="")

tmp_heatmap <- ggplot(temp, aes(x=pgs_full, y=cell_type_full)) +
  geom_tile(aes(fill=linear_beta_prs)) +
  scale_fill_gradient2(low="#466983",
                       mid="white",
                       high="#8b0000",
                       midpoint = 0,
                       guide = "colourbar") +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(limits=rev) +
  geom_text(aes(label=sig), color="black", size=9, nudge_y = -0.2) +
  labs(fill="Effect size", size=12) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement="outside",
        strip.text.y = element_text(size=18, angle=0, face="bold"),
        strip.text.x = element_text(size=22, color="white", face="bold"),
        strip.background.x = element_rect(color="black", fill=c("#8b0000","#466983"), size=1.5, linetype="solid")) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black", size = 18),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 18, colour = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5)) +
  theme(legend.position = "right", legend.title = element_text(size=16), legend.text = element_text(size=16), legend.key.size = unit(1, 'cm')) +
  facet_grid(method_full~tissue_full,scales="free_y",space="free_y") 

pdf("results/plots/main/germline_immune_manuscript_heatmap_all_fullmod_fdrbyprs_0.25_new_color.pdf", width=21, height=15)
tmp_heatmap
dev.off() 

##--------------------------------------------
## Figure 3
##--------------------------------------------

## Forest plot (replication and meta-analysis)

tmp_trait <- "ibd"
tmp_dat <- combined_res %>%
  mutate(beta = as.numeric(beta),
         se = as.numeric(se)) %>%
  mutate(lower = beta - 1.96*se,
         upper = beta + 1.96*se) %>%
  filter(pgs_trait == tmp_trait)

p1 <- ggplot(data=tmp_dat, aes(y=cell_type_full, x=beta, xmin=lower, xmax=upper, fill=study_full, shape=study_full)) +
  geom_linerange(size=0.4,position=position_dodge(width = 0.5)) +
  geom_point(size=2, colour="white", stroke = 0.5, position=position_dodge(width = 0.5), show.legend=F) +
  labs(title=toupper(tmp_trait), x='Effect Size', y='') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  scale_y_discrete(limits=rev) +
  scale_fill_manual(values=c("#5DB1DD","#CC9900","gray29")) +
  scale_shape_manual(values = c(21,21,22)) +
  theme_bw()+
  theme(legend.title= element_blank(),
        plot.title = element_text(size = 10, face = "bold"))

tmp_trait <- "cd"
tmp_dat <- combined_res %>%
  mutate(beta = as.numeric(beta),
         se = as.numeric(se)) %>%
  mutate(lower = beta - 1.96*se,
         upper = beta + 1.96*se) %>%
  filter(pgs_trait == tmp_trait)

p2 <- ggplot(data=tmp_dat, aes(y=cell_type_full, x=beta, xmin=lower, xmax=upper, fill=study_full, shape=study_full)) +
  geom_linerange(size=0.4,position=position_dodge(width = 0.5)) +
  geom_point(size=2, colour="white", stroke = 0.5, position=position_dodge(width = 0.5), show.legend=F) +
  labs(title=toupper(tmp_trait), x='Effect Size', y='') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  scale_y_discrete(limits=rev) +
  scale_fill_manual(values=c("#5DB1DD","#CC9900","gray29")) +
  scale_shape_manual(values = c(21,21,22)) +
  theme_bw()+
  theme(legend.title= element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_blank())

tmp_trait <- "uc"
tmp_dat <- combined_res %>%
  mutate(beta = as.numeric(beta),
         se = as.numeric(se)) %>%
  mutate(lower = beta - 1.96*se,
         upper = beta + 1.96*se) %>%
  filter(pgs_trait == tmp_trait)

p3 <- ggplot(data=tmp_dat, aes(y=cell_type_full, x=beta, xmin=lower, xmax=upper, fill=study_full, shape=study_full)) +
  geom_linerange(size=0.4,position=position_dodge(width = 0.5)) +
  geom_point(size=2, colour="white", stroke = 0.5, position=position_dodge(width = 0.5)) +
  labs(title=toupper(tmp_trait), x='Effect Size', y='') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  scale_y_discrete(limits=rev) +
  scale_fill_manual(values=c("#5DB1DD","#CC9900","gray29")) +
  scale_shape_manual(values = c(21,21,22)) +
  theme_bw()+
  theme(legend.title= element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_blank())

pdf("results/plots/main/germline_immune_manuscript_forest_replication_meta.pdf", width = 10, height = 6, onefile = F) 
egg::ggarrange(p1,p2,p3, nrow=1, ncol=3,
               widths = c(1,1,1), heights=c(10))
dev.off() 