library(data.table)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra) 
library(ggpubr)
library(gplots)
library(tidyverse)

##--------------------------------------------
## Figure S1 and Figure S2
##--------------------------------------------

source("scripts/display_dict.R")
source("scripts/immune_dict.R")

prs_covar_dat <- read.csv("data/processed_1000G_PRS_var_PC_covar.txt")

##-------------------------------
## By methods
##-------------------------------

##-----------------
## CIBERSORTx
##-----------------

### Data

cibersortx <- read.csv("data/immune/anno_cibersortx_all.csv")
cibersortx <- cibersortx[cibersortx$X %in% prs_covar_dat$X,]
cibersortx_immune <- subset(cibersortx, select=lm22.selected)

sum(colnames(cibersortx_immune) %in% cell_dict$cell_type)
tmp <- cell_dict[cell_dict$cell_type %in% colnames(cibersortx_immune),]
setnames(cibersortx_immune, tmp$cell_type, tmp$cell_type_full)

cibersortx_immune.tumor <- cibersortx_immune[cibersortx$tissue=="tumor",]
cibersortx_immune.normal <- cibersortx_immune[cibersortx$tissue=="normal",]

### Heatmap

# Tumor - Normal 
corr1 <- round(cor(cibersortx_immune.tumor), 2)
corr2 <- round(cor(cibersortx_immune.normal), 2)
new <- corr1
new[upper.tri(new)] <- corr2[upper.tri(corr2)]
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/supp/heatmap_cluster_cibersortx_tumor_normal.pdf", width=6, height=6)
heatmap(x = new, col = col, symm = T, margins = c(12, 12), xlab="Tumor", ylab="Normal")
dev.off() 

##-----------------
## TIMER 2.0
##-----------------

### Data

timer2 <- read.csv("data/immune/anno_timer2_all.csv")
timer2 <- timer2[timer2$X %in% prs_covar_dat$X,]
timer2_immune <- subset(timer2, select=timer2.selected)

sum(colnames(timer2_immune) %in% cell_dict$cell_type)
tmp <- cell_dict[cell_dict$cell_type %in% colnames(timer2_immune),]
setnames(timer2_immune, tmp$cell_type, tmp$cell_type_full)

timer2_immune.tumor <- timer2_immune[timer2$tissue=="tumor",]
timer2_immune.normal <- timer2_immune[timer2$tissue=="normal",]

### Heatmap

# All 
corr <- round(cor(timer2_immune), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_timer2_all.pdf", width=6, height=6)
heatmap(x = corr, col = col, symm = TRUE, margins = c(15, 15))
dev.off() 

# Tumor 
corr <- round(cor(timer2_immune.tumor), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_timer2_tumor.pdf", width=6, height=6)
heatmap(x = corr, col = col, symm = TRUE, margins = c(15, 15))
dev.off() 

# Normal 
corr <- round(cor(timer2_immune.normal), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_timer2_normal.pdf", width=6, height=6)
heatmap(x = corr, col = col, symm = TRUE, margins = c(15, 15))
dev.off() 

# Tumor - Normal 
corr1 <- round(cor(timer2_immune.tumor), 2)
corr2 <- round(cor(timer2_immune.normal), 2)
new <- corr1
new[upper.tri(new)] <- corr2[upper.tri(corr2)]
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/supp/heatmap_cluster_timer2_tumor_normal.pdf", width=6, height=6)
heatmap(x = new, col = col, symm = TRUE, margins = c(14, 14), xlab="Tumor", ylab="Normal")
dev.off() 

##-----------------
## xCell
##-----------------

### Data

xcell <- read.csv("data/immune/anno_xcell_all.csv")
xcell <- xcell[xcell$X %in% prs_covar_dat$X,]
xcell_immune <- subset(xcell, select=xcell.selected)

sum(colnames(xcell_immune) %in% cell_dict$cell_type)
tmp <- cell_dict[cell_dict$cell_type %in% colnames(xcell_immune),]
setnames(xcell_immune, tmp$cell_type, tmp$cell_type_full)

xcell_immune.tumor <- xcell_immune[xcell$tissue=="tumor",]
xcell_immune.normal <- xcell_immune[xcell$tissue=="normal",]

### Heatmap

# All 
corr <- round(cor(xcell_immune), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_xcell_all.pdf", width=10, height=10)
heatmap(x = corr, col = col, symm = TRUE, margins = c(15, 15))
dev.off() 

# Tumor 
corr <- round(cor(xcell_immune.tumor), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_xcell_tumor.pdf", width=10, height=10)
heatmap(x = corr, col = col, symm = TRUE, margins = c(15, 15))
dev.off() 

# Normal 
corr <- round(cor(xcell_immune.normal), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_xcell_normal.pdf", width=10, height=10)
heatmap(x = corr, col = col, symm = TRUE, margins = c(15, 15))
dev.off() 

# Tumor - Normal 
corr1 <- round(cor(xcell_immune.tumor), 2)
corr2 <- round(cor(xcell_immune.normal), 2)
new <- corr1
new[upper.tri(new)] <- corr2[upper.tri(corr2)]
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/supp/heatmap_cluster_xcell_tumor_normal.pdf", width=10, height=10)
heatmap(x = new, col = col, symm = TRUE, margins = c(15, 15), xlab="Tumor", ylab="Normal")
dev.off() 

##-----------------
## EPIC
##-----------------

### Data

epic <- read.csv("data/immune/anno_epic_all.csv")
epic <- epic[epic$X %in% prs_covar_dat$X,]
epic_immune <- subset(epic, select=epic.selected)

sum(colnames(epic_immune) %in% cell_dict$cell_type)
tmp <- cell_dict[cell_dict$cell_type %in% colnames(epic_immune),]
setnames(epic_immune, tmp$cell_type, tmp$cell_type_full)

epic_immune.tumor <- epic_immune[epic$tissue=="tumor",]
epic_immune.normal <- epic_immune[epic$tissue=="normal",]

### Heatmap

# All 
corr <- round(cor(epic_immune), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_epic_all.pdf", width=6, height=6)
heatmap(x = corr, col = col, symm = TRUE, margins = c(15, 15))
dev.off() 

# Tumor 
corr <- round(cor(epic_immune.tumor), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_epic_tumor.pdf", width=6, height=6)
heatmap(x = corr, col = col, symm = TRUE, margins = c(15, 15))
dev.off() 

# Normal 
corr <- round(cor(epic_immune.normal), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_epic_normal.pdf", width=6, height=6)
heatmap(x = corr, col = col, symm = TRUE, margins = c(15, 15))
dev.off() 

# Tumor - Normal 
corr1 <- round(cor(epic_immune.tumor), 2)
corr2 <- round(cor(epic_immune.normal), 2)
new <- corr1
new[upper.tri(new)] <- corr2[upper.tri(corr2)]
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/supp/heatmap_cluster_epic_tumor_normal.pdf", width=6, height=6)
heatmap(x = new, col = col, symm = TRUE, margins = c(15, 15), xlab="Tumor", ylab="Normal")
dev.off() 

##-----------------
## MCPcounter
##-----------------

### Data

mcpcounter <- read.csv("data/immune/anno_mcpcounter_all.csv")
mcpcounter <- mcpcounter[mcpcounter$X %in% prs_covar_dat$X,]
mcpcounter_immune <- subset(mcpcounter, select=mcpcounter.selected)

sum(colnames(mcpcounter_immune) %in% cell_dict$cell_type)
tmp <- cell_dict[cell_dict$cell_type %in% colnames(mcpcounter_immune),]
setnames(mcpcounter_immune, tmp$cell_type, tmp$cell_type_full)

mcpcounter_immune.tumor <- mcpcounter_immune[mcpcounter$tissue=="tumor",]
mcpcounter_immune.normal <- mcpcounter_immune[mcpcounter$tissue=="normal",]

### Heatmap

# All 
corr <- round(cor(mcpcounter_immune), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_mcpcounter_all.pdf", width=6, height=6)
heatmap(x = corr, col = col, symm = TRUE, margins = c(12, 12))
dev.off() 

# Tumor 
corr <- round(cor(mcpcounter_immune.tumor), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_mcpcounter_tumor.pdf", width=6, height=6)
heatmap(x = corr, col = col, symm = TRUE, margins = c(12, 12))
dev.off() 

# Normal 
corr <- round(cor(mcpcounter_immune.normal), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_mcpcounter_normal.pdf", width=6, height=6)
heatmap(x = corr, col = col, symm = TRUE, margins = c(12, 12))
dev.off() 

# Tumor - Normal 
corr1 <- round(cor(mcpcounter_immune.tumor), 2)
corr2 <- round(cor(mcpcounter_immune.normal), 2)
new <- corr1
new[upper.tri(new)] <- corr2[upper.tri(corr2)]
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/supp/heatmap_cluster_mcpcounter_tumor_normal.pdf", width=6, height=6)
heatmap(x = new, col = col, symm = TRUE, margins = c(12, 12), xlab="Tumor", ylab="Normal")
dev.off() 

##-----------------
## quantiseq
##-----------------

### Data

quantiseq <- read.csv("data/immune/anno_quantiseq_all.csv")
quantiseq <- quantiseq[quantiseq$X %in% prs_covar_dat$X,]
quantiseq_immune <- subset(quantiseq, select=quantiseq.selected)

sum(colnames(quantiseq_immune) %in% cell_dict$cell_type)
tmp <- cell_dict[cell_dict$cell_type %in% colnames(quantiseq_immune),]
setnames(quantiseq_immune, tmp$cell_type, tmp$cell_type_full)

quantiseq_immune.tumor <- quantiseq_immune[quantiseq$tissue=="tumor",]
quantiseq_immune.normal <- quantiseq_immune[quantiseq$tissue=="normal",]

### Heatmap

# All 
corr <- round(cor(quantiseq_immune), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_quantiseq_all.pdf", width=6, height=6)
heatmap(x = corr, col = col, symm = TRUE, margins = c(12, 12))
dev.off() 

# Tumor 
corr <- round(cor(quantiseq_immune.tumor), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_quantiseq_tumor.pdf", width=6, height=6)
heatmap(x = corr, col = col, symm = TRUE, margins = c(12, 12))
dev.off() 

# Normal 
corr <- round(cor(quantiseq_immune.normal), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_quantiseq_normal.pdf", width=6, height=6)
heatmap(x = corr, col = col, symm = TRUE, margins = c(12, 12))
dev.off() 

# Tumor - Normal 
corr1 <- round(cor(quantiseq_immune.tumor), 2)
corr2 <- round(cor(quantiseq_immune.normal), 2)
new <- corr1
new[upper.tri(new)] <- corr2[upper.tri(corr2)]
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/supp/heatmap_cluster_quantiseq_tumor_normal.pdf", width=6, height=6)
heatmap(x = new, col = col, symm = TRUE, margins = c(12, 12), xlab="Tumor", ylab="Normal")
dev.off() 

##-----------------
## ihc
##-----------------

### Data

ihc <- read.csv("data/merged_ihc.csv")
ihc_immune_tumor <- subset(ihc, select=ihc.selected)

sum(colnames(ihc_immune_tumor) %in% cell_dict$cell_type)
tmp <- cell_dict[cell_dict$cell_type %in% colnames(ihc_immune_tumor),]
setnames(ihc_immune_tumor, tmp$cell_type, tmp$cell_type_full)
ihc_immune_tumor <- subset(ihc_immune_tumor, select=grepl("IHC",colnames(ihc_immune_tumor)))
ihc_immune_tumor <- ihc_immune_tumor[complete.cases(ihc_immune_tumor),]

### Heatmap

# Tumor 
corr <- round(cor(ihc_immune_tumor), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_ihc_tumor.pdf", width=6, height=6)
heatmap(x = corr, col = col, symm = TRUE, margins = c(12, 12))
dev.off() 

##-------------------------------
## Wolf (Figure S1)
##-------------------------------

### Data

wolf <- read.csv("data/immune/anno_wolf_all.csv")
wolf <- wolf[wolf$X %in% prs_covar_dat$X,]
wolf_immune <- subset(wolf, select=wolf.selected)
wolf_immune.tumor <- wolf_immune[wolf$tissue=="tumor",]
wolf_immune.normal <- wolf_immune[wolf$tissue=="normal",]

### Heatmap

# All 
corr <- round(cor(wolf_immune), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_wolf_all.pdf", width=10, height=10)
heatmap(x = corr, col = col, symm = TRUE, margins = c(12, 12))
dev.off() 

# Tumor 
corr <- round(cor(wolf_immune.tumor), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_wolf_tumor.pdf", width=10, height=10)
heatmap(x = corr, col = col, symm = TRUE, margins = c(12, 12))
dev.off() 

# Normal 
corr <- round(cor(wolf_immune.normal), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_wolf_normal.pdf", width=10, height=10)
heatmap(x = corr, col = col, symm = TRUE, margins = c(12, 12))
dev.off() 

# Tumor - Normal 
corr1 <- round(cor(wolf_immune.tumor), 2)
corr2 <- round(cor(wolf_immune.normal), 2)
new <- corr1
new[upper.tri(new)] <- corr2[upper.tri(corr2)]
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/supp/heatmap_cluster_wolf_tumor_normal.pdf", width=12, height=12)
heatmap(x = new, col = col, symm = TRUE, margins = c(14, 14), xlab="Tumor", ylab="Normal")
dev.off() 

##-------------------------------
## Combined (Figure S2)
##-------------------------------

### Data

colnames(cibersortx_immune) <- paste0(colnames(cibersortx_immune)," (CIBERSORTx)")
colnames(timer2_immune) <- paste0(colnames(timer2_immune)," (TIMER 2.0)")
colnames(xcell_immune) <- paste0(colnames(xcell_immune)," (xCell)")
colnames(epic_immune) <- paste0(colnames(epic_immune)," (EPIC)")
colnames(mcpcounter_immune) <- paste0(colnames(mcpcounter_immune)," (MCPcounter)")
colnames(quantiseq_immune) <- paste0(colnames(quantiseq_immune)," (quanTIseq)")

all_immune <- as.data.frame(cbind(cibersortx_immune,
                                  timer2_immune,
                                  xcell_immune,
                                  epic_immune,
                                  mcpcounter_immune,
                                  quantiseq_immune))

tumor_immune <- all_immune[cibersortx$tissue == "tumor",]
normal_immune <- all_immune[cibersortx$tissue == "normal",]

### Heatmap

col_vec <- c(rep("#A7AC36", dim(cibersortx_immune)[2]),
             rep("#282274", dim(timer2_immune)[2]),
             rep("#C75127", dim(xcell_immune)[2]),
             rep("#5DB1DD", dim(epic_immune)[2]),
             rep("#CC9900", dim(mcpcounter_immune)[2]),
             rep("#339900", dim(quantiseq_immune)[2]))

# All 
corr <- round(cor(all_immune), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_all.pdf", width=15, height=15)
heatmap(x = corr, col = col, RowSideColors = col_vec, ColSideColors = col_vec, symm = TRUE, margins = c(12, 12))
dev.off()

# Tumor 
corr <- round(cor(tumor_immune), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_tumor.pdf", width=15, height=15)
heatmap(x = corr, col = col, RowSideColors = col_vec, ColSideColors = col_vec, symm = TRUE, margins = c(12, 12))
dev.off() 

# Normal 
corr <- round(cor(normal_immune), 2)
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/heatmap_cluster_normal.pdf", width=15, height=15)
heatmap(x = corr, col = col, RowSideColors = col_vec, ColSideColors = col_vec, symm = TRUE, margins = c(12, 12))
dev.off() 

# Tumor - Normal 
corr1 <- round(cor(tumor_immune), 2)
corr2 <- round(cor(normal_immune), 2)
new <- corr1
new[upper.tri(new)] <- corr2[upper.tri(corr2)]
col <- colorRampPalette(c("#466983", "white", "#8b0000"))(20)
pdf("results/plots/supp/heatmap_cluster_tumor_normal.pdf", width=15, height=15)
heatmap(x = new, col = col, RowSideColors = col_vec, ColSideColors = col_vec, symm = TRUE, margins = c(12, 12), xlab="Tumor", ylab="Normal")
dev.off() 
pdf("results/plots/supp/heatmap_cluster_tumor_normal_legend.pdf", width=15, height=15)
plot.new()
legend("topleft", title = "Methods", legend=c("CIBERSORTx","TIMER 2.0","xCell","EPIC","MCPcounter","quanTIseq"), 
       fill=c("#A7AC36","#282274","#C75127","#5DB1DD","#CC9900","#339900"), cex=1.5, box.lty=0)
dev.off() 

##--------------------------------------------
## Figure S3
##--------------------------------------------

source("scripts/display_dict.R")
prs_immune_res_all <- merge(prs_immune_res_all, tissue_dict, by="tissue")
prs_immune_res_all <- merge(prs_immune_res_all, pgs_dict, by="pgs_trait")
prs_immune_res_all <- merge(prs_immune_res_all, method_dict, by="method")
prs_immune_res_all <- merge(prs_immune_res_all, cell_dict, by="cell_type")

min(primary_res$fdr_byprs)
table(primary_res$fdr_byprs < 0.25)
primary_res_sig <- primary_res[primary_res$fdr_byprs < 0.25,]
primary_res_sig$method_celltype_prs <- paste0(primary_res_sig$method_celltype,"_",primary_res_sig$pgs_trait)
prs_immune_res_all$method_celltype_prs <- paste0(prs_immune_res_all$method_celltype,"_",prs_immune_res_all$pgs_trait)

tmp <- primary_res_sig[duplicated(primary_res_sig$method_celltype_prs),] # sig in both tumor and normal

for (mcp in unique(primary_res_sig$method_celltype_prs)){
  
  temp <- prs_immune_res_all[prs_immune_res_all$method_celltype_prs==mcp & prs_immune_res_all$model_family=="gaussian",]
  if (grepl("cibersortx",mcp)){temp <- temp[temp$p_threshold!="none",]}
  
  temp_highlight <- primary_res_sig[primary_res_sig$method_celltype_prs==mcp,]$tissue
  if (length(temp_highlight) > 1) {
    temp_highlight <- c("Tumor","Normal")
  } else if (temp_highlight=="tumor"){
    temp_highlight <- "Tumor"
  } else if (temp_highlight=="normal"){
    temp_highlight <- "Normal"
  } 
  
  tumor_res <- temp[temp$tissue=="tumor" & temp$stage=="all" & temp$subtype=="all" & temp$model=="full",]
  normal_res <- temp[temp$tissue=="normal" & temp$stage=="all" & temp$subtype=="all" & temp$model=="full",]
  erp_tumor_res <- temp[temp$tissue=="tumor" & temp$stage=="all" & temp$subtype=="erp" & temp$model=="full",]
  ern_tumor_res <- temp[temp$tissue=="tumor" & temp$stage=="all" & temp$subtype=="ern" & temp$model=="full",]
  erp_normal_res <- temp[temp$tissue=="normal" & temp$stage=="all" & temp$subtype=="erp" & temp$model=="full",]
  ern_normal_res <- temp[temp$tissue=="normal" & temp$stage=="all" & temp$subtype=="ern" & temp$model=="full",]
  stage12_tumor_res <- temp[temp$tissue=="tumor" & temp$stage=="stage12" & temp$subtype=="all" & temp$model=="full",]
  stage34_tumor_res <- temp[temp$tissue=="tumor" & temp$stage=="stage34" & temp$subtype=="all" & temp$model=="full",]
  stage12_normal_res <- temp[temp$tissue=="normal" & temp$stage=="stage12" & temp$subtype=="all" & temp$model=="full",]
  stage34_normal_res <- temp[temp$tissue=="normal" & temp$stage=="stage34" & temp$subtype=="all" & temp$model=="full",]
  
  tmp_fig <- as.data.frame(rbind(tumor_res,normal_res,erp_tumor_res,ern_tumor_res,erp_normal_res,ern_normal_res,stage12_tumor_res,stage12_normal_res))
  tmp_fig$tissue_full <- str_to_title(tmp_fig$tissue)
  tmp_fig$tissue_full <- factor(tmp_fig$tissue_full, levels=c("Tumor","Normal"))
  tmp_fig$model_full <- c("All","All","ER+","ER-","ER+","ER-","Stage I/II","Stage I/II")
  tmp_fig$model_full <- factor(tmp_fig$model_full, levels=c("All","Stage I/II","ER+","ER-"))
  
  tmp_fig$upper_prs <- tmp_fig$linear_beta_prs + 1.96*tmp_fig$linear_se_prs
  tmp_fig$lower_prs <- tmp_fig$linear_beta_prs - 1.96*tmp_fig$linear_se_prs
  
  ggplot(data=tmp_fig, aes(y=model_full, x=linear_beta_prs, xmin=lower_prs, xmax=upper_prs, fill=tissue_full)) +
    geom_linerange(size=0.4,position=position_dodge(width = 0.5)) +
    geom_point(size=2, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
    labs(x='Effect Size', y='') +
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values=c("darkred", "#466983")) +
    theme_bw()+
    theme(legend.title= element_blank(),
          plot.title = element_text(size = 10, face = "bold")) +
    ggtitle(gsub(" \\(Expression signature\\)","",unique(paste0(tmp_fig$pgs_full," PRS\n",tmp_fig$cell_type_full, " (",tmp_fig$method_full,")"))))
  
  ggsave(paste0("results/plots/main/sensitivity_fig/erp_ern_stage12/",mcp,".pdf"), width=4, height=3)
  
}

##--------------------------------------------
## Figure S4
##--------------------------------------------

sum(colnames(wolf_immune) %in% cell_dict$cell_type)
tmp <- cell_dict[cell_dict$cell_type %in% colnames(wolf_immune),]
setnames(wolf_immune, tmp$cell_type, tmp$cell_type_full)
colnames(wolf_immune) <- paste0(colnames(wolf_immune)," (Expression signature)")
all_immune <- as.data.frame(cbind(cibersortx_immune,
                                  timer2_immune,
                                  xcell_immune,
                                  epic_immune,
                                  mcpcounter_immune,
                                  quantiseq_immune,
                                  wolf_immune))
all_immune$tissue <- ifelse(cibersortx$tissue == "tumor", "tumor", "normal")
all_immune$X <- cibersortx$X
all_dat <- merge(all_immune, prs_covar_dat, by="X")

primary_res <- read.csv("results/model/anno_fullmod_immune_prs_results_09192023.csv")
primary_res <- merge(primary_res, method_dict, by="method")
primary_res <- merge(primary_res, cell_dict, by="cell_type")
sig_primary_res <- primary_res[primary_res$fdr_byprs < 0.25, ]
sig_primary_res <- sig_primary_res[!duplicated(sig_primary_res$method_celltype),]
sig_primary_res$cell_type_full_method_full <- paste0(sig_primary_res$cell_type_full," (",sig_primary_res$method_full,")")

for (cell_type in sig_primary_res$cell_type_full_method_full){
  cell_type_short <- gsub(" \\(Expression signature\\)","",cell_type)
  tmp <- all_dat[,c("NHSID","tissue",cell_type)]
  tmp <- tmp %>% group_by(NHSID) %>% filter(n() == 2)
  tmp$tissue <- ifelse(tmp$tissue=="tumor","Tumor","Normal")
  tmp$value <- qnorm((rank(tmp[[cell_type]],na.last="keep")-0.5)/sum(!is.na(tmp[[cell_type]])))
  ggpaired(tmp, x="tissue", y="value", id="NHSID",
           color="tissue", palette=c("darkred","#466983"), line.color = "gray", line.size = 0.1,) +
    labs(x="", y="") +
    stat_compare_means(method="t.test", paired=T, size=5) +
    theme(legend.position="None",
          plot.title = element_text(size = 15, face = "bold")) +
    ggtitle(cell_type_short)
  ggsave(paste0("results/plots/supp/paired_tumor_normal_",cell_type,".pdf"), width=4.5, height=5)
}

##--------------------------------------------
## Figure S5
##--------------------------------------------

all_immune$X <- cibersortx$X
all_dat <- merge(all_immune, prs_covar_dat, by="X")

all_dat$agedx_factor <- ifelse(all_dat$agedx>65, "Age > 65", "Age <= 65")
all_dat$onebef_bmi_factor <- "Normal"
all_dat$onebef_bmi_factor <- ifelse(all_dat$onebef_bmi>=25,"Overweight",all_dat$onebef_bmi_factor)
all_dat$onebef_bmi_factor <- ifelse(all_dat$onebef_bmi>=30,"Obesity",all_dat$onebef_bmi_factor)
all_dat$onebef_bmi_factor <- as.factor(all_dat$onebef_bmi_factor)
all_dat$stage_factor <- ifelse(all_dat$stage %in% c(1,2), "Stage I/II", "Stage III/IV")
all_dat$er_factor <- ifelse(all_dat$er==1,"ER+","ER-")

table(all_dat$agedx_factor)
table(all_dat$onebef_bmi_factor)
table(all_dat$stage_factor)
table(all_dat$er_factor)

table1_res_list <- list()
i <- 1
for (cell_type in sig_primary_res$cell_type_full_method_full){
  for (tissue in c("tumor","normal")){
    for (x in c("agedx_factor","stage_factor","er_factor","onebef_bmi_factor")){
      tmp <- all_dat[,c(x,cell_type,"tissue")]
      tmp <- tmp[tmp$tissue==tissue,]
      tmp$value <- qnorm((rank(tmp[[cell_type]],na.last="keep")-0.5)/sum(!is.na(tmp[[cell_type]])))
      temp.mod <- glm(value ~ tmp[[x]], data=tmp)
      beta1 <- summary(temp.mod)$coefficients[,"Estimate"][2]
      se1 <- summary(temp.mod)$coefficients[,"Std. Error"][2]
      lower1 <- beta1 - 1.96*se1
      upper1 <- beta1 + 1.96*se1
      table1_res_list[[i]] <- c(cell_type, tissue, x,"beta1", beta1, se1, lower1, upper1)
      if (x=="onebef_bmi_factor"){
        beta2 <- summary(temp.mod)$coefficients[,"Estimate"][3]
        se2 <- summary(temp.mod)$coefficients[,"Std. Error"][3]
        lower2 <- beta2 - 1.96*se2
        upper2 <- beta2 + 1.96*se2
        table1_res_list[[i+1]] <- c(cell_type, tissue, x,"beta2", beta2, se2, lower2, upper2)
        i <- i + 1
      }
      i <- i + 1
    }
  }
}

table1_res <- as.data.frame(do.call("rbind", table1_res_list))
colnames(table1_res) <- c("cell_type","tissue","variable","beta","effect","SE","lower","upper")
table1_res$effect <- as.numeric(table1_res$effect)
table1_res$lower <- as.numeric(table1_res$lower)
table1_res$upper <- as.numeric(table1_res$upper)
table1_res$variable_beta <- paste0(table1_res$variable,"_",table1_res$beta)
variable_dict <- data.frame(variable_beta=c("agedx_factor_beta1","stage_factor_beta1","er_factor_beta1",
                                            "onebef_bmi_factor_beta1","onebef_bmi_factor_beta2"),
                            variable_name=c("Age > 65 vs. Age <= 65",
                                            "Stage III/IV vs. Stage I/II",
                                            "ER+ vs. ER-",
                                            "Obesity vs. Normal",
                                            "Overweight vs. Normal"))

table1_res <- merge(table1_res, variable_dict, by="variable_beta")
table1_res$variable_name <- factor(table1_res$variable_name, levels=c("Age > 65 vs. Age <= 65",
                                                                      "Stage III/IV vs. Stage I/II",
                                                                      "ER+ vs. ER-",
                                                                      "Overweight vs. Normal",
                                                                      "Obesity vs. Normal"))

for (cell_type in sig_primary_res$cell_type_full_method_full){
  tmp <- table1_res[table1_res$cell_type==cell_type,]
  tmp$tissue <- ifelse(tmp$tissue=="tumor","Tumor","Normal")
  tmp$tissue <- factor(tmp$tissue, levels=c("Tumor","Normal"))
  cell_type_short <- gsub(" \\(Expression signature\\)","",cell_type)
  ggplot(data=tmp, aes(y=variable_name, x=effect, xmin=lower, xmax=upper, fill=tissue)) +
    geom_linerange(size=0.6,position=position_dodge(width = 0.5)) +
    geom_point(size=2.5, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
    labs(title=cell_type_short, x='Effect Size', y='') +
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values=c("darkred", "#466983")) +
    theme_bw()+
    theme(legend.title= element_blank(),
          plot.title = element_text(size = 10, face = "bold"))
  ggsave(paste0("results/plots/supp/covar_",cell_type,".pdf"), width=5, height=3)
}
