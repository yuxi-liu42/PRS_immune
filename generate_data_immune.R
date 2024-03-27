library(data.table)
library(devtools)
library(xCell)
library(EPIC)
library(quantiseqr)
library(MCPcounter)

##--------------------------------------
## CIBERSORTx 
##--------------------------------------

### Input

lm22 <- read_excel("data/immune/processed_input/cibersortx/LM22.xlsx")
dim(lm22)[1]
sum(lm22$`Gene symbol` %in% rownames(tumor.input.symbol))
sum(lm22$`Gene symbol` %in% rownames(tumor.input.symbol)) / dim(lm22)[1]

### Run 

# Go to https://cibersortx.stanford.edu/runcibersortx.php
# Select 2. Impute Cell Fractions
# Choose Signature matrix file: LM22
# Upload mixture files: processed_symbol_all_expr.tsv

# Do not enable batch correction
# Do not disable quantile normalization
# Run in both absolute and relative mode
# Permutation = 100

##--------------------------------------
## TIMER 2.0
##--------------------------------------

### Input

for (i in 1:round(dim(tumor.input.symbol)[2]/100)){
  temp <- tumor.input.symbol[,((i-1)*100+1):min((i*100),dim(tumor.input.symbol)[2])]
  write.csv(temp,paste0("data/immune/processed_input/timer2/processed_symbol_tumor_",i,".csv"), row.name=rownames(temp))
}

for (i in 1:round(dim(normal.input.symbol)[2]/100)){
  temp <- normal.input.symbol[,((i-1)*100+1):min((i*100),dim(normal.input.symbol)[2])]
  write.csv(temp,paste0("data/immune/processed_input/timer2/processed_symbol_normal_",i,".csv"), row.name=rownames(temp))
}

### Run 

# Go to http://timer.cistrome.org/
# Select Estimation
# Choose Human and BRCA
# Upload mixture files: processed_symbol_tumor_*.csv and processed_symbol_normal_*.csv

##-------------------
## xCell 
##-------------------

### Input

# matrix should contain HUGO gene symbols as row names and the columns are samples.
xcell.input <- all.input.symbol
row.names(xcell.input) <- xcell.input$GeneSymbol
xcell.input <- subset(xcell.input, select=-GeneSymbol)

### Run

xcell.res <- xCellAnalysis(xcell.input, rnaseq=F)
temp <- as.data.frame(t(xcell.res))
# write.csv(temp, "data/immune/raw_output/xcell/output_xcell_all.csv")

##-------------------
## EPIC 
##-------------------

### Input

epic.input <- all.input.symbol
row.names(epic.input) <- epic.input$GeneSymbol
epic.input <- subset(epic.input, select=-GeneSymbol)

### Run

epic.res <- EPIC(bulk = epic.input)
temp <- as.data.frame(epic.res$cellFractions)
# write.csv(temp, "data/immune/raw_output/epic/output_epic_all.csv")

##-------------------
## quanTIseq 
##-------------------

### Input

quantiseq.input <- all.input.symbol
row.names(quantiseq.input) <- quantiseq.input$GeneSymbol
quantiseq.input <- subset(quantiseq.input, select=-GeneSymbol)

### Run

?run_quantiseq

quantiseq.res <- run_quantiseq(
  expression_data = quantiseq.input,
  signature_matrix = "TIL10",
  is_arraydata = TRUE,
  is_tumordata = TRUE,
  scale_mRNA = TRUE
)

# write.csv(quantiseq.res, "data/immune/raw_output/quantiseq/output_quantiseq_all.csv", row.names = F)

##-------------------
## MCPcounter 
##-------------------

### Input

mcpcounter.input <- all.input.symbol
row.names(mcpcounter.input) <- mcpcounter.input$GeneSymbol
mcpcounter.input <- subset(mcpcounter.input, select=-GeneSymbol)

### Run

?MCPcounter.estimate

mcpcounter.res <- MCPcounter.estimate(
  expression = mcpcounter.input,
  featuresType = "HUGO_symbols"
)

temp <- as.data.frame(t(mcpcounter.res))
write.csv(temp, "data/immune/raw_output/mcpcounter/output_mcpcounter_all.csv")