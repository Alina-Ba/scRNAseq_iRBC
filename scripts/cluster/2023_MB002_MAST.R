#library(scater)
#library(scran)
library(ggplot2)
library(patchwork)
library(deMULTIplex)
library(viridis)
library(BiocGenerics)
library(reshape2)
library(DESeq2)
library(glmGamPoi)
library(apeglm)
library(SingleCellExperiment)
library(scran)
library(MAST)

###############################
setwd("/g/bernabeu-hd/batzilla/scRNAseq_final/")
rescaled <- readRDS("./20240103_scMB10x02_integrated_BBBonly.rds")
#####
#filter for human genes
rescaled_human_all <- rescaled
rescaled_human <- rescaled_human_all[!grepl("^PFHG",rownames(rescaled_human_all)),]
rownames(rescaled_human) <- gsub("GRCh38-", "",rownames(rescaled_human))
#subset per cell type
rescaled_human_EC <- rescaled_human[,rescaled_human$celltypes == "Endothelial cells"]
rescaled_human_Pericytes <- rescaled_human[,rescaled_human$celltypes == "Pericytes"]
rescaled_human_Astrocytes <- rescaled_human[,rescaled_human$celltypes == "Astrocytes"]

rescaled_human_Pericytes1 <- rescaled_human[,rescaled_human$label == 2]
rescaled_human_Pericytes2 <- rescaled_human[,rescaled_human$label == 4]

sca_EC = SceToSingleCellAssay(rescaled_human_EC)
sca_Pericytes = SceToSingleCellAssay(rescaled_human_Pericytes)
sca_Astrocytes = SceToSingleCellAssay(rescaled_human_Astrocytes)

sca_Pericytes1 = SceToSingleCellAssay(rescaled_human_Pericytes1)
sca_Pericytes2 = SceToSingleCellAssay(rescaled_human_Pericytes2)

###############################
colData(sca_EC)$MULTIseq_ID_call<-relevel(colData(sca_EC)$MULTIseq_ID_call,"RBC")
colData(sca_Pericytes)$MULTIseq_ID_call<-relevel(colData(sca_Pericytes)$MULTIseq_ID_call,"RBC")
colData(sca_Astrocytes)$MULTIseq_ID_call<-relevel(colData(sca_Astrocytes)$MULTIseq_ID_call,"RBC")


MAST_model_EC <- zlm( ~ MULTIseq_ID_call+batch, sca = sca_EC, exprs_value = 'logcounts')
MAST_model_Pericytes <- zlm( ~ MULTIseq_ID_call+batch, sca = sca_Pericytes, exprs_value = 'logcounts')
MAST_model_Astrocytes <- zlm( ~ MULTIseq_ID_call+batch, sca = sca_Astrocytes, exprs_value = 'logcounts')

#only test the condition coefficient - Trophs
summaryCond_EC_Troph <- summary(MAST_model_EC, doLRT='MULTIseq_ID_callTrophs') 
summaryCond_Pericyte_Troph <- summary(MAST_model_Pericytes, doLRT='MULTIseq_ID_callTrophs')
summaryCond_Astrocyte_Troph <- summary(MAST_model_Astrocytes, doLRT='MULTIseq_ID_callTrophs')

#only test the condition coefficient - Schizonts
summaryCond_EC_Schizont <- summary(MAST_model_EC, doLRT='MULTIseq_ID_callSchizonts') 
summaryCond_Pericyte_Schizont <- summary(MAST_model_Pericytes, doLRT='MULTIseq_ID_callSchizonts')
summaryCond_Astrocyte_Schizont <- summary(MAST_model_Astrocytes, doLRT='MULTIseq_ID_callSchizonts')

#save
saveRDS(summaryCond_EC_Troph, file = "20240104_MAST_EC_Troph.rds")
saveRDS(summaryCond_Pericyte_Troph, file = "20240104_MAST_PC_Troph.rds")
saveRDS(summaryCond_Astrocyte_Troph, file = "20240104_MAST_A_Troph.rds")

saveRDS(summaryCond_EC_Schizont, file = "20240104_MAST_EC_Schiont.rds")
saveRDS(summaryCond_Pericyte_Schizont, file = "20240104_MAST_PC_Schiont.rds")
saveRDS(summaryCond_Astrocyte_Schizont, file = "20240104_MAST_A_Schiont.rds")