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
scMB10x01_filtered <- readRDS("./20231211_scMB10x01_filtered.rds")
#####

#subset per cell type
#subset per cell type
scMB10x01_EC <- scMB10x01_filtered[,scMB10x01_filtered$celltypes == "Endothelial cells"]

scMB10x01_Pericytes <- scMB10x01_filtered[,scMB10x01_filtered$celltypes == "Pericytes"]
scMB10x01_Pericytes1 <- scMB10x01_filtered[,scMB10x01_filtered$label == 2]
scMB10x01_Pericytes2 <- scMB10x01_filtered[,scMB10x01_filtered$label == 3]

scMB10x01_Astrocytes <- scMB10x01_filtered[,scMB10x01_filtered$celltypes == "Astrocytes"]

sca_EC = SceToSingleCellAssay(scMB10x01_EC)
sca_Pericytes = SceToSingleCellAssay(scMB10x01_Pericytes)
sca_Pericytes1 = SceToSingleCellAssay(scMB10x01_Pericytes1)
sca_Pericytes2 = SceToSingleCellAssay(scMB10x01_Pericytes2)
sca_Astrocytes = SceToSingleCellAssay(scMB10x01_Astrocytes)

###############################
colData(sca_EC)$MULTIseq_ID_call2<-relevel(colData(sca_EC)$MULTIseq_ID_call2,"control_extract")
colData(sca_Pericytes)$MULTIseq_ID_call2<-relevel(colData(sca_Pericytes)$MULTIseq_ID_call2,"control_extract")
colData(sca_Pericytes1)$MULTIseq_ID_call2<-relevel(colData(sca_Pericytes1)$MULTIseq_ID_call2,"control_extract")
colData(sca_Pericytes2)$MULTIseq_ID_call2<-relevel(colData(sca_Pericytes2)$MULTIseq_ID_call2,"control_extract")
colData(sca_Astrocytes)$MULTIseq_ID_call2<-relevel(colData(sca_Astrocytes)$MULTIseq_ID_call2,"control_extract")

MAST_model_EC <- zlm( ~ MULTIseq_ID_call2, sca = sca_EC, exprs_value = 'logcounts')
MAST_model_Pericytes <- zlm( ~ MULTIseq_ID_call2, sca = sca_Pericytes, exprs_value = 'logcounts')
MAST_model_Pericytes1 <- zlm( ~ MULTIseq_ID_call2, sca = sca_Pericytes1, exprs_value = 'logcounts')
MAST_model_Pericytes2 <- zlm( ~ MULTIseq_ID_call2, sca = sca_Pericytes2, exprs_value = 'logcounts')
MAST_model_Astrocytes <- zlm( ~ MULTIseq_ID_call2, sca = sca_Astrocytes, exprs_value = 'logcounts')

#only test the condition coefficient.
summaryCond_EC <- summary(MAST_model_EC, doLRT='MULTIseq_ID_call2iRBC_extract') 
summaryCond_Pericyte <- summary(MAST_model_Pericytes, doLRT='MULTIseq_ID_call2iRBC_extract')
summaryCond_Pericyte1 <- summary(MAST_model_Pericytes1, doLRT='MULTIseq_ID_call2iRBC_extract')
summaryCond_Pericyte2 <- summary(MAST_model_Pericytes2, doLRT='MULTIseq_ID_call2iRBC_extract')
summaryCond_Astrocyte <- summary(MAST_model_Astrocytes, doLRT='MULTIseq_ID_call2iRBC_extract')

#save

saveRDS(summaryCond_EC, file = "20231211_MAST_EC.rds")
saveRDS(summaryCond_Pericyte, file = "20231211_MAST_PC.rds")
saveRDS(summaryCond_Pericyte1, file = "20231211_MAST_PC1.rds")
saveRDS(summaryCond_Pericyte2, file = "20231211_MAST_PC2.rds")
saveRDS(summaryCond_Astrocyte, file = "20231211_MAST_A.rds")