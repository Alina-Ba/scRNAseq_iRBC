library(BiocGenerics)
library(reshape2)
library(apeglm)
library(SingleCellExperiment)
library(scran)
library(MAST)
library(stats)

###############################
scMB10x01_filtered <- readRDS("./data/pre-processed_objects/20231211_scMB10x01_filtered.rds")
#####

#subset per cell type
#subset per cell type
scMB10x01_EC <- scMB10x01_filtered[,scMB10x01_filtered$celltypes == "Endothelial cells"]

scMB10x01_Pericytes <- scMB10x01_filtered[,scMB10x01_filtered$celltypes == "Pericytes"]

scMB10x01_Astrocytes <- scMB10x01_filtered[,scMB10x01_filtered$celltypes == "Astrocytes"]

sca_EC = SceToSingleCellAssay(scMB10x01_EC)
sca_Pericytes = SceToSingleCellAssay(scMB10x01_Pericytes)
sca_Astrocytes = SceToSingleCellAssay(scMB10x01_Astrocytes)

###############################
colData(sca_EC)$MULTIseq_ID_call2<-relevel(colData(sca_EC)$MULTIseq_ID_call2,"control_extract")
colData(sca_Pericytes)$MULTIseq_ID_call2<-relevel(colData(sca_Pericytes)$MULTIseq_ID_call2,"control_extract")
colData(sca_Astrocytes)$MULTIseq_ID_call2<-relevel(colData(sca_Astrocytes)$MULTIseq_ID_call2,"control_extract")

MAST_model_EC <- zlm( ~ MULTIseq_ID_call2, sca = sca_EC, exprs_value = 'logcounts')
MAST_model_Pericytes <- zlm( ~ MULTIseq_ID_call2, sca = sca_Pericytes, exprs_value = 'logcounts')
MAST_model_Astrocytes <- zlm( ~ MULTIseq_ID_call2, sca = sca_Astrocytes, exprs_value = 'logcounts')

#only test the condition coefficient.
summaryCond_EC <- summary(MAST_model_EC, doLRT='MULTIseq_ID_call2iRBC_extract') 
summaryCond_Pericyte <- summary(MAST_model_Pericytes, doLRT='MULTIseq_ID_call2iRBC_extract')
summaryCond_Astrocyte <- summary(MAST_model_Astrocytes, doLRT='MULTIseq_ID_call2iRBC_extract')

#save

saveRDS(summaryCond_EC, file = "./data/MAST_DE_results/20231211_MAST_EC.rds")
saveRDS(summaryCond_Pericyte, file = "./data/MAST_DE_results/20231211_MAST_PC.rds")
saveRDS(summaryCond_Astrocyte, file = "./data/MAST_DE_results/20231211_MAST_A.rds")