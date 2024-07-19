# Weighted correlation network analysis using WGCNA

Data preparation:  
1: Create DESeq2 object  
2: Keep genes with at least 20 samples with a count of 10 or higher  
3: Remove globin-encoding genes  
4: Variance stabilised transformation  
5: Correct for batch-effects / co-variates using limma 

```
#create DESeq2 object
dds <- DESeqDataSet(gse, design = ~ Sub_Group)
using counts and average transcript lengths from tximeta

#Count number of rows in dds
nrow(dds)

#Filter lowly expressed genes
keep <- rowSums(counts(dds) >= 10) >= 20  
dds <- dds[keep,]
nrow(dds)

#Find gene IDs of globin genes
globin_ID <- grep("ENSG00000244734.*|ENSG00000213931.*|ENSG00000213934.*|ENSG00000196565.*|ENSG00000223609.*|ENSG00000206172.*|ENSG00000188536.*|ENSG00000130656.*|ENSG00000206177.*|ENSG00000086506.*", row.names(dds), value=TRUE) 

#Remove globin genes from dataset
dds <- dds[setdiff(rownames(dds), globin_ID),]

#verify they have been removed
globin_ID %in% rownames(dds)
nrow(dds)
dds <- DESeq(dds)

#variance stabilised transformation
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

#Export file
write.csv(vsd_mat, file= "vsd_mat_WGCNA_project_10082_v2.csv")![image](https://github.com/user-attachments/assets/7b355788-54e9-47f9-bc9c-f55e1f210a12)

#Use limma to remove batch effects / confounding co-variates
library(limma)
mm <- model.matrix(~ Sub_Group, colData(vsd))

#Use limma to remove batch effects - RNA_Batch and RIN
vsd_mat_batch <- limma::removeBatchEffect(vsd_mat, batch=vsd$Batch, batch2=vsd$RIN.cat, design=mm)

#correct for sequencing depth and batch
vsd_mat_batch <- limma::removeBatchEffect(vsd_mat_batch,batch=as.factor(colData$seqdep.cat), batch2=as.factor(colData$Sequencing_batch),design=mm)

#correct for sex and age
vsd_mat_batch <- limma::removeBatchEffect(vsd_mat_batch, batch=as.factor(colData$Sex), covariates=colData$AgeBL, design=mm)![image](https://github.com/user-attachments/assets/df7366f2-354f-41f5-8ad9-3ab8d25f43b2)
```
