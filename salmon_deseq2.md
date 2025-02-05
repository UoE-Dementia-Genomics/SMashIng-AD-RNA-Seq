# RNA-seq analysis workflow using Salmon and DESeq2

Workflow using Salmon and DESeq2 to analyse short-read RNA-seq data looking at differential gene expression.

Workflow adapted from https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

Salmon documentation: https://salmon.readthedocs.io/en/latest/

DESeq2 documentation: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

**Install salmon in a conda environment**

```
conda create -n salmon salmon
conda activate salmon
conda install -c bioconda salmon
```

**Build Salmon Index using whole genome as a decoy**

https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

1: Download reference genome and transcriptome from Gencode (version 38 in this example)
```
wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
```


2: Extract name of genome targets

```
grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
```


3: Concatenate transcript and genome sequences (transcripts first)

`cat gencode.v38.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz`

4: Create salmon index

`salmon index -t gentrome.fa.gz -d decoys.txt -p 8 -i salmon_index --gencode -k 21`

In this example I am using kmer length (-k) of 21 because my reads were 50 bases. Use kmer length of 31 for reads of 75 bases or longer.
NOTE: --gencode flag is for removing extra metadata in the target header separated by | from the Gencode reference. You can skip it if using other references.

**Transcript quantification using Salmon (paired-end reads)**

`salmon quant -i salmon_index -l A -1 read1.fq.gz -2 read2.fq.gz -p 8 --validateMappings --gcBias --numGibbsSamples 20 -o output_directory`

see https://salmon.readthedocs.io/en/latest/salmon.html#description-of-important-options for description of parameters





Example slurm array script (runs through all read-pairs in current directory)

```
#!/bin/bash
#SBATCH -A Research_Project-XXXXX # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=0:20:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=user@exeter.ac.uk # email address
#SBATCH --array=0-14
 
samples=($(ls *r1*))
sample=${samples[${SLURM_ARRAY_TASK_ID}]}
stem=${sample::-20}
 
salmon quant -i salmon_index -l A -1 ${stem}_trimmed_51_r1.fq.gz -2 ${stem}_trimmed_51_r2.fq.gz -p 16 --validateMappings --gcBias --numGibbsSamples 20 -o ../salmon/$stem
```


One directory is produced per sample (containing quant.sf file amongst other data)

# R code for importing data and analysing using DESeq2
```
library(tximeta)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
```


**Importing data into R using tximeta**


```
#load file of sample data (one row per sample)
colData <- read.table("pheno_file.txt",stringsAsFactors=FALSE,sep="\t", header=TRUE,check.names = FALSE)

#convert co-variates to factors as needed e.g.
colData$Sex <- as.factor(colData$Sex)

#add columns for salmon import: name of folder for each sample and path to #“quant.sf” file
colData$names <- colData$RNASeq_sequencing_sample
colData$files <- file.path(colData$RNASeq_sequencing_sample,"quant.sf")
file.exists(colData$files)

#import salmon data
se <- tximeta(colData)
```


If the reference transcriptome checksum was recognized by tximeta, and if you have a working internet connection, tximeta will locate and download the relevant annotation data from various sources.

```
#summarise transcript-level quantifications to gene-level
gse <- summarizeToGene(se)

#optional – subset gse based on sample information in colData
#e.g. subset gse to include only visit 3
gse_visit_3 <- gse[,gse$Visit==3]
```



**Creating DESeq2 object and analysis**

```
#Create DESeq2 data object
#design is linear model – last expression is tested for differential expression
dds <- DESeqDataSet(gse, design = ~ Sex + Batch + Sub_Group)

#Count number of rows in dds
nrow(dds)

# Filter lowly expressed genes
# keep genes with at least 20 samples with a count of 10 or higher
keep <- rowSums(counts(dds) >= 10) >= 20
dds <- dds[keep,]
nrow(dds)

#optional steps when using blood samples and polyA library
#remove globin encoding genes
globin_ID <- grep("ENSG00000244734.*|ENSG00000213931.*|ENSG00000213934.*|ENSG00000196565.*|ENSG00000223609.*|ENSG00000206172.*|ENSG00000188536.*|ENSG00000130656.*|ENSG00000206177.*|ENSG00000086506.*", row.names(dds), value=TRUE)
dds <- dds[setdiff(rownames(dds), globin_ID),]
#verify they have been removed
globin_ID %in% rownames(dds)
nrow(dds)

# variance stabilised transformation
vsd <- vst(dds, blind = FALSE)

# clustering of samples + heatmap
library(pheatmap)
library(RColorBrewer)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#PCA plot of samples (coloured by sex, symbol based on subgroup)
library(ggplot2)
pcaData <- plotPCA(vsd, intgroup = c("Sex", "Sub_Group"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Sex, shape = Sub_Group)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
p + geom_text(label=p$data$name,cex=2,nudge_y = 1)

#Run differential expression pipeline
dds <- DESeq(dds)

#If you get error 
#20 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest
#Use code below to remove these rows

dds <- dds[which(mcols(dds)$betaConv),]
nrow(dds)
```
**Building results table**
```
res <- results(dds)
summary(res)

#out of 17122 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 50, 0.29%
#LFC < 0 (down)     : 32, 0.19%
#outliers [1]       : 0, 0%
#low counts [2]     : 3320, 19%
#(mean count < 20)


#Annotating and exporting result, in this case selecting results with padj < 0.05
library("AnnotationDbi")
library("org.Hs.eg.db")
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
resSig <- subset(res, padj < 0.05)
write.csv(as.data.frame(resSig),file="results.csv")

res

#log2 fold change (MLE): Sub Group 1 vs 0 
#Wald test p-value: Sub Group 1 vs 0 
#DataFrame with 17122 rows and 8 columns
#                      baseMean log2FoldChange     lfcSE       stat    pvalue      padj      symbol      entrez
#                     <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric> <character> <character>
#ENSG00000000003.15      6.9527    -0.14720357 0.2782935  -0.528951  0.596840  0.999999      TSPAN6        7105
#ENSG00000000419.14    230.2236     0.08877606 0.1524382   0.582374  0.560315  0.999999        DPM1        8813

#volcano plot of results
library(bioplotr)
plot_volcano(res)
```
**Using SVA with DESeq2**

Surrogate variable analysis identifies surrogate variables to adjust for unknown, unmodeled, or latent sources of noise. Surrogate variables can be used in linear models. 

sva package https://bioconductor.org/packages/release/bioc/html/sva.html 

The svaseq function is designed to be used with count data.

```
library(sva)

#start with a DESeq2 object with lowly expressed genes filtered out

#run DESeq command to estimate size factors
dds <- DESeq(dds)

#extract normalised gene counts
dat  <- counts(dds, normalized = TRUE)

#full model matrix containing variable of interest
mod  <- model.matrix(~ Sub_Group, colData(dds))

#null model matrix with intercept term
mod0 <- model.matrix(~ 1, colData(dds))

#estimate 2 surrogate variables 
#number defined by n.sv - leaving out this term means that svaseq determines the number of significant surrogate variables to calculate
svseq <- svaseq(dat, mod, mod0, n.sv = 2)

#create duplicate DESeq2 object
ddssva <- dds

#add surrogate variables to DESeq2 object
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]

#define design containing surrogate variables and variable of interest
design(ddssva) <- ~ SV1 + SV2 + Sub_Group

#run differential expression analysis
ddssva <- DESeq(ddssva)
```

Known adjustment variables e.g. sex, batch can also be used in the model matrices
e.g.

```
#full model matrix containing adjustment variables and variable of interest
mod  <- model.matrix(~ Sex + Batch + Sub_Group, colData(dds))

#null model matrix containing adjustment variables
mod0 <- model.matrix(~ Sex + Batch, colData(dds))

#estimate 2 surrogate variables
svseq <- svaseq(dat, mod, mod0, n.sv = 2)

#create duplicate DESeq2 object
ddssva <- dds

#add surrogate variables to DESeq2 object
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]

#define design containing adjustment variables, surrogate variables and variable of interest
design(ddssva) <- ~ Sex + Batch + SV1 + SV2 + Sub_Group

#run differential expression analysis
ddssva <- DESeq(ddssva)
```
