# Weighted correlation network analysis using WGCNA

**Data preparation**

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

#Use limma to remove batch effects / confounding co-variates
library(limma)
mm <- model.matrix(~ Sub_Group, colData(vsd))

#Use limma to remove batch effects - RNA_Batch and RIN
vsd_mat_batch <- limma::removeBatchEffect(vsd_mat, batch=vsd$Batch, batch2=vsd$RIN.cat, design=mm)

#correct for sequencing depth and batch
vsd_mat_batch <- limma::removeBatchEffect(vsd_mat_batch,batch=as.factor(colData$seqdep.cat), batch2=as.factor(colData$Sequencing_batch),design=mm)

#correct for sex and age
vsd_mat_batch <- limma::removeBatchEffect(vsd_mat_batch, batch=as.factor(colData$Sex), covariates=colData$AgeBL, design=mm)

#Export corrected vsd file
write.csv(vsd_mat_batch, file= "vsd_mat_batch.csv")
```
**Needed for WGCNA analysis**  

Matrix of expression levels - rows = samples, columns = genes  
Traits - rows = samples, columns = traits  

```
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
vsd_mat_batch <- read.csv("vsd_mat_batch.csv", header=TRUE, check.names = FALSE, row.names = 1)
colData <- read.csv("colData.csv", header=TRUE, check.names = FALSE, row.names = 1)

#transpose vsd_mat_batch
datExpr0 = as.data.frame(t(vsd_mat_batch_sex_age))

#Check of genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

#cluster samples to look for outliers
sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#remove outliers
abline(h = 110, col = "red")
clust = cutreeStatic(sampleTree, cutHeight = 110, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#cluster new set of data
sampleTree = hclust(dist(datExpr), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#create trait data table
allTraits <- colData[, c(3,14,17:24,53,57:66,68)]
samples <- rownames(datExpr)
datTraits <- allTraits[c(samples),]

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap", cex.dendroLabels = 0.5)
# save data expression and traits data
save(datExpr,datTraits, file="Project_10082_WGCNA_data_both__batch_sex_age_paper.RData")

# Choosing the soft-thresholding power: analysis of network topology
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#Constructing the gene network and identifying modules
net = blockwiseModules(datExpr, power = 3, networkType = "signed", 
                       TOMType = "signed", minModuleSize = 20,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3, maxBlockSize = 20000, corType = "bicor", maxPOutliers = 0.05)
#how many modules and what sizes
table(net$colors)

#dendrogram
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Project_10082_WGCNA_data_both_batch_sex_age_paper_SF3.RData")

# Relating modules to external information
# Load the expression and trait data saved in the first part (if needed)
lnames = load(file = "Project_10082_WGCNA_data_both__batch_sex_age_paper.RData")
# Load network data saved in the second part (if needed)
lnames = load(file = "networkConstruction_WGCNA_project_10082_both_batch_sex_paper.RData")
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#graphical representation
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#correlation restricted data_traits (no cell types) + colour-blind friendly colours
datTraits_2 <- datTraits[,1:11]
moduleTraitCor = cor(MEs, datTraits_2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#graphical representation
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits_2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab = 0.7,
               cex.lab.x = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#Gene relationship to trait and important modules: Gene Significance and Module Membership


# Define variable Sub_Group containing the Sub_Group column of datTrait
Sub_Group <- as.data.frame(datTraits$Sub_Group)
names(Sub_Group) = "Sub_Group"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, Sub_Group, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(Sub_Group), sep="")
names(GSPvalue) = paste("p.GS.", names(Sub_Group), sep="")

#Intramodular analysis: identifying genes with high GS and MM
module = "grey"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sub_Group",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

#load in annotation file
annot = read.csv(file="gene_annotation_gene_entrez.csv", head=TRUE)
dim(annot)
names(annot)

genes = names(datExpr)
genes2annot = match(genes,annot$gene_id)
# The following is the number or probes without annotation:
sum(is.na(genes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(gene_id = genes,
                       gene_name = annot$gene_name[genes2annot],
                       gene_type = annot$gene_type[genes2annot],
                       entrez_id = annot$entrez[genes2annot],
                       symbol = annot$symbol[genes2annot],
                       description = annot$description[genes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for Sub_Group
modOrder = order(-abs(cor(MEs, Sub_Group, use = "p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Sub_Group))
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo.csv")

#GOEnrichment
# Get the corresponding Locus Link IDs (entrez)
allLLIDs = annot$entrez[genes2annot]

#Enrichment within R
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "human", nBestP = 10)
tab = GOenr$bestPTerms[[4]]$enrichment
write.table(tab, file = "GO_enrichment_top_10_unsigned.txt", sep = "\t", quote = FALSE)

#Adjusting P-values for multiple sampling (FDR, BH)
subGroupPvalue <- moduleTraitPvalue[,2]
subGroupQvalue <- p.adjust(subGroupPvalue, method="BH")
sub_group_stats <- cbind(subGroupPvalue,subGroupQvalue)
write.csv(sub_group_stats, file = "sub_group_stats.csv")
HeightPvalue <- moduleTraitPvalue[,5]
HeightQvalue <- p.adjust(HeightPvalue, method="BH")
Height_stats <- cbind(HeightPvalue,HeightQvalue)
write.csv(Height_stats, file = "Height_stats.csv")
WeightPvalue <- moduleTraitPvalue[,6]
WeightQvalue <- p.adjust(WeightPvalue, method="BH")
Weight_stats <- cbind(WeightPvalue,WeightQvalue)
write.csv(Weight_stats, file = "Weight_stats.csv")
BMIPvalue <- moduleTraitPvalue[,7]
BMIQvalue <- p.adjust(BMIPvalue, method="BH")
BMI_stats <- cbind(BMIPvalue,BMIQvalue)
write.csv(BMI_stats, file = "BMI_stats.csv")
EducPvalue <- moduleTraitPvalue[,8]
EducQvalue <- p.adjust(EducPvalue, method="BH")
Educ_stats <- cbind(EducPvalue,EducQvalue)
write.csv(Educ_stats, file = "Educ_stats.csv")
SmokePvalue <- moduleTraitPvalue[,9]
SmokeQvalue <- p.adjust(SmokePvalue, method="BH")
Smoke_stats <- cbind(SmokePvalue,SmokeQvalue)
write.csv(Smoke_stats, file = "Smoke_stats.csv")

#Visualising the network of eigengenes
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate traits from the clinical traits
Sub_Group = as.data.frame(datTraits$Sub_Group);
names(Sub_Group) = "Sub_Group"

Height = as.data.frame(datTraits$Height);
names(Height) = "Height"

Weight = as.data.frame(datTraits$WeightKg);
names(Weight) = "Weight"

BMI = as.data.frame(datTraits$BMI);
names(BMI) = "BMI"

YearsEduc = as.data.frame(datTraits$YearsEduc);
names(YearsEduc) = "YearsEduc"

Smoking = as.data.frame(datTraits$SmokeHx);
names(Smoking) = "Smoking"

# Add the weight to existing module eigengenes
#MET = orderMEs(cbind(MEs, Sub_Group, Height, Weight, BMI, YearsEduc, Smoking))
MET = orderMEs(cbind(MEs, Sub_Group))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
```
