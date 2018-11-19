#WGCNA_scripts courtesy of Kevin Johnson (2018)

###Load the counts and generate normalized count files for WGCNA input. These should be filtered so that 90% of samples or so have at least 10 cpm.

library(limma) 
library(edgeR)
library(Rsubread)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(dynOmics)
library(cluster)
library(VennDiagram)
library(HTSFilter)
library(flashClust)
library(lubridate)
library(ggplot2)
library(MASS)
library(vegan)
library(rgl)
library(ape)
library(beepr)
library(GO.db)

setwd("~/Desktop/WGCNA/")

##load data for tp2 only
load("~/Desktop/Expression_phys_correlation/PCA_LabMeeting_Data.RData")
PhenData$pop <- c(1,1,0,0,1,1,0,0,1,0,1,0) #G is 1 and P is 0
PhenData$treatment <- c(1,1,1,1,0,0,0,0,1,1,0,0) #L is 1 and H is 0
PhenData2 <- PhenData[, c(3,4,5,6,7,8,9,10,11,12,13)]
PhenData2 = PhenData2[c(1,2,9,3,4,10,5,6,11,7,8,12),]

data_input = read.table("Orthoblast_RSEM_merged_matrix", header=T, row.names="library")
data_input <- round(data_input)
data_input = data_input[rowSums(cpm(data_input)>=5) >= 6,] #reduce to only genes with high expression some people restrict to 10.
data_input = data_input[, c(13,14,15,16,17,18,19,20,21,22,23,24)] #select on timepoint two (day 29)
data_input = data_input[, c(1,2,9,3,4,10,5,6,11,7,8,12)] #reorder to all the groups are clustered together--this will make the heatmaps look better)]

#l###load data for tp1 and tp2--all data
data_input = read.table("Orthoblast_RSEM_merged_matrix", header=T, row.names="library")
data_input <- round(data_input)
data_input = data_input[rowSums(cpm(data_input)>=5) >= 6,] #reduce to only genes with high expression some people restrict to 10.
data_input = data_input[, c(1,2,9,3,4,10,5,6,11,7,8,12,13,14,21,15,16,22,17,18,23,19,20,24)] #reorder to all the groups are clustered together--this will make the heatmaps look better)]

PhenData2 = read.delim("Traits_WGCNA.txt")


dim(data_input)
#traitData = read.table("WGCNA_traits.txt",sep="\t",stringsAsFactors = F,header=T);
#head(traitData)

#data_input <- rld
#load("PCA_LabMeeting_Data.RData")
traitData <- PhenData2

#Group <- factor(paste(traitData$Treatment,traitData$Temperature,sep="_"))
#cbind(traitData,Group=Group)
#Group
Group <- traitData$Condition
Group

y <- DGEList(counts=data_input, group=Group)
y <- calcNormFactors(y)
WGCNA_input <- cpm(y,log = T,prior.count = 2)#THis gives normalized log transformed counts that will be used in WGCNA
dim(WGCNA_input)
# Load the WGCNA package
source("http://bioconductor.org/biocLite.R")
biocLite("GO.db")
biocLite("impute")
biocLite("preprocessCore")
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

allowWGCNAThreads(nThreads=2)# This can allow WGNCA to use more threads depending on how many CPUs you have access to I use 3 on my personal computer

head(WGCNA_input)
dim(WGCNA_input);
names(WGCNA_input) <- row.names(WGCNA_input)

datExpr = as.data.frame(t(WGCNA_input));
#names of datExpr should print all gene names/IDs
names(datExpr)
#rownames of datExpr should print all sample names
rownames(datExpr)

gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


#check traitData info again
dim(traitData)
names(traitData)

# remove columns that hold information we do not need if any...I use all mine
allTraits = traitData
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the physiological traits.

Belegans_samples = rownames(WGCNA_input);
traitRows = match(Belegans_samples, allTraits$Sample);

#datTraits = allTraits[traitRows, -1];
#rownames(datTraits) = allTraits[traitRows, 1];

#
datTraits <- allTraits[, c(2,3,4,5,6,7,8,9,10,11)] #for tp2 data
datTraits <- allTraits[, c(2,3,4)] #for all data

collectGarbage();
 # Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")



save(datExpr, datTraits, file = "Belegans_ALL_01-dataInput.RData")

# Load the data saved in the first part if starting from here
lnames = load(file = "Belegans_ALL_01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
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

##Usually it is recomended to choose the lowest softPower that has an R2 at 90%
softPower = 7;
adjacency = adjacency(datExpr, power = softPower);


# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

beep(0)
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# Adjust the minimum module size, 30 is common as it is what is suggested in the official tutorials.
minModuleSize = 30;
# Module identification using dynamic tree cut:
#setting cutHeith =.9 means the maximum dissimilarity is (90%) that qualifies modules for merging.
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# set a threshold for merging modules
dynamicMergeCut(18) # calculates the threshold for module merging using the number of samples (in this case, 18)
MEDissThres = 0.3009439
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "Belegans_ALL-02-networkConstruction-stepByStep.RData")

load("Belegans_AM-02-networkConstruction-stepByStep.RData") #load("Belegans_ALL-02-networkConstruction-stepByStep.RData")
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p"); #####*********
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

moduleTraitCorAndPvalue=merge(moduleTraitCor, moduleTraitPvalue, by="row.names")

View(moduleTraitPvalue)
View(moduleTraitCorAndPvalue)

#tp2 only below
moduleTraitPvalue_sig_0.6 <- moduleTraitPvalue[c("MEbisque4","MEblue2","MEbrown4","MEchocolate2","MEcyan","MEdarkgoldenrod3","MEdarkolivegreen4","MEdarkseagreen2","MEdarkseagreen3","MEgreen4","MEhoneydew","MEmediumpurple4","MEmidnightblue","MEburlywood","MEdarkviolet","MElavenderblush3","MEmediumorchid3","MEmediumpurple","MEmistyrose4","MEorange2","MEorange4","MErosybrown4","MEroyalblue","MEviolet"),] #keeping these rows that have a significant pvalue and a correlation above 0.7
moduleTraitCor_sig_0.6 <- moduleTraitCor[c("MEbisque4","MEblue2","MEbrown4","MEchocolate2","MEcyan","MEdarkgoldenrod3","MEdarkolivegreen4","MEdarkseagreen2","MEdarkseagreen3","MEgreen4","MEhoneydew","MEmediumpurple4","MEmidnightblue","MEburlywood","MEdarkviolet","MElavenderblush3","MEmediumorchid3","MEmediumpurple","MEmistyrose4","MEorange2","MEorange4","MErosybrown4","MEroyalblue","MEviolet"),] #keeping these rows that have a significant pvalue and a correlation above 0.7

moduleTraitPvalue_sig_0.5 <- moduleTraitPvalue[c("MEbisque4","MEblue2","MEbrown4","MEchocolate2","MEcyan","MEdarkgoldenrod3","MEdarkolivegreen4","MEdarkseagreen2","MEdarkseagreen3","MEgreen4","MEhoneydew","MEmediumpurple4","MEmidnightblue","MEburlywood","MEdarkviolet","MElavenderblush3","MEmediumorchid3","MEmediumpurple","MEmistyrose4","MEorange2","MEorange4","MErosybrown4","MEroyalblue","MEviolet","MEcoral2","MEgreen1","MElightblue2","MEpink3"),] #keeping these rows that have a significant pvalue and a correlation above 0.7
moduleTraitCor_sig_0.5 <- moduleTraitCor[c("MEbisque4","MEblue2","MEbrown4","MEchocolate2","MEcyan","MEdarkgoldenrod3","MEdarkolivegreen4","MEdarkseagreen2","MEdarkseagreen3","MEgreen4","MEhoneydew","MEmediumpurple4","MEmidnightblue","MEburlywood","MEdarkviolet","MElavenderblush3","MEmediumorchid3","MEmediumpurple","MEmistyrose4","MEorange2","MEorange4","MErosybrown4","MEroyalblue","MEviolet","MEcoral2","MEgreen1","MElightblue2","MEpink3"),] #keeping these rows that have a significant pvalue and a correlation above 0.7

View(moduleTraitPvalue_sig_0.6)
View(moduleTraitCorAndPvalue_sig_0.6)


moduleTraitPvalue_sig2 <- moduleTraitPvalue[c(4,8,12,14,18,19,22,24,25,33,35,49,50),] #keeping these rows that have a significant pvalue and a correlation above 0.7
moduleTraitCor_sig2 <- moduleTraitCor[c(4,8,12,14,18,19,22,24,25,33,35,49,50),] #keeping these rows that have a significant pvalue and a correlation above 0.7
moduleTraitPvalue_sig2 <- moduleTraitPvalue[c("MEbisque4","MEblue2","MEbrown4","MEchocolate2","MEcyan","MEdarkgoldenrod3","MEdarkolivegreen4","MEdarkseagreen2","MEdarkseagreen3","MEgreen4","MEhoneydew","MEmediumpurple4","MEmidnightblue"),] #keeping these rows that have a significant pvalue and a correlation above 0.7
moduleTraitCor_sig2 <- moduleTraitCor[c("MEbisque4","MEblue2","MEbrown4","MEchocolate2","MEcyan","MEdarkgoldenrod3","MEdarkolivegreen4","MEdarkseagreen2","MEdarkseagreen3","MEgreen4","MEhoneydew","MEmediumpurple4","MEmidnightblue"),] #keeping these rows that have a significant pvalue and a correlation above 0.7
View(moduleTraitCor_sig2)

#######This is all the modules
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
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
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#######This displays all modules with significant pvalue
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor_sig, 2), "\n(",
                    signif(moduleTraitPvalue_sig, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor_sig)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor_sig,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#######This displays all modules with significant pvalue and correlation above 0.7
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor_sig_0.5, 2), "\n(",
                    signif(moduleTraitPvalue_sig_0.5, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor_sig_0.5)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor_sig_0.5,
               xLabels = names(datTraits),
               yLabels = rownames(moduleTraitCor_sig_0.5),
               ySymbols = rownames(moduleTraitCor_sig_0.5),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# Define variable weight containing the weight column of datTrait
#pCO2 = as.data.frame(datTraits$sterol);
#names(pCO2) = "pCO2 (uatm)"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

#geneTraitSignificance = as.data.frame(cor(datExpr, pCO2, use = "p"));
#GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

#names(geneTraitSignificance) = paste("GS.", names(pCO2), sep="");
#names(GSPvalue) = paste("p.GS.", names(pCO2), sep="");


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


module = "MElightblue2"
column = match(module, modNames);
moduleGenes = moduleColors==module;


sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for pCO2 (uatm)",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

names(datExpr)

##select one from below, change line 317 above to pick a different one below:
module_paleovioletred2 <- names(datExpr)[moduleColors=="palevioletred2"]
module_paleovioletred2
write.csv(module_paleovioletred2, "Module_paleovioletred2.csv")

################
# eigengene-heatmap plot (sanity check - is the whole module driven by just one crazy sample?)
# note: this part does not make much sense for unsigned modules


which.module="royalblue" 
datME=MEs

#quartz()
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,clabels=T,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")

length(datExpr[1,moduleColors==which.module ]) # number of genes in chosen module

#################
# saving selected modules for GO and KOG analysis (two-parts: Fisher test, MWU test within-module)

# calculating modul memberships for all genes for all modules
allkME =as.data.frame(signedKME(datExpr, MEs)) 
names(allkME)=gsub("kME","",names(allkME))
vsd.wg = t(datExpr)

whichModule="grey"
table(moduleColors==whichModule) # how many genes are in it?

# Saving data for Fisher-MWU combo test (GO_MWU)
inModuleBinary=as.numeric(moduleColors==whichModule)
write.csv(inModuleBinary, file=paste("inModuleBinary"))
combo=data.frame("gene"=row.names(vsd.wg),"Fish_kME"=inModuleBinary)  #edited this line, see original in Matz lab script
write.csv(combo,file=paste(whichModule,"-ALL.csv",sep=""),row.names=F,quote=F)

####################
annot = read.delim("~/Dropbox/BOX/Limacina_GE/WGCNA/Limacina_transcriptome_seqeunce_name.txt", row.names=1)
annot = Limacina_transcriptome_seqeunce_name
head(annot)
dim(annot)
annot <- annot[1:2]
names(annot)

annot$Sequence <- row.names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$Sequence)

# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(Limacina = probes,
                       geneDescription = annot$Description[probes2annot],
                       EnzymeID = annot$Enzyme[probes2annot],
                       GOterms = annot$GO[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for pCO2
modOrder = order(-abs(cor(MEs, pCO2, use = "p")));
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
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.pCO2..uatm.));
geneInfo = geneInfo0[geneOrder, ]
View(geneInfo)

write.csv(geneInfo, file = "geneInfo.csv")


###PLOTS
####
###
###The following are codes to make various heatmaps, TOMPlots, and top gene views### use as needed


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate temperature from the clinical traits
Respiration = as.data.frame(datTraits$mean_respiration);
names(Respiration) = "Respiration"

Temperature = as.data.frame(datTraits$Temperature);
names(Temperature) = "Temp"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs))
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

#check multi-dimensional scaling plot using two dimensions. The calculation may take some time
cmd1=cmdscale(as.dist(dissTOM),2)
sizeGrWindow(7, 6)
par(mfrow=c(1,1))
plot(cmd1, col=as.character(colorh1), main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

###Now restrict the visualization to the 30 most significant genes identified by network screening 
power=20 #use whatever you used the first time up above
color1=colorDynamicTOM
restGenes= (color1 != "grey")
diss1=1-TOMsimilarityFromExpr( datExpr[, restGenes], power = power )
hier1=hclust(as.dist(diss1), method="average" )
diag(diss1) = NA;
sizeGrWindow(7,7)
TOMplot(diss1^4, hier1, as.character(color1[restGenes]),
        main = "TOM heatmap plot, module genes" )
#In the heatmap, rows and columns correspond to single genes, 
#light colors represent low topological overlap, and 
#progressively darker orange and red colors represent higher topological overlap.

power=20
color1=colorDynamicTOM
restGenes= (color1 != "grey")
diss1=1-adjacency( datExpr[, restGenes], power = 20 )
hier1=hclust(as.dist(diss1), method="average" )
diag(diss1) = NA;
sizeGrWindow(7,7)
TOMplot(diss1^4, hier1, as.character(color1[restGenes]),
        main = "Adjacency heatmap plot, module genes" )

# In the heatmap, rows and columns correspond to single genes, light
# colors represent low adjacency, and progressively darker orange and 
# red colors represent higher adjacency. The corresponding gene dendrograms 
# and module assignment are shown on the left and top.


sizeGrWindow(7,7)
topList=rank(NS1$PCO2,ties.method="first")<=50#This sets how many genes to restrict to
gene.names= names(datExpr)[topList]
# The following shows the correlations between the 50 top genes
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="signed", useTOM=FALSE,
                   power=1, main="signed correlations")


sizeGrWindow(7,7)
# The following shows the correlations between the top genes
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="unsigned", useTOM=FALSE,
                   power=1, main="signed correlations")


sizeGrWindow(7,7)
# The following shows the TOM heatmap in a signed network
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="signed", useTOM=TRUE,
                   power=12, main="C. TOM in a signed network")
# The following shows the TOM heatmap in a unsigned network
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="unsigned", useTOM=TRUE,
                   power=6, main="D. TOM in an unsigned network")

