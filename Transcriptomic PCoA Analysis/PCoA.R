source('http://bioconductor.org/biocLite.R')
biocLite('DESeq2')
biocLite("edgeR")
biocLite("limma")
biocLite("statmod")
biocLite("affycoretools")
biocLite("ReportingTools")
install.packages("pheatmap")
install.packages("vegan")
install.packages("ape")
install.packages("rgl")
library("edgeR")
library("statmod")
library('DESeq2')
library('vegan')
library('ape')

setwd("~/Desktop/DESeq2/Orthoblast_transcriptome")

countdata_1 <- read.delim("Orthoblast_RSEM_merged_matrix", header=TRUE, row.names="library") #read in the file of the count data and call it countdata, row.names tells the name in the top left cell, the gene names

countdata_1 = round(countdata_1)

keep <- rowSums(cpm(countdata_1)>1) >=4
# keeps rows (genes) where at least 4 columns (libraries) have at least 1 count per million. This means that if a gene is only expressed in say one treatment (which has three replicates), this gene will not be thrown out of the analysis

countdata<- countdata_1[keep,] #formatting for organizing the kept rows that summed to at least 1 cpm in the step above

countdata=round(countdata)

countdata = countdata[, c(1,2,3,4,5,6,7,8,9,10,11,12)] #select on timepoint one (day 9)
countdata = countdata[, c(13,14,15,16,17,18,19,20,21,22,23,24)] #select on timepoint two (day 29)
countdata = countdata[, c(13,14,17,18,21,23)] #select on timepoint two (day 29) and LU (GOL) only

coldata <- read.delim("../column_for_PCA_tp2.txt", header=TRUE, row.names=1) #use all_v3 for all data

ddsFullCountTable <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design = ~ POP*TREATMENT) #use POP*TREATMENT for tp 1 and 2, POP*TREATMENT*TP for all data

ddsFull <- DESeq(ddsFullCountTable) # this is the analysis!
head(ddsFull)
res=results(ddsFull)
res
resOrdered=res[order(res$padj),]
head(resOrdered)
sum(res$padj<0.05, na.rm=TRUE)


rld <-rlogTransformation(ddsFull)

# assembling table of conditions to lable PCoA plot:
# (in the chunk below, replace factor1 and factor2 with your actual factor names from myConditions table)
factor1=as.character(colData(ddsFull)$TREATTP) #select for all data
factor2=as.character(colData(ddsFull)$POP) #select for all data and tp1 and tp2
factor3=as.character(colData(ddsFull)$condition) #select for all data and tp 1 and tp2
factor4=as.character(colData(ddsFull)$TREATMENT) #select for all data and tp 1 and tp2
factor5=as.character(colData(ddsFull)$TP) #select for all data



# actual PCoA analysis
vsd = assay(rld)
dds.pcoa=pcoa(vegdist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors
percent <- dds.pcoa$values$Eigenvalues
percent / sum(percent) #percent for each axes

#plot PC axis 1 and 2 for all data
plot(scores[,1], scores[,2], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(0, 15, 1, 16)[as.numeric(as.factor(factor1))], xlab = "PC1 (42.47%)", ylab = "PC2 (17.22%)")
ordispider(scores,factor3, col=c("blue", "blue", "blue", "blue", "orange", "orange", "orange", "orange"))


#plot PC axis 3 and 4 for all data
plot(scores[,3], scores[,4], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(0, 15, 1, 16)[as.numeric(as.factor(factor1))], xlab = "PC1 (42.47%)", ylab = "PC2 (17.22%)")


#testing significance between group variances (need to do this because adonis will give sig p value even if groups overlap, but its sig b/c the variance is different)
#https://github.com/vegandevs/vegan/issues/233
input <- vegdist(vsd, method="manhattan")
mod <- betadisper(input, group = factor2, type = "median")
mod

oneByTwo=paste(factor1,factor2,sep=".")
conditions=data.frame(cbind(factor1,factor2,oneByTwo))
adonis2(t(vsd)~factor2*factor4*factor5, data = conditions, permutations = 1000000, method = "manhattan")

##test adonis order doesn't matter for my data:
adonis2(t(vsd)~factor4*factor5*factor2, data = conditions, permutations = 1000000, method = "manhattan") #almost exactly the same as above, don't need to worry about order!


#plot PC axis 1 and 2 for tp1
plot(scores[,1], scores[,2], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(0, 1)[as.numeric(as.factor(factor4))], xlab = "PC1 (31.17%)", ylab = "PC2 (14.65%)")
ordispider(scores,factor3, col=c("blue", "blue", "orange", "orange"))

oneByTwo=paste(factor2,factor4,sep=".")
conditions=data.frame(cbind(factor2,factor4,oneByTwo))
adonis2(t(vsd)~factor2*factor4, data = conditions, permutations = 1000000, method = "manhattan")

#plot PC axis 3 and 4 for tp1
plot(scores[,3], scores[,4], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(0, 1)[as.numeric(as.factor(factor4))], xlab = "PC1 (11.79%)", ylab = "PC2 (9.46%)")
ordispider(scores,factor3, col=c("blue", "blue", "orange", "orange"))


#plot PC1 and 2 for tp2
plot(scores[,1], scores[,2], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(15, 16)[as.numeric(as.factor(factor4))], xlab = "PC1 (27.01%)", ylab = "PC2 (15.27%)")
ordispider(scores,factor3, col=c("blue", "blue", "orange", "orange"))

oneByTwo=paste(factor2,factor4,sep=".")
conditions=data.frame(cbind(factor2,factor4,oneByTwo))
adonis2(t(vsd)~factor2+factor4, data = conditions, permutations = 1000000, method = "manhattan")
adonis2(t(vsd)~factor3, data = conditions, permutations = 1000000, method = "manhattan")

adonis2(t(vsd)~factor4, data = conditions, permutations = 1000000, method = "manhattan") ### for GOL tp2 only!

#plot PC axis 3 and 4 for tp2
plot(scores[,3], scores[,4], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(0, 1)[as.numeric(as.factor(factor4))], xlab = "PC1 (12.72%)", ylab = "PC2 (9.19%)")
ordispider(scores,factor3, col=c("blue", "blue", "orange", "orange"))


#plot PC1 and 2 for tp2 with phenotpyic data
PhenData=read.table("../../lipid_protein/lipid_protein_data_forPCA.txt", header = TRUE, row.names = "library") #import and attach phenotypic data
PhenColData <- read.delim("../../lipid_protein/lipid_protein_coldata_forPCA.txt", header=TRUE, row.names=1)

vare.pca <- rda(PhenData, scale = TRUE)
vare.pca
fit2 <- envfit(vare.pca, PhenColData, perm = 999)
fit2
plot(fit2)
plot(vare.pca, scaling = 3)
biplot(vare.pca, scaling = 3)
dev.off()

plot(scores[,1], scores[,2], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(15, 16)[as.numeric(as.factor(factor4))], xlab = "PC1 (27.01%)", ylab = "PC2 (15.27%)")
ordispider(scores,factor3, col=c("blue", "blue", "orange", "orange"))
plot(fit2)
plot(vare.pca, scaling = 3)

ord2 <- cca(PhenData ~ protein + lipid + phospholipid + sterol + fattyacid + triacylglycerol + waxester + respiration, PhenColData)
#plot(ord2, type="p")
fit <- envfit(ord2, coldata, perm = 999, display = "lc")
fit
plot(fit, col = "red")

