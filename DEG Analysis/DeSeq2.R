source('http://bioconductor.org/biocLite.R')
biocLite('DESeq2')
biocLite("edgeR")
biocLite("limma")
biocLite("statmod")
biocLite("affycoretools")
biocLite("ReportingTools")
library("edgeR")
library("statmod")
library('DESeq2')

countdata_1 <- read.delim("Orthoblast_RSEM_merged_matrix", header=TRUE, row.names="library") #read in the file of the count data and call it countdata, row.names tells the name in the top left cell, the gene names

countdata_1 = round(countdata_1)

keep <- rowSums(cpm(countdata_1)>1) >=4
# keeps rows (genes) where at least 4 columns (libraries) have at least 1 count per million. This means that if a gene is only expressed in say one treatment (which has three replicates), this gene will not be thrown out of the analysis

countdata<- countdata_1[keep,] #formatting for organizing the kept rows that summed to at least 1 cpm in the step above

countdata=round(countdata)

coldata <- read.delim("../column_2transcriptomes.txt", header=TRUE, row.names=1)

#select one from below:
countdata = countdata[, c(1,2,3,4,5,6,7,8,9,10,11,12)] #select on timepoint one (day 9)
countdata = countdata[, c(1,2,5,6,9,11)] #select timepoint one GOL only
countdata = countdata[, c(3,4,7,8,10,12)] #select timepoint one PAC only
countdata = countdata[, c(13,14,15,16,17,18,19,20,21,22,23,24)] #select on timepoint two (day 29)
countdata = countdata[, c(15,16,19,20,22,24)] #select timepoint two PAC only
countdata = countdata[, c(13,14,17,18,21,23)] #select timpoint two GOL only
countdata = countdata[, c(1,2,5,6,9,11,13,14,17,18,21,23)] #select all GOL (TP one and two, L and H)
countdata = countdata[, c(3,4,7,8,10,12,15,16,19,20,22,24)] #select all PAC (Tp1&2, H&L)



ddsFullCountTable <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design = ~ TREATMENT)

ddsFull <- DESeq(ddsFullCountTable) # this is the analysis!
head(ddsFull)
res=results(ddsFull)
res
resOrdered=res[order(res$padj),]
head(resOrdered)
sum(res$padj<0.05, na.rm=TRUE)


rld <-rlogTransformation(ddsFull)

write.table(assay(rld), "DEG_logtransformed_Orthoblast_matrix", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

head(assay(rld))

hist(assay(rld))

#install.packages("gplots")
library(gplots) 
#install.packages("RColorBrewer")
library("RColorBrewer") 
library( "genefilter" )


topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 40)



heatmap.2(assay(rld)[topVarGenes, ], scale="row",
          trace="none", dendrogram="both", key=TRUE, keysize = 1.5, margins =c(3,11), density.info = "density",
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = heatmap.2( assay(rld)[ topVarGenes, ], # category labels
                           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),  # color key
                           lty= 1,             # line style
                           lwd = 10))          # line width


write.table(resOrdered, "output_DEG_GOLtp2.txt")   
