## INSTALL NECESSARY PACKAGES #######
source('http://bioconductor.org/biocLite.R')
biocLite('DESeq2')
install.packages("adegenet")
install.packages("vegan")
install.packages("pca3d")
install.packages("rgl")
library("rgl")
library("adegenet")
library("DESeq2")
library("vegan")
library("pca3d")

setwd("~/Desktop/Expression_phys_correlation") 


##### READ IN DATA ######
load("PCA_LabMeeting_Data.RData")
View(PhenData)
View(rld)


######## Determining how many PCs there are total
vsd = assay(rld)
dds.pcoa=pcoa(vegdist(t(rld),method="manhattan")/1000)
dds.pcoa
scores=dds.pcoa$vectors
scores #I think there are 11 PCs?
percent <- dds.pcoa$values$Eigenvalues
percent / sum(percent) #percent for each axes

xval <- xvalDapc(t(rld), PhenData$Pop, n.pca.max = 20, training.set = 0.9, #made n.pca.max 20 here, just a few PCs more than results above
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval
xval$DAPC #this should be the same output as dp2 below

######## establishing discriminant function for Population (Low upwelling site vs high upwelling site)
dp2=dapc(t(rld),PhenData$Pop, n.pca=3, n.da=1)
dp2
pred2=predict.dapc(dp2,newdata=(t(rld))) 

plot(density(pred2$ind.scores),col="black",bty="n",yaxt="n",ylab="",xlab="discriminant function",ylim=c(0,1.1),main="",mgp=c(2.3,1,0))

polygon(density(pred2$ind.scores[PhenData$Condition=="GH"]),col=rgb(225/225,128/225,0/225,alpha = 0.5),border="orange")
polygon(density(pred2$ind.scores[PhenData$Condition=="GL"]),col=rgb(252/255,228/255,10/255,alpha=0.1),border="orange")
polygon(density(pred2$ind.scores[PhenData$Condition=="PH"]),col=rgb(0,0,1,alpha=0.5),border="blue")
polygon(density(pred2$ind.scores[PhenData$Condition=="PL"]),col=rgb(10/255,195/255,252/255,alpha=0.1),border="blue")
#ran the above 4 lines twice to get darker colors on plots

######## establishing discriminant function for pH treatment (low vs high)
xval <- xvalDapc(t(rld), PhenData$Treatment, n.pca.max = 20, training.set = 0.9, #made n.pca.max 20 here, just a few PCs more than results above
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval

dp3=dapc(t(rld),PhenData$Treatment,n.pca=6, n.da=1)
dp3
pred3=predict.dapc(dp3,newdata=(t(rld))) 

plot(density(pred3$ind.scores),col="black",bty="n",yaxt="n",ylab="",xlab="discriminant function",ylim=c(0,0.8),main="",mgp=c(2.3,1,0))

polygon(density(pred3$ind.scores[PhenData$Condition=="GH"]),col=rgb(225/225,128/225,0/225,alpha = 0.5),border="orange")
polygon(density(pred3$ind.scores[PhenData$Condition=="GL"]),col=rgb(252/255,228/255,10/255,alpha=0.1),border="orange")
polygon(density(pred3$ind.scores[PhenData$Condition=="PH"]),col=rgb(0,0,1,alpha=0.5),border="blue")
polygon(density(pred3$ind.scores[PhenData$Condition=="PL"]),col=rgb(10/255,195/255,252/255,alpha=0.1),border="blue")

########### CONNECT DAPC RESULTS WITH PHENOTYPIC DATA ############
pop.dapc=data.frame(dp2$ind.coord)
treat.dapc=data.frame(dp3$ind.coord)
rownames(pop.dapc) #check rownames are correct
rownames(treat.dapc)

all.dapc=merge(pop.dapc, treat.dapc, by="row.names")
View(all.dapc)
row.names(all.dapc) <- all.dapc$Row.names
rownames(all.dapc)
all.dapc$Row.names = NULL #this tells it to ignore row.names as data
View(all.dapc)


#ADD DAPC RESULTS TO TRAIT DATA
m=merge(PhenData, all.dapc, by="row.names")
View(m)
row.names(m) <- m$Row.names
m$Row.names = NULL #this tells it to ignore row.names as data
View(m)
m$GE_pop <- m$LD1.x
m$GE_pH <- m$LD1.y

m.pca <- m[, c(4,5,6,7,8,9,10,11,14,15)] #Just keeping columns with numbers (not variables, ie pop or treatment) and getting rid of duplicate columns LD1.x and LD1.y
vare.pca <- rda(m.pca, scale = TRUE) #this calculates PCA axis for the phenotypic data
vare.pca #view it quickly
fit2 <- envfit(vare.pca ~m.pca$waxester + m.pca$GE_pop + m.pca$GE_pH + m.pca$triacylglycerol + m.pca$fattyacid + m.pca$sterol + m.pca$phospholipid + m.pca$lipid + m.pca$protein + m.pca$respiration, PhenData, perm = 999) #run envfit to get a p-value to see if there are any patterns in data
fit2 #view p-values
scores(vare.pca)
percent <- vare.pca$CA$eig
percent
percent / sum(percent)



biplot(vare.pca, scaling=3, type = "p", xlab = "PC1 (44%)", ylab = "PC2 (19%)") #add vectors to plot
points(vare.pca, scaling= 3, col=c("orange", "orange", "orange", "orange", "orange", "orange", "blue", "blue", "blue", "blue", "blue", "blue"), pch=c(15,16,15,16,15,16,15,16,15,16,15,16))
text(vare.pca, ordiArrowTextXY(bp))
dev.off()

biplot(vare.pca, scaling=3) #run this to get labels for arrows/vectors
ggplot(vare.pca)










#### trying it with contrasts between GL and GH ### this looks weird, but you can't do it with just half the data because then there is missing data in PCA and it won't work
"""
rld.GOL=rld[,PhenData$Pop=="G"]
tr.GOL=PhenData[PhenData$Pop=="G",]

dp4=dapc(t(rld.GOL),tr.GOL$Treatment,perc.pca=80, n.da=1)
dp4
pred4=predict.dapc(dp4,newdata=(t(rld))) 

plot(density(pred4$ind.scores),col="black",bty="n",yaxt="n",ylab="",xlab="discriminant function",ylim=c(0,0.8),main="",mgp=c(2.3,1,0))

polygon(density(pred4$ind.scores[PhenData$Condition=="GL"]),col=rgb(0,0,1,alpha=0.5),border="blue")
polygon(density(pred4$ind.scores[PhenData$Condition=="GH"]),col=rgb(10/255,195/255,252/255,alpha=0.5),border="skyblue")
polygon(density(pred4$ind.scores[PhenData$Condition=="PL"]),col=rgb(1,0,0,alpha=0.5),border="red")
polygon(density(pred4$ind.scores[PhenData$Condition=="PH"]),col=rgb(252/255,228/255,10/255,alpha=0.5),border="gold")

#### trying it with contrasts between PL and PH ### this looks weird, but you can't do it with just half the data because then there is missing data in PCA and it won't work

rld.PAC=rld[,PhenData$Pop=="P"]
tr.PAC=PhenData[PhenData$Pop=="P",]

dp5=dapc(t(rld.PAC),tr.PAC$Treatment,perc.pca=80, n.da=1)
dp5
pred5=predict.dapc(dp5,newdata=(t(rld))) 

plot(density(pred5$ind.scores),col="black",bty="n",yaxt="n",ylab="",xlab="discriminant function",ylim=c(0,0.8),main="",mgp=c(2.3,1,0))

polygon(density(pred5$ind.scores[PhenData$Condition=="GL"]),col=rgb(0,0,1,alpha=0.5),border="blue")
polygon(density(pred5$ind.scores[PhenData$Condition=="GH"]),col=rgb(10/255,195/255,252/255,alpha=0.5),border="skyblue")
polygon(density(pred5$ind.scores[PhenData$Condition=="PL"]),col=rgb(1,0,0,alpha=0.5),border="red")
polygon(density(pred5$ind.scores[PhenData$Condition=="PH"]),col=rgb(252/255,228/255,10/255,alpha=0.5),border="gold")

########### CONNECT DAPC RESULTS WITH PHENOTYPIC DATA ############
pop.dapc=data.frame(pred2$ind.scores)
colnames(pop.dapc)[1] <- "Pop"
treat.dapc=data.frame(pred4$ind.scores)
colnames(treat.dapc)[1] <- "pH"
popG.dapc=data.frame(pred4$ind.scores)
colnames(popG.dapc)[1] <- "popG"
popP.dapc=data.frame(pred5$ind.scores)
colnames(popP.dapc)[1] <- "PopP"
rownames(pop.dapc) #check rownames are correct
rownames(treat.dapc)
rownames(popG.dapc)
rownames(popP.dapc)

temp.dapc=merge(pop.dapc, popG.dapc, by="row.names")
row.names(temp.dapc) <- temp.dapc$Row.names
rownames(temp.dapc)
temp.dapc$Row.names = NULL #this tells it to ignore row.names as data
View(temp.dapc)

temp2.dapc=merge(temp.dapc, popP.dapc, by="row.names")
row.names(temp2.dapc) <- temp2.dapc$Row.names
rownames(temp2.dapc)
temp2.dapc$Row.names = NULL #this tells it to ignore row.names as data
View(temp2.dapc)

all.dapc=merge(temp2.dapc, treat.dapc, by="row.names")
View(all.dapc)
row.names(all.dapc) <- all.dapc$Row.names
rownames(all.dapc)
all.dapc$Row.names = NULL #this tells it to ignore row.names as data
View(all.dapc)

#ADD DAPC RESULTS TO TRAIT DATA
m=merge(PhenData, all.dapc, by="row.names")
View(m)
row.names(m) <- m$Row.names
m$Row.names = NULL #this tells it to ignore row.names as data
View(m)

m.pca <- m[, c(4,5,6,7,8,9,10,11,12,13,14,15)] #Just keeping columns with numbers (not variables, ie pop or treatment)
vare.pca <- rda(m.pca, scale = TRUE) #this calculates PCA axis for the phenotypic data
vare.pca #view it quickly
fit2 <- envfit(vare.pca ~m.pca$waxester + m.pca$Pop.y + m.pca$popG + m.pca$PopP + m.pca$pH + m.pca$triacylglycerol + m.pca$fattyacid + m.pca$sterol + m.pca$phospholipid + m.pca$lipid + m.pca$protein + m.pca$respiration, PhenData, perm = 999) #run envfit to get a p-value to see if there are any patterns in data
fit2 #view p-values
plot(vare.pca, scaling = 3) #plot phenotpyic and gene expression PCA 
biplot(vare.pca, scaling = 3) #add vectors to plot

### trying again with GH and GL but without adding back in all the new data, just GH and GL data ####
rld.GOL=rld[,PhenData$Pop=="G"]
tr.GOL=PhenData[PhenData$Pop=="G",]

dp5=dapc(t(rld.GOL),tr.GOL$Treatment,perc.pca=80, n.da=1)
dp5
pred5=predict.dapc(dp4,newdata=(t(rld.GOL))) 

plot(density(pred5$ind.scores),col="black",bty="n",yaxt="n",ylab="",xlab="discriminant function",ylim=c(0,0.8),main="",mgp=c(2.3,1,0))

polygon(density(pred5$ind.scores[tr.GOL$Treatment=="L"]),col=rgb(0,0,1,alpha=0.5),border="blue")
polygon(density(pred5$ind.scores[tr.GOL$Treatment=="H"]),col=rgb(1,0,0,alpha=0.5),border="red")

"""
