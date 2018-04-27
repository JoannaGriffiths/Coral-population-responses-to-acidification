setwd("~/Desktop/Griffiths2018_data")

library(lme4)

Data=read.table("lipid_protein_data_manuscript.csv", header = TRUE) #import and attach data
Data
attach(Data)
str(Data) #quicky view data
head(Data)

#making the Jar and Coral numbers factors and not integers
Data$Jar <- as.factor(Data$Jar)
Data$Coral <- as.factor(Data$Coral)
str(Data)
head(Data)

############################### total lipids

tot_lipid1=(tot_lipid/Dry_weight) #normalizing by dry weight
tot_lipid2=sqrt(tot_lipid1) #normalizing by square root

tot_lipid.model1 = lme(tot_lipid2 ~ pH*Location,  random = ~ 1|Jar/Coral, data = Data) 

aov1 = aov(tot_lipid.model1) 
summary(aov1) 

posthoc1 = TukeyHSD(aov1)
posthoc1
plot(posthoc1)

############################### wax ester

wax_ester1=(waxester/Dry_weight)
wax_ester2=sqrt(wax_ester1)

wax_ester.model1 = lme(wax_ester2 ~ pH*Location,  random = ~ 1|Jar/Coral, data = Data)

aov1 = aov(wax_ester.model1) 
summary(aov1) 

posthoc1 = TukeyHSD(aov1)
posthoc1
plot(posthoc1)

############################### triacylglycerol

triacylglycerol1=(triacylglycerol/Dry_weight)
triacylglycerol2=sqrt(triacylglycerol1)

triacylglycerol.model1 = lme(triacylglycerol2 ~ pH*Location,  random = ~ 1|Jar/Coral, data = Data)

aov1 = aov(triacylglycerol.model1) 
summary(aov1)

posthoc1 = TukeyHSD(aov1)
posthoc1
plot(posthoc1)

############################### sterol

sterol1=(sterol/Dry_weight)
sterol2=sqrt(sterol1)

sterol.model1 = lme(sterol2 ~ pH*Location,  random = ~ 1|Jar/Coral, data = Data)

aov1 = aov(sterol.model1) 
summary(aov1) 

posthoc1 = TukeyHSD(aov1)
posthoc1
plot(posthoc1)

############################### fattyacid

fattyacid1=(fattyacid/Dry_weight)
fattyacid2=sqrt(fattyacid1)

fattyacid.model1 = lme(fattyacid2 ~ pH*Location,  random = ~ 1|Jar/Coral, data = Data)

aov1 = aov(fattyacid.model1) 
summary(aov1)

posthoc1 = TukeyHSD(aov1)
posthoc1
plot(posthoc1)

############################### phospholipid

phospholipid1=(phospholipid/Dry_weight)
phospholipid2=sqrt(phospholipid1)

phospholipid.model1 = lme(phospholipid2 ~ pH*Location,  random = ~ 1|Jar/Coral, data = Data)

aov1 = aov(phospholipid.model1) 
summary(aov1)

posthoc1 = TukeyHSD(aov1)
posthoc1
plot(posthoc1)

############################### protein

protein1=(protein/Dry_weight)
protein2=sqrt(protein1)

protein.model1 = lme(protein2 ~ pH*Location,  random = ~ 1|Jar/Coral, data = Data)

aov1 = aov(protein.model1) 
summary(aov1) 

posthoc1 = TukeyHSD(aov1)
posthoc1
plot(posthoc1)
