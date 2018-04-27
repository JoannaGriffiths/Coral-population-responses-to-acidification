setwd("~/Desktop/Griffiths2018_data")

library(lme4)

Resp1=read.table("Resp_day9_manuscript.csv", header = TRUE) #read in data
str(Resp1) #quickly view data
head(Resp1)

#making the Jar and Coral numbers facotrs and not integers
Resp1$Jar <- as.factor(Resp1$Jar)
Resp1$Coral <- as.factor(Resp1$Coral)
str(Resp1)
head(Resp1)

#normalizaing the data
attach(Resp1)
abResp1=abs(umolghr)
abResp1
lgResp1=log(abResp1)
lgResp1

hist(lgResp1, main="Histogram of Respiration Week 1")
shapiro.test(lgResp1)

resp.model1 = lme(lgResp1 ~ pH*Location,  random = ~ 1|Jar/Coral, data = Resp1)

aov1 = aov(resp.model1) 
summary(aov1)


############################################ same as above but with week 4

Resp4=read.table("Resp_day29_manuscript.csv", header = TRUE)
str(Resp4)
head(Resp4)

#making the Jar and Coral numbers facotrs and not integers
Resp4$Jar <- as.factor(Resp4$Jar)
Resp4$Coral <- as.factor(Resp4$Coral)
str(Resp4)
head(Resp4)

#normalizaing the data
attach(Resp4)
abResp4=abs(umolghr)
abResp4
lgResp4=log(abResp4)
lgResp4


resp.model4 = lme(lgResp4 ~ pH*Location,  random = ~ 1|Jar/Coral, data = Resp4)

aov4 = aov(resp.model4) 
aov4
summary(aov4)

