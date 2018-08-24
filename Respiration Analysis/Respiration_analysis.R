setwd("~/Desktop/Griffiths2018_data")


library(lme4)
library(car)
library(emmeans)

Resp14=read.table("Resp1-4_umol_v2.csv", header = TRUE)
Resp14

Resp14$Jar <- as.factor(Resp14$Jar)
Resp14$Coral <- as.factor(Resp14$Coral)
Resp14$Week <- as.factor(Resp14$Week)
attach(Resp14)
str(Resp14)
head(Resp14)

abResp14=abs(Resp14$umolghr) #normalizing data
abResp14
Resp14$lgResp14=log(abResp14) 
lgResp14 = log(abResp14)

Resp4HU_L <- c(-0.50080676, -1.24811931, -1.20419333, -1.4574839, -0.9196009, -0.60189205, -0.08969823, -1.27003099, -0.74468992)
Resp4HU_H <- c(0.64724495, 0.85766221, -0.8422587, -1.25309985, -0.75145199, -0.8163816, -1.31263915, -0.98143827, -1.06324171)
Resp4LU_L <- c(-1.6467178, -1.71574565, -1.75187842, -4.63752707, -0.28805273, -1.85820188, -1.22042155, -1.09866861, -2.76825158)
Resp4LU_H <- c(-2.25851883, -3.19195594, -0.69002833, -0.60867025, -0.49310615, -0.84649393, -1.67568734, -1.09543746, -1.24130478)
shapiro.test(Resp4HU_L)
shapiro.test(Resp4HU_H)
shapiro.test(Resp4LU_H)
shapiro.test(Resp4LU_L)

resp.model14 = lmer(lgResp14 ~ pH*Location*Week + (1|Jar/Coral), data = Resp14)
summary(resp.model14)

Anova(resp.model14,test.statistic = "F") 

lsmeans(resp.model14, list(pairwise ~ pH*Location*Week), adjust = "tukey") 
