########################################################################################################################################################################################

# Article title: Investigating the biodiversity benefits of different kinds of heterogeneity in agricultural landscapes: Comparing crop versus landscape and compositional versus configurational heterogeneity

# Journal: Ecology Letters

# Authors: Tharaka S. Priyadarshana (ORCID: 0000-0003-3962-5465), Emily A. Martin (0000-0001-5785-9105), Clélia Sirami (0000-0003-1741-3082), Ben A. Woodcock (0000-0003-0300-9951), Eben Goodale (0000-0003-3403-2847), Carlos Martínez-Núñez (0000-0001-7814-4985), Myung-Bok Lee (0000-0003-2680-5707), Emilio Pagani-Núñez (0000-0001-8839-4005), Chloé A. Raderschall (0000-0003-2005-1705), Lluís Brotons (0000-0002-4826-4457), Anushka Rege (0000-0002-8383-0258), Annie Ouin (0000-0001-7046-2719), Teja Tscharntke (0000-0002-4482-3178), Eleanor M. Slade (0000-0002-6108-1196)

# contact: Tharaka S. Priyadarshana; tharakas001@e.ntu.edu.sg; tharakas.priyadarshana@gmail.com

# testing the effects of mean other land-use area (i.e. total area - crop + semi-natural area) on the estimated average relationships between spatial heterogeneity and biodiversity 

dat<- read.csv ("Data/datOtherArea.csv", header=TRUE, sep=",",na.strings=TRUE)

# load packages
library(metafor)
library(jtools)
library(multcomp)
library(snow)

# calculate the effect size (r to z)
dat <- escalc(measure="ZCOR", ri=Pearson_Correlation, ni=Sample_Size, data=dat)
head(dat)

# tabulate number of studies/effect sizes in biodiversity metrics by taxa groups
table(dat$Verts_Inverts_Plants_Without_Pests, dat$Biodiversity_Category)

# mean and SE/SD of the crop cover across all the studies 
dat$Other_Cover_Mean <- as.numeric(dat$Other_Cover_Mean) 
mean(dat$Other_Cover_Mean)
sd(dat$Other_Cover_Mean)
standard_error <- function(x) sd(x) / sqrt(length(x))
x <- c(dat$Other_Cover_Mean) 
standard_error(x) 

# testing the effects of mean other land-use area (i.e. total area - crop + semi-natural area) on the relationship between spatial heterogeneity and animal (i.e. invertebrate and vertebrate without pests) richness 
mod1 <- rma.mv(yi, vi, mods = ~ Other_Cover_Mean, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness")&(Animals_Plants=="Animals"), method="REML", test="t", level=95)
summary(mod1)

# back-transform the results of correlations and obtain the prediction intervals
summary(mod1, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod1, sigma2=1)
profile(mod1, sigma2=2)

# check for outliers 
(x1 <- cooks.distance(mod1, parallel = "snow", ncpus=40)) #this may take some time, depend on the power of the computer. instadded "snow" one can use parallel="multicore". Here ncpus means the integer to specify the number of processes to use in the parallel processing.
plot(x1, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# intraclass correlation (ICC) of the true effects
# compute the ICC (i.e. estimated correlation of the effects belonging to the same study, or how strong the dependency of true effect coming from the same study, or how strong the correlation of within study  heterogeneity)
round(mod1$sigma2[1] / sum(mod1$sigma2), 4)

# total amount of heterogeneity 
round(sum(mod1$sigma2), 4)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod1$vi)
X <- model.matrix(mod1)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod1$sigma2) / (sum(mod1$sigma2) + (mod1$k-mod1$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod1$sigma2 / (sum(mod1$sigma2) + (mod1$k-mod1$p)/sum(diag(P)))

########################################################################################################################################################################################

# testing the effects of mean other land-use area (i.e. total area - crop + semi-natural area) on the relationship between spatial heterogeneity and animal (i.e. invertebrate and vertebrate without pests) diversity 
mod2 <- rma.mv(yi, vi, mods = ~ Other_Cover_Mean, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Animals_Plants=="Animals"), method="REML", test="t", level=95)
summary(mod2)

# back-transform the results of correlations and obtain the prediction interval
summary(mod2, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod2, sigma2=1)
profile(mod2, sigma2=2)

# check for outliers 
(x2 <- cooks.distance(mod2, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x2, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# intraclass correlation (ICC) of the true effects
# compute the ICC (i.e. estimated correlation of the effects belonging to the same study, or how strong the dependency of true effect coming from the same study, or how strong the correlation of within study  heterogeneity)
round(mod2$sigma2[1] / sum(mod2$sigma2), 4)

# total amount of heterogeneity 
round(sum(mod2$sigma2), 4)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod2$vi)
X <- model.matrix(mod2)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod2$sigma2) / (sum(mod2$sigma2) + (mod2$k-mod2$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod2$sigma2 / (sum(mod2$sigma2) + (mod2$k-mod2$p)/sum(diag(P)))

########################################################################################################################################################################################

# testing the effects of mean other land-use area (i.e. total area - crop + semi-natural area) on the relationship between spatial heterogeneity and animal (i.e. invertebrate and vertebrate without pests) abundance 
mod3 <-  rma.mv(yi, vi, mods = ~ Other_Cover_Mean, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance")&(Animals_Plants=="Animals"), method="REML", test="t", level=95)
summary(mod3)

# back-transform the results of correlations and obtain the prediction interval
summary(mod3, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod3, sigma2=1)
profile(mod3, sigma2=2)

# check for outliers 
(x3 <- cooks.distance(mod3, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x3, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# intraclass correlation (ICC) of the True Effects
# compute the ICC (i.e. estimated correlation of the effects belonging to the same study, or how strong the dependency of true effect coming from the same study, or how strong the correlation of within study  heterogeneity)
round(mod3$sigma2[1] / sum(mod3$sigma2), 4)

# total amount of heterogeneity 
round(sum(mod3$sigma2), 4)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod3$vi)
X <- model.matrix(mod3)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod3$sigma2) / (sum(mod3$sigma2) + (mod3$k-mod3$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod3$sigma2 / (sum(mod3$sigma2) + (mod3$k-mod3$p)/sum(diag(P)))

########################################################################################################################################################################################

# testing the effects of mean other land-use area (i.e. total area - crop + semi-natural area) on the relationship between spatial heterogeneity and plant biodiversity

dat<- read.csv ("Data/datOtherArea.csv", header=TRUE, sep=",",na.strings=TRUE)

# load packages
library(metafor)
library(jtools)
library(multcomp)
library(snow)

# calculate the effect size (r to z)
dat <- escalc(measure="ZCOR", ri=Pearson_Correlation, ni=Sample_Size, data=dat)
head(dat)

# tabulate number of studies/effect sizes in biodiversity metrics by taxa groups
table(dat$Verts_Inverts_Plants_Without_Pests, dat$Biodiversity_Category)

dat$Seminatural_Cover_Mean <- as.numeric(dat$Seminatural_Cover_Mean) 

# testing the effects of mean other land-use area (i.e. total area - crop + semi-natural area) on the relationship between spatial heterogeneity and plant richness 
mod1 <- rma.mv(yi, vi, mods = ~ Other_Cover_Mean, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness")&(Animals_Plants=="Plants"), method="REML", test="t", level=95)
summary(mod1)

# back-transform the results of correlations and obtain the prediction interval
summary(mod1, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod1, sigma2=1)
profile(mod1, sigma2=2)

# check for outliers 
(x1 <- cooks.distance(mod1, parallel = "snow", ncpus=40)) #this may take some time, depend on the power of the computer. instadded "snow" one can use parallel="multicore". Here ncpus means the integer to specify the number of processes to use in the parallel processing.
plot(x1, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# intraclass correlation (ICC) of the true effects
# compute the ICC (i.e. estimated correlation of the effects belonging to the same study, or how strong the dependency of true effect coming from the same study, or how strong the correlation of within study  heterogeneity)
round(mod1$sigma2[1] / sum(mod1$sigma2), 4)

# total amount of heterogeneity 
round(sum(mod1$sigma2), 4)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod1$vi)
X <- model.matrix(mod1)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod1$sigma2) / (sum(mod1$sigma2) + (mod1$k-mod1$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod1$sigma2 / (sum(mod1$sigma2) + (mod1$k-mod1$p)/sum(diag(P)))

########################################################################################################################################################################################

# testing the effects of mean other land-use area (i.e. total area - crop + semi-natural area) on the relationship between spatial heterogeneity and plant diversity 
mod2 <- rma.mv(yi, vi, mods = ~ Other_Cover_Mean, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Animals_Plants=="Plants"), method="REML", test="t", level=95)
summary(mod2)

# back-transform the results of correlations and obtain the prediction interval
summary(mod2, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod2, sigma2=1)
profile(mod2, sigma2=2)

# check for outliers 
(x2 <- cooks.distance(mod2, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x2, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# intraclass correlation (ICC) of the true effects
# compute the ICC (i.e. estimated correlation of the effects belonging to the same study, or how strong the dependency of true effect coming from the same study, or how strong the correlation of within study  heterogeneity)
round(mod2$sigma2[1] / sum(mod2$sigma2), 4)

# total amount of heterogeneity 
round(sum(mod2$sigma2), 4)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod2$vi)
X <- model.matrix(mod2)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod2$sigma2) / (sum(mod2$sigma2) + (mod2$k-mod2$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod2$sigma2 / (sum(mod2$sigma2) + (mod2$k-mod2$p)/sum(diag(P)))

########################################################################################################################################################################################

# testing the effects of mean other land-use area (i.e. total area - crop + semi-natural area) on the relationship between spatial heterogeneity and plant abundance 
mod3 <-  rma.mv(yi, vi, mods = ~ Other_Cover_Mean, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance")&(Animals_Plants=="Plants"), method="REML", test="t", level=95)
summary(mod3)

# back-transform the results of correlations and obtain the prediction interval
summary(mod3, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod3, sigma2=1)
profile(mod3, sigma2=2)

# check for outliers 
(x3 <- cooks.distance(mod3, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x3, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# intraclass correlation (ICC) of the true effects
# compute the ICC (i.e. estimated correlation of the effects belonging to the same study, or how strong the dependency of true effect coming from the same study, or how strong the correlation of within study  heterogeneity)
round(mod3$sigma2[1] / sum(mod3$sigma2), 4)

# total amount of heterogeneity 
round(sum(mod3$sigma2), 4)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod3$vi)
X <- model.matrix(mod3)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod3$sigma2) / (sum(mod3$sigma2) + (mod3$k-mod3$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod3$sigma2 / (sum(mod3$sigma2) + (mod3$k-mod3$p)/sum(diag(P)))
