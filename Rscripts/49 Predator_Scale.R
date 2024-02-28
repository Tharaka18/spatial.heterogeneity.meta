########################################################################################################################################################################################

# Article title: Investigating the biodiversity benefits of different kinds of heterogeneity in agricultural landscapes: Comparing crop versus landscape and compositional versus configurational heterogeneity

# Journal: Ecology Letters

# Authors: Tharaka S. Priyadarshana (ORCID: 0000-0003-3962-5465), Emily A. Martin (0000-0001-5785-9105), Clélia Sirami (0000-0003-1741-3082), Ben A. Woodcock (0000-0003-0300-9951), Eben Goodale (0000-0003-3403-2847), Carlos Martínez-Núñez (0000-0001-7814-4985), Myung-Bok Lee (0000-0003-2680-5707), Emilio Pagani-Núñez (0000-0001-8839-4005), Chloé A. Raderschall (0000-0003-2005-1705), Lluís Brotons (0000-0002-4826-4457), Anushka Rege (0000-0002-8383-0258), Annie Ouin (0000-0001-7046-2719), Teja Tscharntke (0000-0002-4482-3178), Eleanor M. Slade (0000-0002-6108-1196)

# contact: Tharaka S. Priyadarshana; tharakas001@e.ntu.edu.sg; tharakas.priyadarshana@gmail.com

dat<- read.csv ("Data/datPredators.csv", header=TRUE, sep=",",na.strings=TRUE)

# load packages
library(metafor)
library(jtools)
library(multcomp)
library(snow)

# calculate the effect size (r to z)
dat <- escalc(measure="ZCOR", ri=Pearson_Correlation, ni=Sample_Size, data=dat)
head(dat)

# tabulate number of studies/effect sizes in biodiversity metrics by taxa groups
table(dat$Scale_Cat)

########################################################################################################################################################################################

# the effect of spatial heterogeneity (i.e. all crop and landscape heterogeneity components) on predator richness
mod1 <- rma.mv(yi, vi, mods = ~ factor(Scale_Cat) - 1, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness"), method="REML",level=95)
summary(mod1)

# back-transform the results of correlations and obtain the prediction interval
summary(mod1, atransf=transf.ztor)

# use cluster-robust inference methods 
rob1 <- robust(mod1, cluster=dat$Study_ID)
summary(rob1, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod1, sigma2=1)
profile(mod1, sigma2=2)

# check for outliers 
(x1 <- cooks.distance(mod1, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x1, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

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

# test difference between the outcomes
summary(glht(rob1, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("none"))
# with Benjamini & Hochberg (1995) correction
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
com1 <- summary(glht(rob1, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("BH"))
summary(com1)
confint(com1)

# graphical display 
t1 <- summary(rob1, atransf=transf.ztor)
plot_summs(t1,
           inner_ci_level=0.90,
           point.shape = FALSE,
           point.size = 12.5,
           colors = "purple")

########################################################################################################################################################################################

# sub-group models
# this is just to count how many studies and effect sizes included in each test
table(dat$Scale_Cat)
A_LocalScale <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness")&(Scale_Cat=="A_LocalScale"), method="REML",level=95)
summary (A_LocalScale)

LandscapeScale_1 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness")&(Scale_Cat=="LandscapeScale_1"), method="REML",level=95)
summary (LandscapeScale_1)

LandscapeScale_2 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness")&(Scale_Cat=="LandscapeScale_2"), method="REML",level=95)
summary (LandscapeScale_2)

########################################################################################################################################################################################

# the effect of spatial heterogeneity (i.e. all crop and landscape heterogeneity components) on predator abundance
mod2 <- rma.mv(yi, vi, mods = ~ factor(Scale_Cat) - 1, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance"), method="REML",level=95)
summary(mod2)

# back-transform the results of correlations and obtain the prediction interval
summary(mod2, atransf=transf.ztor)

# use cluster-robust inference methods 
rob2 <- robust(mod2, cluster=dat$Study_ID)
summary(rob2, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod2, sigma2=1)
profile(mod2, sigma2=2)

# check for outliers 
(x2 <- cooks.distance(mod2, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x2, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

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

# test difference between the outcomes
summary(glht(rob2, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("none"))
# with Benjamini & Hochberg (1995) correction
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
com2 <- summary(glht(rob2, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("BH"))
summary(com2)
confint(com2)

# graphical display 
t2 <- summary(rob2, atransf=transf.ztor)
plot_summs(t2,
           inner_ci_level=0.90,
           point.shape = FALSE,
           point.size = 12.5,
           colors = "purple")

########################################################################################################################################################################################

# sub-group models
# this is just to count how many studies and effect sizes included in each test
table(dat$Scale_Cat)
A_LocalScale <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance")&(Scale_Cat=="A_LocalScale"), method="REML",level=95)
summary (A_LocalScale)

LandscapeScale_1 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance")&(Scale_Cat=="LandscapeScale_1"), method="REML",level=95)
summary (LandscapeScale_1)

LandscapeScale_2 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance")&(Scale_Cat=="LandscapeScale_2"), method="REML",level=95)
summary (LandscapeScale_2)

########################################################################################################################################################################################

# the effect of spatial heterogeneity (i.e. all crop and landscape heterogeneity components) on predator diversity
mod3 <- rma.mv(yi, vi, mods = ~ factor(Scale_Cat) - 1, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity"), method="REML",level=95)
summary(mod3)

# back-transform the results of correlations and obtain the prediction interval
summary(mod3, atransf=transf.ztor)

# use cluster-robust inference methods 
rob3 <- robust(mod3, cluster=dat$Study_ID)
summary(rob3, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod3, sigma2=1)
profile(mod3, sigma2=2)

# check for outliers 
(x3 <- cooks.distance(mod1, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x3, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

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

# test difference between the outcomes
summary(glht(rob3, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("none"))
# with Benjamini & Hochberg (1995) correction
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
com3 <- summary(glht(rob3, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("BH"))
summary(com3)
confint(com3)

# graphical display 
t3 <- summary(rob3, atransf=transf.ztor)
plot_summs(t3,
           inner_ci_level=0.90,
           point.shape = FALSE,
           point.size = 12.5,
           colors = "purple")

########################################################################################################################################################################################

# sub-group models
# this is just to count how many studies and effect sizes included in each test
table(dat$Scale_Cat)
A_LocalScale <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Scale_Cat=="A_LocalScale"), method="REML",level=95)
summary (A_LocalScale)

LandscapeScale_1 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Scale_Cat=="LandscapeScale_1"), method="REML",level=95)
summary (LandscapeScale_1)

LandscapeScale_2 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Scale_Cat=="LandscapeScale_2"), method="REML",level=95)
summary (LandscapeScale_2)

########################################################################################################################################################################################

# the effect of spatial compositional heterogeneity (i.e. crop and landscape compositional heterogeneity components) on predator richness
table(dat$Heterogeneity)
mod4 <- rma.mv(yi, vi, mods = ~ factor(Scale_Cat) - 1, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness")&(Heterogeneity=="Compositional_het"), method="REML",level=95)
summary(mod4)

# back-transform the results of correlations and obtain the prediction interval
summary(mod4, atransf=transf.ztor)

# use cluster-robust inference methods 
rob4 <- robust(mod4, cluster=dat$Study_ID)
summary(rob4, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod4, sigma2=1)
profile(mod4, sigma2=2)

# check for outliers 
(x4 <- cooks.distance(mod4, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x4, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# total amount of heterogeneity 
round(sum(mod4$sigma2), 4)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod4$vi)
X <- model.matrix(mod4)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod4$sigma2) / (sum(mod4$sigma2) + (mod4$k-mod4$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod4$sigma2 / (sum(mod4$sigma2) + (mod4$k-mod4$p)/sum(diag(P)))

# test difference between the outcomes
summary(glht(rob4, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("none"))
# with Benjamini & Hochberg (1995) correction
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
com4 <- summary(glht(rob4, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("BH"))
summary(com4)
confint(com4)

# graphical display 
t4 <- summary(rob4, atransf=transf.ztor)
plot_summs(t4,
           inner_ci_level=0.90,
           point.shape = FALSE,
           point.size = 12.5,
           colors = "purple")

########################################################################################################################################################################################

# sub-group models
# this is just to count how many studies and effect sizes included in each test
table(dat$Scale_Cat)
A_LocalScale <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness")&(Heterogeneity=="Compositional_het")&(Scale_Cat=="A_LocalScale"), method="REML",level=95)
summary (A_LocalScale)

LandscapeScale_1 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness")&(Heterogeneity=="Compositional_het")&(Scale_Cat=="LandscapeScale_1"), method="REML",level=95)
summary (LandscapeScale_1)

LandscapeScale_2 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness")&(Heterogeneity=="Compositional_het")&(Scale_Cat=="LandscapeScale_2"), method="REML",level=95)
summary (LandscapeScale_2)

########################################################################################################################################################################################

# the effect of spatial compositional heterogeneity (i.e. crop and landscape compositional heterogeneity components) on predator diversity
table(dat$Scale_Cat)
mod5 <- rma.mv(yi, vi, mods = ~ factor(Scale_Cat) - 1, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Heterogeneity=="Compositional_het"), method="REML",level=95)
summary(mod5)

# back-transform the results of correlations and obtain the prediction interval
summary(mod5, atransf=transf.ztor)

# use cluster-robust inference methods 
rob5 <- robust(mod5, cluster=dat$Study_ID)
summary(rob5, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod5, sigma2=1)
profile(mod5, sigma2=2)

# check for outliers 
(x5 <- cooks.distance(mod5, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x5, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# total amount of heterogeneity 
round(sum(mod5$sigma2), 4)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod5$vi)
X <- model.matrix(mod5)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod5$sigma2) / (sum(mod5$sigma2) + (mod5$k-mod5$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod5$sigma2 / (sum(mod5$sigma2) + (mod5$k-mod5$p)/sum(diag(P)))

# test difference between the outcomes
summary(glht(rob5, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("none"))
# with Benjamini & Hochberg (1995) correction
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
com5 <- summary(glht(rob5, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("BH"))
summary(com5)
confint(com5)

# graphical display 
t5 <- summary(rob5, atransf=transf.ztor)
plot_summs(t5,
           inner_ci_level=0.90,
           point.shape = FALSE,
           point.size = 12.5,
           colors = "purple")

########################################################################################################################################################################################

# sub-group models
# this is just to count how many studies and effect sizes included in each test
table(dat$Scale_Cat)
A_LocalScale <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Heterogeneity=="Compositional_het")&(Scale_Cat=="A_LocalScale"), method="REML",level=95)
summary (A_LocalScale)

LandscapeScale_1 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Heterogeneity=="Compositional_het")&(Scale_Cat=="LandscapeScale_1"), method="REML",level=95)
summary (LandscapeScale_1)

LandscapeScale_2 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Heterogeneity=="Compositional_het")&(Scale_Cat=="LandscapeScale_2"), method="REML",level=95)
summary (LandscapeScale_2)

########################################################################################################################################################################################

# the effect of spatial compositional heterogeneity (i.e. crop and landscape compositional heterogeneity components) on predator abundance
table(dat$Scale_Cat)
mod6 <- rma.mv(yi, vi, mods = ~ factor(Scale_Cat) - 1, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance")&(Heterogeneity=="Compositional_het"), method="REML",level=95)
summary(mod6)

# back-transform the results of correlations and obtain the prediction interval
summary(mod6, atransf=transf.ztor)

# use cluster-robust inference methods 
rob6 <- robust(mod6, cluster=dat$Study_ID)
summary(rob6, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod6, sigma2=1)
profile(mod6, sigma2=2)

# check for outliers 
(x6 <- cooks.distance(mod6, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x6, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# total amount of heterogeneity 
round(sum(mod6$sigma2), 4)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod6$vi)
X <- model.matrix(mod6)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod6$sigma2) / (sum(mod6$sigma2) + (mod6$k-mod6$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod6$sigma2 / (sum(mod6$sigma2) + (mod6$k-mod6$p)/sum(diag(P)))

# test difference between the outcomes
summary(glht(rob6, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("none"))
# with Benjamini & Hochberg (1995) correction
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
com6 <- summary(glht(rob6, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("BH"))
summary(com6)
confint(com6)

# graphical display 
t6 <- summary(rob6, atransf=transf.ztor)
plot_summs(t6,
           inner_ci_level=0.90,
           point.shape = FALSE,
           point.size = 12.5,
           colors = "purple")

########################################################################################################################################################################################

# sub-group models
# this is just to count how many studies and effect sizes included in each test
table(dat$Scale_Cat)
A_LocalScale <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance")&(Heterogeneity=="Compositional_het")&(Scale_Cat=="A_LocalScale"), method="REML",level=95)
summary (A_LocalScale)

LandscapeScale_1 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance")&(Heterogeneity=="Compositional_het")&(Scale_Cat=="LandscapeScale_1"), method="REML",level=95)
summary (LandscapeScale_1)

LandscapeScale_2 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance")&(Heterogeneity=="Compositional_het")&(Scale_Cat=="LandscapeScale_2"), method="REML",level=95)
summary (LandscapeScale_2)

########################################################################################################################################################################################

# the effect of spatial configurational heterogeneity (i.e. crop and landscape configurational heterogeneity components) on predator richness
table(dat$Heterogeneity)
mod7 <- rma.mv(yi, vi, mods = ~ factor(Scale_Cat) - 1, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness")&(Heterogeneity=="Configurational_het"), method="REML",level=95)
summary(mod7)

# back-transform the results of correlations and obtain the prediction interval
summary(mod7, atransf=transf.ztor)

# use cluster-robust inference methods 
rob7 <- robust(mod7, cluster=dat$Study_ID)
summary(rob7, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod7, sigma2=1)
profile(mod7, sigma2=2)

# check for outliers 
(x7 <- cooks.distance(mod7, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x7, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# total amount of heterogeneity 
round(sum(mod7$sigma2), 4)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod7$vi)
X <- model.matrix(mod7)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod7$sigma2) / (sum(mod7$sigma2) + (mod7$k-mod7$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod7$sigma2 / (sum(mod7$sigma2) + (mod7$k-mod7$p)/sum(diag(P)))

# test difference between the outcomes
summary(glht(rob7, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("none"))
# with Benjamini & Hochberg (1995) correction
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
com7 <- summary(glht(rob7, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("BH"))
summary(com7)
confint(com7)

# graphical display 
t7 <- summary(rob7, atransf=transf.ztor)
plot_summs(t7,
           inner_ci_level=0.90,
           point.shape = FALSE,
           point.size = 12.5,
           colors = "purple")

########################################################################################################################################################################################

# sub-group models
# this is just to count how many studies and effect sizes included in each test
table(dat$Scale_Cat)
A_LocalScale <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness")&(Heterogeneity=="Configurational_het")&(Scale_Cat=="A_LocalScale"), method="REML",level=95)
summary (A_LocalScale)

LandscapeScale_1 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness")&(Heterogeneity=="Configurational_het")&(Scale_Cat=="LandscapeScale_1"), method="REML",level=95)
summary (LandscapeScale_1)

LandscapeScale_2 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Richness")&(Heterogeneity=="Configurational_het")&(Scale_Cat=="LandscapeScale_2"), method="REML",level=95)
summary (LandscapeScale_2)

########################################################################################################################################################################################

# the effect of spatial configurational heterogeneity (i.e. crop and landscape configurational heterogeneity components) on predator diversity
table(dat$Scale_Cat)
mod8 <- rma.mv(yi, vi, mods = ~ factor(Scale_Cat) - 1, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Heterogeneity=="Configurational_het"), method="REML",level=95)
summary(mod8)

# back-transform the results of correlations and obtain the prediction interval
summary(mod8, atransf=transf.ztor)

# use cluster-robust inference methods 
rob8 <- robust(mod8, cluster=dat$Study_ID)
summary(rob8, atransf=transf.ztor)

# Profile Likelihood Plots
par(mfrow=c(2,1))
profile(mod8, sigma2=1)
profile(mod8, sigma2=2)

# check for outliers 
(x8 <- cooks.distance(mod8, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x8, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# total amount of heterogeneity 
round(sum(mod8$sigma2), 4)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod8$vi)
X <- model.matrix(mod8)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod8$sigma2) / (sum(mod8$sigma2) + (mod8$k-mod8$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod8$sigma2 / (sum(mod8$sigma2) + (mod8$k-mod8$p)/sum(diag(P)))

# test difference between the outcomes
summary(glht(rob8, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("none"))
# with Benjamini & Hochberg (1995) correction
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
com8 <- summary(glht(rob8, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("BH"))
summary(com8)
confint(com8)

# graphical display 
t8 <- summary(rob8, atransf=transf.ztor)
plot_summs(t8,
           inner_ci_level=0.90,
           point.shape = FALSE,
           point.size = 12.5,
           colors = "purple")

########################################################################################################################################################################################

# sub-group models
# this is just to count how many studies and effect sizes included in each test
table(dat$Scale_Cat)
A_LocalScale <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Heterogeneity=="Configurational_het")&(Scale_Cat=="A_LocalScale"), method="REML",level=95)
summary (A_LocalScale)

LandscapeScale_1 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Heterogeneity=="Configurational_het")&(Scale_Cat=="LandscapeScale_1"), method="REML",level=95)
summary (LandscapeScale_1)

LandscapeScale_2 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Heterogeneity=="Configurational_het")&(Scale_Cat=="LandscapeScale_2"), method="REML",level=95)
summary (LandscapeScale_2)

########################################################################################################################################################################################

# the effect of spatial configurational heterogeneity (i.e. crop and landscape configurational heterogeneity components) on predator abundance
table(dat$Scale_Cat)
mod9 <- rma.mv(yi, vi, mods = ~ factor(Scale_Cat) - 1, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance")&(Heterogeneity=="Configurational_het"), method="REML",level=95)
summary(mod9)

# back-transform the results of correlations and obtain the prediction interval
summary(mod9, atransf=transf.ztor)

# use cluster-robust inference methods 
rob9 <- robust(mod9, cluster=dat$Study_ID)
summary(rob9, atransf=transf.ztor)

# profile likelihood plots
par(mfrow=c(2,1))
profile(mod9, sigma2=1)
profile(mod9, sigma2=2)

# check for outliers 
(x9 <- cooks.distance(mod9, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x9, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# total amount of heterogeneity 
round(sum(mod9$sigma2), 4)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod9$vi)
X <- model.matrix(mod9)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod9$sigma2) / (sum(mod9$sigma2) + (mod9$k-mod9$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod9$sigma2 / (sum(mod9$sigma2) + (mod9$k-mod9$p)/sum(diag(P)))

# test difference between the outcomes
summary(glht(rob9, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("none"))
# with Benjamini & Hochberg (1995) correction
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
com9 <- summary(glht(rob9, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("BH"))
summary(com9)
confint(com9)

# graphical display 
t9 <- summary(rob9, atransf=transf.ztor)
plot_summs(t9,
           inner_ci_level=0.90,
           point.shape = FALSE,
           point.size = 12.5,
           colors = "purple")

########################################################################################################################################################################################

# sub-group models
# this is just to count how many studies and effect sizes included in each test
table(dat$Scale_Cat)
A_LocalScale <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance")&(Heterogeneity=="Configurational_het")&(Scale_Cat=="A_LocalScale"), method="REML",level=95)
summary (A_LocalScale)

LandscapeScale_1 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance")&(Heterogeneity=="Configurational_het")&(Scale_Cat=="LandscapeScale_1"), method="REML",level=95)
summary (LandscapeScale_1)

LandscapeScale_2 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="Abundance")&(Heterogeneity=="Configurational_het")&(Scale_Cat=="LandscapeScale_2"), method="REML",level=95)
summary (LandscapeScale_2)

########################################################################################################################################################################################
