########################################################################################################################################################################################

# Article title: Crop and landscape heterogeneity increase biodiversity in agricultural landscapes: A global review and meta-analysis

# Journal: Ecology Letters

# Authors: Tharaka S. Priyadarshana (ORCID: 0000-0003-3962-5465), Emily A. Martin (0000-0001-5785-9105), Clélia Sirami (0000-0003-1741-3082), Ben A. Woodcock (0000-0003-0300-9951), Eben Goodale (0000-0003-3403-2847), Carlos Martínez-Núñez (0000-0001-7814-4985), Myung-Bok Lee (0000-0003-2680-5707), Emilio Pagani-Núñez (0000-0001-8839-4005), Chloé A. Raderschall (0000-0003-2005-1705), Lluís Brotons (0000-0002-4826-4457), Anushka Rege (0000-0002-8383-0258), Annie Ouin (0000-0001-7046-2719), Teja Tscharntke (0000-0002-4482-3178), Eleanor M. Slade (0000-0002-6108-1196)

# contact: Tharaka S. Priyadarshana; tharakas001@e.ntu.edu.sg; tharakas.priyadarshana@gmail.com

dat<- read.csv ("Data/dat.csv", header=TRUE, sep=",",na.strings=TRUE)

# load packages
library(metafor)
library(jtools)
library(multcomp)
library(snow)

# calculate the effect size (r to z)
dat <- escalc(measure="ZCOR", ri=Pearson_Correlation, ni=Sample_Size, data=dat)
head(dat)

# tabulate number of studies/effect sizes in biodiversity metrics by taxa groups
table(dat$Heterogeneity)
table(dat$Biodiversity_Category)
table(dat$Heterogeneity, dat$Biodiversity_Category)
table(dat$Verts_Inverts_Plants_Without_Pests, dat$Biodiversity_Category)
table(dat$Study_ID, dat$Biodiversity_Category, dat$Heterogeneity_Type)

# the effect of spatial heterogeneity (i.e. all crop and landscape heterogeneity components) on plants diversity
mod1 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Verts_Inverts_Plants_Without_Pests=="Plants"), method="REML",level=95)
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
(x1 <- cooks.distance(mod1, parallel = "snow", ncpus=40)) #this may take some time, depend on the power of the computer. instadded "snow" one can use parallel="multicore". Here ncpus means the integer to specify the number of processes to use in the parallel processing.
plot(x1, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# intraclass correlation (ICC) of the true effects
# compute the ICC (i.e. estimated correlation of the effects belonging to the same study, or how strong the dependency of true effect coming from the same study, or how strong the correlation of within study  heterogeneity)
round(mod1$sigma2[1] / sum(mod1$sigma2), 4)

# total amount of heterogeneity 
round(sum(mod1$sigma2), 0)

# forest plot 

#forest(mod1, addpred=TRUE, header=TRUE, atransf=transf.ztor, annotate=TRUE, order="obs", xlim=c(-1.5 , 1.5), addfit=TRUE, slab=dat$Author_With_Scale)

#check for publication bias 

# standard funnel plot with the SEs on the y-axis
#funnel(mod1)

# contour-enhanced fuel plot
#funnel(mod1, xlim=c(-1.5,1.5), ylim=c(0,0.5), refline=0, level=c(90,95,99),
       #shade=c("gray90","gray65","gray55"), back="white", hlines=NULL, las=1, 
       #legend=TRUE, digits=list(1L,1))

# put 1/sqrt(ni) on the y-axis
#funnel(mod1, yaxis="sqrtninv", xlim=c(-1.5,1.5), ylim=c(0,0.4), las=2)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod1$vi)
X <- model.matrix(mod1)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod1$sigma2) / (sum(mod1$sigma2) + (mod1$k-mod1$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod1$sigma2 / (sum(mod1$sigma2) + (mod1$k-mod1$p)/sum(diag(P)))

# graphical display 
t1 <- summary(rob1, atransf=transf.ztor)
plot_summs(t1,
           inner_ci_level=0.90,
           colors = "red")

########################################################################################################################################################################################

# the effects of spatial compositional heterogeneity (i.e. crop and landscape compositional heterogeneity components) and spatial configurational heterogeneity (i.e. crop and landscape configurational heterogeneity components) on plant diversity
mod2 <- rma.mv(yi, vi, mods = ~ factor(Heterogeneity) - 1, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Verts_Inverts_Plants_Without_Pests=="Plants"), method="REML",level=95)
summary(mod2)

# back-transform the results of correlations and obtain the prediction interval
summary(mod2, atransf=transf.ztor)

# use cluster-robust inference methods
rob2 <- robust(mod2, cluster=dat$Study_ID)
summary(rob2, atransf=transf.ztor)

# Profile Likelihood Plots
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

# forest plot 

#forest(mod1, addpred=TRUE, header=TRUE, atransf=transf.ztor, annotate=TRUE, order="obs", xlim=c(-1.5 , 1.5), addfit=TRUE, slab=dat$Author_With_Scale)

#check for publication bias 

# standard funnel plot with the SEs on the y-axis
#funnel(mod2)

# contour-enhanced fuel plot
#funnel(mod2, xlim=c(-1.5,1.5), ylim=c(0,0.5), refline=0, level=c(90,95,99),
       #shade=c("gray90","gray65","gray55"), back="white", hlines=NULL, las=1, 
       #legend=TRUE, digits=list(1L,1))

# put 1/sqrt(ni) on the y-axis
#funnel(mod2, yaxis="sqrtninv", xlim=c(-1,1), ylim=c(0,0.5), las=1)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod2$vi)
X <- model.matrix(mod2)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod2$sigma2) / (sum(mod2$sigma2) + (mod2$k-mod2$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod2$sigma2 / (sum(mod2$sigma2) + (mod2$k-mod2$p)/sum(diag(P)))

# test difference between the two outcomes (crop/landscape compositional and crop/landscape configurational heterogeneity)
anova(mod2, L=c(1, -1),digits = 4)

# graphical display 
t2 <- summary(rob2, atransf=transf.ztor)
plot_summs(t2,
           inner_ci_level=0.90,
           colors = "purple")

########################################################################################################################################################################################
# sub-group models
# this is just to count how many studies and effect sizes included in each test
table(dat$Heterogeneity)
Compo <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Verts_Inverts_Plants_Without_Pests=="Plants")&(Heterogeneity=="Compositional_het"), method="REML",level=95)
summary (Compo)

table(dat$Heterogeneity)
Config <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Verts_Inverts_Plants_Without_Pests=="Plants")&(Heterogeneity=="Configurational_het"), method="REML",level=95)
summary (Config)

########################################################################################################################################################################################

# the effects of crop compositional heterogeneity, crop configurational heterogeneity, landscape compositional and landscape configurational heterogeneity on plant diversity
mod3 <- rma.mv(yi, vi, mods = ~ factor(Heterogeneity_Type) - 1, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Verts_Inverts_Plants_Without_Pests=="Plants"), method="REML",level=95)
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
(x3 <- cooks.distance(mod3, parallel="multicore")) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x3, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) #yes all below 1 

# intraclass correlation (ICC) of the true effects
# compute the ICC (i.e. estimated correlation of the effects belonging to the same study, or how strong the dependency of true effect coming from the same study, or how strong the correlation of within study  heterogeneity)
round(mod3$sigma2[1] / sum(mod3$sigma2), 4)

# total amount of heterogeneity 
round(sum(mod3$sigma2), 4)

# forest plot 

#forest(mod1, addpred=TRUE, header=TRUE, atransf=transf.ztor, annotate=TRUE, order="obs", xlim=c(-1.5 , 1.5), addfit=TRUE, slab=dat$Author_With_Scale)

#check for publication bias 

# standard funnel plot with the SEs on the y-axis
#funnel(mod3)

# contour-enhanced fuel plot
#funnel(mod3, xlim=c(-1.5,1.5), ylim=c(0,0.5), refline=0, level=c(90,95,99),
       #shade=c("gray90","gray65","gray55"), back="white", hlines=NULL, las=1, 
       #legend=TRUE, digits=list(1L,1))

# put 1/sqrt(ni) on the y-axis
#funnel(mod3, yaxis="sqrtninv", xlim=c(-1,1), ylim=c(0,0.5), las=1)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod3$vi)
X <- model.matrix(mod3)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod3$sigma2) / (sum(mod3$sigma2) + (mod3$k-mod3$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod3$sigma2 / (sum(mod3$sigma2) + (mod3$k-mod3$p)/sum(diag(P)))

# test difference between the four outcomes (crop/landscape compositional and crop/landscape configurational heterogeneity)
summary(glht(mod3, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("none"))
# with Benjamini & Hochberg (1995) correction
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
com1 <- summary(glht(mod3, linfct=cbind(contrMat(rep(1,3), type="Tukey"))), test=adjusted("BH"))
summary(com1)
confint(com1)

# graphical display 
t3 <- summary(rob3, atransf=transf.ztor)
plot_summs(t3,
           inner_ci_level=0.90,
           colors = "yellow")

########################################################################################################################################################################################
# sub-group models
# this is just to count how many studies and effect sizes included in each test
table(dat$Heterogeneity_Type)
CropCompo <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Verts_Inverts_Plants_Without_Pests=="Plants")&(Heterogeneity_Type=="Crop_compositional_het"), method="REML",level=95)
summary (CropCompo)

CropConfig <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Verts_Inverts_Plants_Without_Pests=="Plants")&(Heterogeneity_Type=="Crop_Configurational_het"), method="REML",level=95)
summary (CropConfig)

LandscapeCompo <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Verts_Inverts_Plants_Without_Pests=="Plants")&(Heterogeneity_Type=="Landscape_compositional_het"), method="REML",level=95)
summary (LandscapeCompo)

LandscapeConfig <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Verts_Inverts_Plants_Without_Pests=="Plants")&(Heterogeneity_Type=="Landscape_Configurational_het"), method="REML",level=95)
summary (LandscapeConfig)

########################################################################################################################################################################################

# the effects of crop heterogeneity (i.e. crop compositional and configurational heterogeneity components) and landscape heterogeneity (i.e. landscape compositional and configurational heterogeneity components) on plant diversity
table(dat$Heterogeneity_Type_L_C)
mod4 <- rma.mv(yi, vi, mods = ~ factor(Heterogeneity_Type_L_C) - 1, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Verts_Inverts_Plants_Without_Pests=="Plants"), method="REML",level=95)
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

# intraclass correlation (ICC) of the true effects
# compute the ICC (i.e. estimated correlation of the effects belonging to the same study, or how strong the dependency of true effect coming from the same study, or how strong the correlation of within study  heterogeneity)
round(mod4$sigma2[1] / sum(mod4$sigma2), 4)

# total amount of heterogeneity 
round(sum(mod4$sigma2), 4)

# forest plot 

#forest(mod1, addpred=TRUE, header=TRUE, atransf=transf.ztor, annotate=TRUE, order="obs", xlim=c(-1.5 , 1.5), addfit=TRUE, slab=dat$Author_With_Scale)

#check for publication bias 

# standard funnel plot with the SEs on the y-axis
#funnel(mod4)

# contour-enhanced fuel plot
#funnel(mod4, xlim=c(-1.5,1.5), ylim=c(0,0.5), refline=0, level=c(90,95,99),
       #shade=c("gray90","gray65","gray55"), back="white", hlines=NULL, las=1, 
       #legend=TRUE, digits=list(1L,1))

# put 1/sqrt(ni) on the y-axis
#funnel(mod4, yaxis="sqrtninv", xlim=c(-1,1), ylim=c(0,0.5), las=1)

# I2 for multilevel and multivariate models
# https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/mod4$vi)
X <- model.matrix(mod4)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(mod4$sigma2) / (sum(mod4$sigma2) + (mod4$k-mod4$p)/sum(diag(P)))

# we can also break things down to estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately

100 * mod4$sigma2 / (sum(mod4$sigma2) + (mod4$k-mod4$p)/sum(diag(P)))

# test difference between the four outcomes (crop and landscape heterogeneity)
anova(mod4, L=c(1, -1),digits = 4)

# graphical display 
t4 <- summary(rob4, atransf=transf.ztor)
plot_summs(t4,
           inner_ci_level=0.90,
           colors = "pink")

########################################################################################################################################################################################
# sub-group models
# this is just to count how many studies and effect sizes included in each test
table(dat$Heterogeneity_Type_L_C)
CropHet <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Verts_Inverts_Plants_Without_Pests=="Plants")&(Heterogeneity_Type_L_C=="Crop_het"), method="REML",level=95)
summary (CropHet)

LandsHet <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, subset=(Biodiversity_Category=="SWDiversity")&(Verts_Inverts_Plants_Without_Pests=="Plants")&(Heterogeneity_Type_L_C=="Landscape_het"), method="REML",level=95)
summary (LandsHet)

########################################################################################################################################################################################

# graphical display of all models 
#https://jtools.jacob-long.com/reference/plot_summs.html
plot_summs(t1,
           t2,
           t4,
           t3,
           inner_ci_level=0.90,
           point.shape = FALSE,
           point.size = 12.5,
           colors = "CUD Bright")

########################################################################################################################################################################################

