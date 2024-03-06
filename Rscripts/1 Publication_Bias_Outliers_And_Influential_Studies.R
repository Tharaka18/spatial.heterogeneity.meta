########################################################################################################################################################################################

# Article title: Crop and landscape heterogeneity increase biodiversity in agricultural landscapes: A global review and meta-analysis

# Journal: Ecology Letters

# Authors: Tharaka S. Priyadarshana (ORCID: 0000-0003-3962-5465), Emily A. Martin (0000-0001-5785-9105), Clélia Sirami (0000-0003-1741-3082), Ben A. Woodcock (0000-0003-0300-9951), Eben Goodale (0000-0003-3403-2847), Carlos Martínez-Núñez (0000-0001-7814-4985), Myung-Bok Lee (0000-0003-2680-5707), Emilio Pagani-Núñez (0000-0001-8839-4005), Chloé A. Raderschall (0000-0003-2005-1705), Lluís Brotons (0000-0002-4826-4457), Anushka Rege (0000-0002-8383-0258), Annie Ouin (0000-0001-7046-2719), Teja Tscharntke (0000-0002-4482-3178), Eleanor M. Slade (0000-0002-6108-1196)

# contact: Tharaka S. Priyadarshana; tharakas001@e.ntu.edu.sg; tharakas.priyadarshana@gmail.com

# test for publication bias and influential studies 

dat<- read.csv ("Data/dat.csv", header=TRUE, sep=",",na.strings=TRUE)

# load metafor package
library(metafor)
library(snow)

# calculate the effect size (r to z)
dat <- escalc(measure="ZCOR", ri=Pearson_Correlation, ni=Sample_Size, data=dat)
head(dat)

mod1 <- rma.mv(yi, vi, random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, method="REML",level=95)
summary(mod1)

# back-transform the results of correlations and obtain the prediction interval

summary(mod1, atransf=transf.ztor)

# use cluster-robust inference methods
rob1 <- robust(mod1, cluster=dat$Study_ID)
summary(rob1, atransf=transf.ztor)

########################################################################################################################################################################################

# check for publication bias using funnel plot

# standard funnel plot with the SEs on the y-axis

funnel(mod1)

# contour-enhanced fuel plot
funnel(mod1, xlim=c(-1.5,1.5), ylim=c(0,0.5), refline=0, level=c(90,95,99),
       shade=c("gray95","gray75","gray65"), back="white", hlines=NULL, las=1, 
       legend=TRUE, digits=list(1L,1), xlab = "Fisher's z transformed correlation coefficient", ylab = "Standard error")

# put 1/sqrt(ni) on the y-axis
funnel(mod1, yaxis="sqrtninv", xlim=c(-1.5,1.5), ylim=c(0,0.5), las=1,back="white")

########################################################################################################################################################################################

# test for publication bias statistically 
# run the regression test (with SEs as predictor) manually 
# add the standard errors to the dataset 
dat$sei <- sqrt(dat$vi)
mod1pb <- rma.mv(yi, vi, mods = ~ sei , random = ~ 1 | Study_ID / Effect_Size_ID, data=dat, method="REML",level=95)
summary(mod1pb)

########################################################################################################################################################################################

# check for outliers
(x <- cooks.distance(mod1, parallel="multicore", scientific=FALSE)) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x, type="o", pch=19, xlab="Observed outcome", ylab="Cook's distance", las=1) #yes all below 1 

########################################################################################################################################################################################

# test for influential studies 

# frist let's aggregate the effect sizes to random mixed-effect models 
# https://wviechtb.github.io/metafor/reference/aggregate.escalc.html?q=aggregate#MathJax-script

agg <- aggregate(dat, cluster = dat$Study_ID, struct="ID", addk=TRUE)

res1 <- rma(yi, vi, method="DL", level=95, data=agg)
summary(res1)

plot(res1, qqplot=TRUE) # normality is also fine 

# now examine the funnel plot at study level 
funnel(res1, las=1)

# contour-enhanced fuel plot
funnel(res1, xlim=c(-2,2), ylim=c(0,0.6), refline=0, level=c(90,95,99),
       shade=c("gray90","gray75","gray65"), back="white", hlines=NULL, las=1, 
       legend=TRUE, digits=list(1L,1))

########################################################################################################################################################################################

# baujat plots to check for influential studies 
baujat(res1, ylab="Influence on overall result", xlab="Contribution to overall heterogeneity", ylim=c(0,0.06), xlim=c(0,2.3),las=1, grid = FALSE, bty="l")

# there is no influential studies. Study 66 has the most influence on overall results, but it is still just above 0.05. To be sure, let's see removing study 66 and creat "GOSH" plot 

########################################################################################################################################################################################

# check for outliers
# https://www.metafor-project.org/doku.php/plots:gosh_plot
## GOSH plots (note: these take a while to run! might want to decrease 'subsets' value)
res.GOSH <- gosh(res1, subsets=1000000, parallel = "snow", ncpus=40) 
res.GOSH
plot(res.GOSH, out=66, cex=1, breaks=c(100,50))

########################################################################################################################################################################################

# check for outliers (suing the aggregated model)
(x1 <- cooks.distance(res1, parallel="multicore", scientific=FALSE)) #this may take some time, depend on the power of the computer. To run this line faster use parallel="multicore"
plot(x1, type="o", pch=19, xlab="Observed outcome", ylab="Cook's distance", las=1, bty="l") #yes all below 1 

########################################################################################################################################################################################
