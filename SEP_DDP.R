# Codes accompanying "Separate Exchangeability as 
# Modeling Principle in Bayesian Nonparametrics"

# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.15
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(DirichletReg)
library(invgamma)
library(tidyr)
library(dplyr)
library(truncnorm)
library(scales)
library(viridisLite)
library(ggplot2)
theme_set(theme_bw(base_size = 14))
library(RColorBrewer)
library(salso)
library(ggpubr)
library(salso)
library(reshape2)
library(mvtnorm)
library(invgamma)
library(splines)
rm(list=ls())
set.seed(135)

load("Data-and-Results/data_protein.RData")

# cleaning data and construct design matrix 
rownames(PL)=c()
PL=cbind(c(rep(0,21),rep(1,21)),PL)
colnames(PL)[1]="z"

PL2=as.data.frame(PL) %>% group_by(z,age) %>% 
  summarise_each(funs(mean(., na.rm = TRUE))) %>% as.matrix()
PL2=PL2[,-which(sapply(1:dim(PL2)[2], function(i) all(is.na(PL[,i]))))]

age=PL2[,2]
age_std=(PL2[,2]-mean(PL2[,2]))/sd(PL2[,2])
x=cbind(PL2[,1],age_std,PL2[,1]*age_std)
y= PL2[,3:dim(PL2)[2]]

#total number of proteins
n=dim(y)[2]
str(y)
str(x); x

# construct splines using bs package
bsx2=bs(x[,2],degree = 3,knots = 
          c(quantile(x[,2],0.33), quantile(x[,2],0.67)),intercept = T)
x=cbind(bsx2,x[,1]*bsx2)

