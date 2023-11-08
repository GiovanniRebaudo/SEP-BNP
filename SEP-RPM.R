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
library(RColorBrewer)
library(salso)


###################################################N
## create global variables for data & field dimensions
# Read and clean the data
###################################################N
## read in data
read.dta  <- function()
{##read-in and clean data
  Dietswap_dataset <- readRDS("Data-and-Results/Dietswap_dataset.RDS")
  
  data = Dietswap_dataset[Dietswap_dataset$timepoint==1,]
  data2=pivot_wider(data[,1:3],names_from = Sample,values_from = Abundance)
  na_in_rows=sapply(1:dim(data2)[1],function(i) sum(data2[i,]==0))
  data3=data.frame(data2[na_in_rows!=38,]) #remove all na rows
  
  
  Y0=data3[,2:39]                ## separate OTU counts
  rownames(Y0)=data3[,1]         ## and OTU names as row names..
  
  n <- dim(Y0)[2]                ## n patients
  drop = apply(Y0==0,1,sum)>n/4  ## drop all OTUS with more than n/4 0's
  Y  <- Y0[!drop,]               ## dropping 20 of the originally 119 OTU's
  ## now total # 0's is only 5.
  gammas=colMeans(Y)
  Y=sapply(1:n,function(i) as.matrix(Y)[,i]/gammas[i])
  y=log(Y+0.05)                  ## log transform (we only have 5 0's)
  return(list(y=y,Y=Y)) ## return both - the log transf "y" and raw "Y"
}
## create global variables for data & fied dimensions
out <- read.dta() ## maybe as.matrix(..)?? See if needed..
y <-   out$y    ## log transf
Y  <-  out$Y    ## abs scale

N <- dim(y)[1]*dim(y)[2]  # global var's
n <- dim(y)[2] # # patients
B <- dim(y)[1] # # OTU's

##### max size of the subject (K) and OTU (L) partition
K=10
L=15
## hyperprior pars
alpha= 1      ## for GEM prior on pi
beta = 1      ## GEM prior on w[k], k=1...K

## hyperpars (use "empirical bayes" i.e., match sample moments)
mu0=mean(as.matrix(y))                  ## mu_l ~ N(mu0,sig0)
sig0=3  ## SD!

a0=3                                    ## sig2_l ~ Ga(a0/2, b0/2)
b0=var(c(y))/4

source("SEP_fcts.R")

ex()
